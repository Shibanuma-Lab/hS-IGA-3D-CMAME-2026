"""
NURBS utility functions
Direct translation from Mathematica code in function.txt
"""

import numpy as np


def find_span(n, p, u, U):
    """
    Find the knot span index
    
    Parameters:
    -----------
    n : int
        Number of control points - 1 (nPts - 1)
    p : int
        Degree
    u : float
        Parameter value
    U : array
        Knot vector
    
    Returns:
    --------
    span : int
        Knot span index
    """
    # Special case: u equals last knot
    if u >= U[n + 1]:
        return n
    
    # Binary search
    low = p
    high = n + 1
    mid = (low + high) // 2
    
    while u < U[mid] or u >= U[mid + 1]:
        if u < U[mid]:
            high = mid
        else:
            low = mid
        mid = (low + high) // 2
    
    return mid


def basis_funs(span, u, p, U):
    """
    Compute non-zero basis functions
    
    Parameters:
    -----------
    span : int
        Knot span index from find_span
    u : float
        Parameter value
    p : int
        Degree
    U : array
        Knot vector
    
    Returns:
    --------
    N : array
        Array of (p+1) basis function values
    """
    N = np.zeros(p + 1)
    left = np.zeros(p + 1)
    right = np.zeros(p + 1)
    
    N[0] = 1.0
    
    for j in range(1, p + 1):
        left[j] = u - U[span + 1 - j]
        right[j] = U[span + j] - u
        saved = 0.0
        
        for r in range(j):
            temp = N[r] / (right[r + 1] + left[j - r])
            N[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        
        N[j] = saved
    
    return N


def der_basis_funs(span, u, p, order, U):
    """
    Compute basis functions and derivatives
    
    Parameters:
    -----------
    span : int
        Knot span index
    u : float
        Parameter value
    p : int
        Degree
    order : int
        Derivative order (0 or 1)
    U : array
        Knot vector
    
    Returns:
    --------
    ders : array
        Array of shape (order+1, p+1) containing basis functions and derivatives
        ders[0, :] = basis functions
        ders[1, :] = first derivatives
    """
    # 1. Compute basis functions
    nMat = np.zeros((p + 1, p + 1))
    left = np.zeros(p + 1)
    right = np.zeros(p + 1)
    
    nMat[0, 0] = 1.0
    
    for j in range(1, p + 1):
        left[j] = u - U[span + 1 - j]
        right[j] = U[span + j] - u
        saved = 0.0
        
        for r in range(j):
            # Lower triangle
            nMat[j, r] = right[r + 1] + left[j - r]
            temp = nMat[r, j - 1] / nMat[j, r]
            
            # Upper triangle
            nMat[r, j] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        
        nMat[j, j] = saved
    
    # 2. Load zero-th derivatives (basis values)
    ders = np.zeros((2, p + 1))
    for j in range(p + 1):
        ders[0, j] = nMat[j, p]
    
    if order == 0:
        return ders
    
    # 3. Compute derivatives
    aMat = np.zeros((2, p + 1))
    
    for j in range(p + 1):
        aMat[0, j] = nMat[j, p]
    
    # Compute first derivative
    for k in range(1, order + 1):
        d = 0.0
        rk = k
        pk = p - k
        
        if rk <= p:
            aMat[1, 0] = aMat[0, 0] / nMat[pk + 1, rk]
            d = aMat[1, 0] * nMat[rk, pk]
        
        j1 = 1 if rk > 0 else rk + 1
        j2 = pk if rk < p else pk + 1
        
        for j in range(j1, j2 + 1):
            aMat[1, j] = (aMat[0, j] - aMat[0, j - 1]) / nMat[pk + 1, rk + j]
            d += aMat[1, j] * nMat[rk + j, pk]
        
        if rk < p:
            aMat[1, pk + 1] = -aMat[0, pk] / nMat[pk + 1, pk + 1]
            d += aMat[1, pk + 1] * nMat[pk + 1, pk]
        
        ders[k, :] = d * np.math.factorial(p) / np.math.factorial(p - k)
        
        # Swap rows for next iteration
        for j in range(p + 1):
            aMat[0, j] = aMat[1, j]
    
    return ders


def nurbs_3d_basis_ders(spanU, spanV, spanW, p, q, r, 
                         knotU, knotV, knotW, 
                         xi, eta, zeta, 
                         weightsGlobal, nU, nV, nW):
    """
    Compute NURBS basis functions and derivatives in 3D
    
    Parameters:
    -----------
    spanU, spanV, spanW : int
        Knot span indices
    p, q, r : int
        Degrees in each direction
    knotU, knotV, knotW : array
        Knot vectors
    xi, eta, zeta : float
        Parameter space coordinates
    weightsGlobal : array
        Global weights vector
    nU, nV, nW : int
        Number of control points - 1 in each direction
    
    Returns:
    --------
    R : array
        NURBS basis functions
    dRdxi, dRdeta, dRdzeta : array
        Derivatives of NURBS basis functions
    """
    # 1. Extract local weights
    uind = spanU - p
    vind = spanV - q
    wind = spanW - r
    
    # Build local indices for this element
    localIdx = []
    for k in range(r + 1):
        wk = spanW - r + k
        for j in range(q + 1):
            vj = spanV - q + j
            for i in range(p + 1):
                ui = spanU - p + i
                # Global index (0-based)
                globalIdx = wk * (nU + 1) * (nV + 1) + vj * (nU + 1) + ui
                localIdx.append(globalIdx)
    
    weightLocal = weightsGlobal[localIdx]
    
    # 2. Compute B-spline basis functions and derivatives
    Nx = basis_funs(spanU, xi, p, knotU)
    Ny = basis_funs(spanV, eta, q, knotV)
    Nz = basis_funs(spanW, zeta, r, knotW)
    
    DNx = der_basis_funs(spanU, xi, p, 1, knotU)[1, :]
    DNy = der_basis_funs(spanV, eta, q, 1, knotV)[1, :]
    DNz = der_basis_funs(spanW, zeta, r, 1, knotW)[1, :]
    
    # 3. Flatten in (k,j,i) order
    basis = []
    derivX = []
    derivY = []
    derivZ = []
    
    for k in range(r + 1):
        for j in range(q + 1):
            for i in range(p + 1):
                basis.append(Nx[i] * Ny[j] * Nz[k])
                derivX.append(DNx[i] * Ny[j] * Nz[k])
                derivY.append(Nx[i] * DNy[j] * Nz[k])
                derivZ.append(Nx[i] * Ny[j] * DNz[k])
    
    basis = np.array(basis)
    derivX = np.array(derivX)
    derivY = np.array(derivY)
    derivZ = np.array(derivZ)
    
    # 4. Weight and normalize
    num = basis * weightLocal
    wTot = np.sum(num)
    
    dwdxi = np.sum(derivX * weightLocal)
    dwdeta = np.sum(derivY * weightLocal)
    dwdzeta = np.sum(derivZ * weightLocal)
    
    R = num / wTot
    dRdxi = (derivX * weightLocal * wTot - num * dwdxi) / (wTot**2)
    dRdeta = (derivY * weightLocal * wTot - num * dwdeta) / (wTot**2)
    dRdzeta = (derivZ * weightLocal * wTot - num * dwdzeta) / (wTot**2)
    
    return R, dRdxi, dRdeta, dRdzeta


def parent_to_parametric_space(xi_range, xi_bar):
    """
    Map from parent element [-1, 1] to parametric space [xi_min, xi_max]
    
    Parameters:
    -----------
    xi_range : tuple
        (xi_min, xi_max)
    xi_bar : float
        Coordinate in parent space
    
    Returns:
    --------
    xi : float
        Coordinate in parametric space
    """
    xi_min, xi_max = xi_range
    xi = 0.5 * ((xi_max - xi_min) * xi_bar + (xi_max + xi_min))
    return xi
