#!/usr/bin/env python3
"""
Calculate L2 norm for S-IGA verification results
Translated from Mathematica code with optimizations

This script computes the relative L2 norm between numerical and analytical (Sneddon) solutions:
    relative_L2_norm = sqrt(∫|u_num - u_ana|²dV / ∫|u_ana|²dV)
"""

import numpy as np
import scipy.io as sio
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay, cKDTree
from scipy.optimize import fsolve, minimize
import json
import os
from pathlib import Path
from datetime import datetime
import argparse
from multiprocessing import Pool, cpu_count
import time
from numba import jit

# Global variables for worker processes
_worker_interpolators = None
_worker_bg_mesh = None
_worker_local_region = None


# ============================================================================
# Numba-accelerated numerical functions for interpolation
# ============================================================================

@jit(nopython=True, cache=True, fastmath=True)
def shape_functions_numba(xi, eta, zeta):
    """
    Numba-compiled H8 hexahedral shape functions.
    About 2-3x faster than pure NumPy version.
    
    Args:
        xi, eta, zeta: Isoparametric coordinates [-1, 1]
    
    Returns:
        N: (8,) array of shape function values
    """
    N = np.empty(8, dtype=np.float64)
    N[0] = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta)
    N[1] = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta)
    N[2] = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta)
    N[3] = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta)
    N[4] = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta)
    N[5] = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta)
    N[6] = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta)
    N[7] = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta)
    return N


@jit(nopython=True, cache=True, fastmath=True)
def shape_derivatives_numba(xi, eta, zeta):
    """
    Numba-compiled derivatives of H8 shape functions.
    Returns [dN/dxi, dN/deta, dN/dzeta] as (3, 8) array.
    
    Args:
        xi, eta, zeta: Isoparametric coordinates [-1, 1]
    
    Returns:
        dN: (3, 8) array of shape function derivatives
    """
    dN = np.empty((3, 8), dtype=np.float64)
    
    # dN/dxi for each node
    dN[0, 0] = -0.125 * (1 - eta) * (1 - zeta)
    dN[0, 1] =  0.125 * (1 - eta) * (1 - zeta)
    dN[0, 2] =  0.125 * (1 + eta) * (1 - zeta)
    dN[0, 3] = -0.125 * (1 + eta) * (1 - zeta)
    dN[0, 4] = -0.125 * (1 - eta) * (1 + zeta)
    dN[0, 5] =  0.125 * (1 - eta) * (1 + zeta)
    dN[0, 6] =  0.125 * (1 + eta) * (1 + zeta)
    dN[0, 7] = -0.125 * (1 + eta) * (1 + zeta)
    
    # dN/deta for each node
    dN[1, 0] = -0.125 * (1 - xi) * (1 - zeta)
    dN[1, 1] = -0.125 * (1 + xi) * (1 - zeta)
    dN[1, 2] =  0.125 * (1 + xi) * (1 - zeta)
    dN[1, 3] =  0.125 * (1 - xi) * (1 - zeta)
    dN[1, 4] = -0.125 * (1 - xi) * (1 + zeta)
    dN[1, 5] = -0.125 * (1 + xi) * (1 + zeta)
    dN[1, 6] =  0.125 * (1 + xi) * (1 + zeta)
    dN[1, 7] =  0.125 * (1 - xi) * (1 + zeta)
    
    # dN/dzeta for each node
    dN[2, 0] = -0.125 * (1 - xi) * (1 - eta)
    dN[2, 1] = -0.125 * (1 + xi) * (1 - eta)
    dN[2, 2] = -0.125 * (1 + xi) * (1 + eta)
    dN[2, 3] = -0.125 * (1 - xi) * (1 + eta)
    dN[2, 4] =  0.125 * (1 - xi) * (1 - eta)
    dN[2, 5] =  0.125 * (1 + xi) * (1 - eta)
    dN[2, 6] =  0.125 * (1 + xi) * (1 + eta)
    dN[2, 7] =  0.125 * (1 - xi) * (1 + eta)
    
    return dN


@jit(nopython=True, cache=True, fastmath=True)
def inverse_isoparametric_numba(point, elem_coords, tol, initial_guess):
    """
    Numba-compiled Newton-Raphson solver for inverse isoparametric mapping.
    This is the performance-critical function - about 60-70% of total time.
    
    Solves: point = sum_i N_i(xi, eta, zeta) * node_i
    for (xi, eta, zeta) using Newton-Raphson iteration.
    
    Args:
        point: (3,) array - physical coordinates to find
        elem_coords: (8, 3) array - element node coordinates
        tol: Convergence tolerance
        initial_guess: (3,) array - starting point [xi, eta, zeta] or None
    
    Returns:
        result: (4,) array - [xi, eta, zeta, success_flag]
                success_flag: 1.0 if converged, 0.0 otherwise
    """
    # Initial guess
    if initial_guess is not None and len(initial_guess) == 3:
        xi = initial_guess[0]
        eta = initial_guess[1]
        zeta = initial_guess[2]
    else:
        xi = 0.0
        eta = 0.0
        zeta = 0.0
    
    max_iter = 50  # Increased from 20 for better accuracy
    for iteration in range(max_iter):
        # Compute shape functions and derivatives
        N = shape_functions_numba(xi, eta, zeta)
        dN = shape_derivatives_numba(xi, eta, zeta)
        
        # Current physical position: x = N · coords
        x_current = np.zeros(3, dtype=np.float64)
        for i in range(8):
            for j in range(3):
                x_current[j] += N[i] * elem_coords[i, j]
        
        # Jacobian matrix: J = dN · coords (3x3)
        J = np.zeros((3, 3), dtype=np.float64)
        for i in range(3):  # dN row
            for j in range(3):  # spatial dimension
                for k in range(8):  # node
                    J[i, j] += dN[i, k] * elem_coords[k, j]
        
        # Residual: r = point - x_current
        residual = point - x_current
        
        # Check convergence
        residual_norm = np.sqrt(residual[0]**2 + residual[1]**2 + residual[2]**2)
        if residual_norm < tol:
            return np.array([xi, eta, zeta, 1.0], dtype=np.float64)
        
        # Newton update: solve J^T · delta = residual
        # Using Gaussian elimination with partial pivoting
        JT = J.T.copy()
        b = residual.copy()
        
        # Forward elimination with partial pivoting
        for i in range(3):
            # Find pivot
            max_val = abs(JT[i, i])
            max_row = i
            for k in range(i + 1, 3):
                if abs(JT[k, i]) > max_val:
                    max_val = abs(JT[k, i])
                    max_row = k
            
            # Check for singular matrix
            if max_val < 1e-14:
                return np.array([xi, eta, zeta, 0.0], dtype=np.float64)
            
            # Swap rows if needed
            if max_row != i:
                for j in range(3):
                    JT[i, j], JT[max_row, j] = JT[max_row, j], JT[i, j]
                b[i], b[max_row] = b[max_row], b[i]
            
            # Eliminate
            for k in range(i + 1, 3):
                factor = JT[k, i] / JT[i, i]
                for j in range(i, 3):
                    JT[k, j] -= factor * JT[i, j]
                b[k] -= factor * b[i]
        
        # Back substitution
        delta = np.zeros(3, dtype=np.float64)
        for i in range(2, -1, -1):
            delta[i] = b[i]
            for j in range(i + 1, 3):
                delta[i] -= JT[i, j] * delta[j]
            delta[i] /= JT[i, i]
        
        # Update
        xi += delta[0]
        eta += delta[1]
        zeta += delta[2]
    
    # Check final convergence
    N = shape_functions_numba(xi, eta, zeta)
    x_final = np.zeros(3, dtype=np.float64)
    for i in range(8):
        for j in range(3):
            x_final[j] += N[i] * elem_coords[i, j]
    
    residual = point - x_final
    residual_norm = np.sqrt(residual[0]**2 + residual[1]**2 + residual[2]**2)
    
    # Relaxed final check: allow up to 100x tolerance (was 10x)
    if residual_norm < tol * 100:
        return np.array([xi, eta, zeta, 1.0], dtype=np.float64)
    else:
        return np.array([xi, eta, zeta, 0.0], dtype=np.float64)


def _init_worker(ugx, ugy, ugz, uglx, ugly, uglz, bg_mesh, sneddon):
    """Initialize worker process with interpolators"""
    global _worker_interpolators, _worker_bg_mesh, _worker_local_region
    _worker_interpolators = {
        'ugx': ugx, 'ugy': ugy, 'ugz': ugz,
        'uglx': uglx, 'ugly': ugly, 'uglz': uglz,
        'sneddon': sneddon
    }
    _worker_bg_mesh = bg_mesh
    # Local region parameters
    _worker_local_region = (0.751, 1.249, 0.0, 0.25)


def _worker_process_batch(args):
    """Worker function that uses global interpolators"""
    elem_indices, quad_points_3d, quad_weights_3d = args
    
    global _worker_interpolators, _worker_bg_mesh, _worker_local_region
    
    ugx = _worker_interpolators['ugx']
    ugy = _worker_interpolators['ugy']
    ugz = _worker_interpolators['ugz']
    uglx = _worker_interpolators['uglx']
    ugly = _worker_interpolators['ugly']
    uglz = _worker_interpolators['uglz']
    sneddon = _worker_interpolators['sneddon']
    
    r_min, r_max, z_min, z_max = _worker_local_region
    
    batch_error = 0.0
    batch_exact = 0.0
    
    for ie in elem_indices:
        elem_nodes = _worker_bg_mesh.elements[ie]
        node_coords = _worker_bg_mesh.nodes[elem_nodes]
        
        # Process all quadrature points for this element
        for (xi, eta, zeta), w in zip(quad_points_3d, quad_weights_3d):
            # Evaluate shape functions
            N = HexahedronElement.shape_functions(xi, eta, zeta)
            dN = HexahedronElement.shape_derivatives(xi, eta, zeta)
            
            # Compute Jacobian
            J = HexahedronElement.jacobian(dN, node_coords)
            det_J = np.abs(np.linalg.det(J))  # Use absolute value for integration
            
            # Physical coordinates of quadrature point
            xyz = N @ node_coords
            x, y, z = xyz
            
            # Get numerical displacement (inline to avoid method call overhead)
            r = np.sqrt(x**2 + y**2)
            if (r_min <= r <= r_max) and (z_min <= z <= z_max):
                ux = uglx(x, y, z)
                uy = ugly(x, y, z)
                uz = uglz(x, y, z)
            else:
                ux = ugx(x, y, z)
                uy = ugy(x, y, z)
                uz = ugz(x, y, z)
            
            u_num = np.array([ux, uy, uz])
            
            # Get analytical displacement
            u_ana = np.array(sneddon.evaluate_cartesian(x, y, z))
            
            # Compute error
            error = u_num - u_ana
            
            # Add to integrals
            batch_error += w * np.dot(error, error) * det_J
            batch_exact += w * np.dot(u_ana, u_ana) * det_J
    
    return batch_error, batch_exact


class ElementMeshInterpolator:
    """
    Mimics Mathematica's ElementMeshInterpolation function.
    
    Uses actual finite element mesh (nodes + elements) for interpolation,
    preserving the mesh topology rather than creating a Delaunay triangulation.
    
    This is equivalent to Mathematica's:
        mesh = ToElementMesh[nodes]
        interp = ElementMeshInterpolation[{mesh}, values]
    """
    def __init__(self, nodes, elements, values, name="FEM", fill_value=0.0):
        """
        Args:
            nodes: (N, 3) array of nodal coordinates
            elements: (M, 8) array of hexahedral element connectivity (0-based)
            values: (N,) array of nodal values
            name: Name for debugging
            fill_value: Value to return for points outside mesh (default: 0.0)
        """
        self.nodes = np.array(nodes)
        self.elements = np.array(elements, dtype=int)
        self.values = np.array(values)
        self.name = name
        self.fill_value = fill_value
        
        # Cache for last successful element search (spatial locality optimization)
        self._last_elem_cache = None
        self._last_xi_eta_zeta_cache = None
        self._cache_hits = 0
        self._cache_misses = 0
        
        print(f"  Building ElementMeshInterpolator for {name}...")
        print(f"    {len(self.nodes)} nodes, {len(self.elements)} elements")
        print(f"    Using Numba JIT-compiled Newton-Raphson solver")
        
        # Pre-compute element bounding boxes for fast search
        self._build_element_bbox()
        
        # Build KD-tree for O(log n) spatial search
        self._build_kdtree()
    
    def _build_element_bbox(self):
        """Pre-compute axis-aligned bounding boxes for all elements"""
        n_elem = len(self.elements)
        self.bbox_min = np.zeros((n_elem, 3))
        self.bbox_max = np.zeros((n_elem, 3))
        
        for i, elem_nodes_idx in enumerate(self.elements):
            elem_coords = self.nodes[elem_nodes_idx]
            self.bbox_min[i] = np.min(elem_coords, axis=0)
            self.bbox_max[i] = np.max(elem_coords, axis=0)
    
    def _build_kdtree(self):
        """Build KD-tree for fast spatial queries using element centroids"""
        n_elem = len(self.elements)
        
        # Compute element centroids
        self.elem_centroids = np.zeros((n_elem, 3))
        for i, elem_nodes_idx in enumerate(self.elements):
            elem_coords = self.nodes[elem_nodes_idx]
            self.elem_centroids[i] = np.mean(elem_coords, axis=0)
        
        # Build KD-tree for O(log n) nearest neighbor search
        self.kdtree = cKDTree(self.elem_centroids)
        print(f"    KD-tree built for {n_elem} elements")
    
    def __call__(self, x, y, z):
        """Interpolate at point (x, y, z)"""
        point = np.array([x, y, z])
        
        # Find element containing the point
        elem_idx, xi, eta, zeta = self._find_containing_element(point)
        
        if elem_idx is None:
            return self.fill_value
        
        # Interpolate using H8 shape functions
        return self._interpolate_in_element(elem_idx, xi, eta, zeta)
    
    def _find_containing_element(self, point, tol=1e-10):  # Increased from 1e-6
        """
        Find which element contains the point.
        Returns: (elem_idx, xi, eta, zeta) or (None, None, None, None)
        """
        x, y, z = point
        
        # Step 1: Check cache first (spatial locality optimization)
        # Adjacent quadrature points often fall in the same element
        if self._last_elem_cache is not None:
            elem_coords = self.nodes[self.elements[self._last_elem_cache]]
            result = self._inverse_isoparametric(point, elem_coords, tol, 
                                                 initial_guess=self._last_xi_eta_zeta_cache)
            
            if result is not None:
                xi, eta, zeta = result
                if (abs(xi) <= 1.0 + tol and 
                    abs(eta) <= 1.0 + tol and 
                    abs(zeta) <= 1.0 + tol):
                    # Cache hit! Update cache and return
                    self._cache_hits += 1
                    self._last_xi_eta_zeta_cache = (xi, eta, zeta)
                    return self._last_elem_cache, xi, eta, zeta
        
        # Step 2: Cache miss - use KD-tree for O(log n) search
        self._cache_misses += 1
        
        # Query KD-tree for nearest elements (typically only need to check 5-10)
        # Use larger k for safety, but check in distance order
        k_neighbors = min(20, len(self.elements))  # Check at most 20 nearest elements
        distances, candidates = self.kdtree.query(point, k=k_neighbors)
        
        # Ensure candidates is iterable even if k=1
        if k_neighbors == 1:
            candidates = [candidates]
            distances = [distances]
        
        # Check candidates in order of increasing distance from centroid
        for elem_idx in candidates:
            elem_coords = self.nodes[self.elements[elem_idx]]
            
            # Try to find isoparametric coordinates
            result = self._inverse_isoparametric(point, elem_coords, tol)
            
            if result is not None:
                xi, eta, zeta = result
                # Check if point is inside element (with tolerance)
                if (abs(xi) <= 1.0 + tol and 
                    abs(eta) <= 1.0 + tol and 
                    abs(zeta) <= 1.0 + tol):
                    # Update cache with successful result
                    self._last_elem_cache = elem_idx
                    self._last_xi_eta_zeta_cache = (xi, eta, zeta)
                    return elem_idx, xi, eta, zeta
        
        # Not found - clear cache
        self._last_elem_cache = None
        self._last_xi_eta_zeta_cache = None
        return None, None, None, None
    
    def _inverse_isoparametric(self, point, elem_coords, tol=1e-10, initial_guess=None):  # Increased from 1e-6
        """
        Solve inverse isoparametric mapping: find (xi, eta, zeta) such that
        point = sum_i N_i(xi, eta, zeta) * node_i
        
        Uses Numba-compiled Newton-Raphson iteration for 2-3x speedup.
        
        Args:
            point: Physical coordinates to find
            elem_coords: Element node coordinates
            tol: Convergence tolerance
            initial_guess: Optional (xi, eta, zeta) to start iteration (improves convergence)
        """
        # Call Numba-compiled version for performance
        result = inverse_isoparametric_numba(
            np.asarray(point, dtype=np.float64),
            np.asarray(elem_coords, dtype=np.float64),
            tol,
            np.asarray(initial_guess, dtype=np.float64) if initial_guess is not None else None
        )
        
        # result = [xi, eta, zeta, success_flag]
        if result[3] > 0.5:  # success
            return result[:3]
        else:
            return None
    
    def _shape_functions(self, xi, eta, zeta):
        """H8 hexahedral shape functions"""
        return 0.125 * np.array([
            (1 - xi) * (1 - eta) * (1 - zeta),  # Node 0
            (1 + xi) * (1 - eta) * (1 - zeta),  # Node 1
            (1 + xi) * (1 + eta) * (1 - zeta),  # Node 2
            (1 - xi) * (1 + eta) * (1 - zeta),  # Node 3
            (1 - xi) * (1 - eta) * (1 + zeta),  # Node 4
            (1 + xi) * (1 - eta) * (1 + zeta),  # Node 5
            (1 + xi) * (1 + eta) * (1 + zeta),  # Node 6
            (1 - xi) * (1 + eta) * (1 + zeta),  # Node 7
        ])
    
    def _shape_derivatives(self, xi, eta, zeta):
        """Derivatives of H8 shape functions: [dN/dxi, dN/deta, dN/dzeta]"""
        return 0.125 * np.array([
            # dN/dxi for each node
            [-(1 - eta) * (1 - zeta), (1 - eta) * (1 - zeta),
             (1 + eta) * (1 - zeta), -(1 + eta) * (1 - zeta),
             -(1 - eta) * (1 + zeta), (1 - eta) * (1 + zeta),
             (1 + eta) * (1 + zeta), -(1 + eta) * (1 + zeta)],
            # dN/deta for each node
            [-(1 - xi) * (1 - zeta), -(1 + xi) * (1 - zeta),
             (1 + xi) * (1 - zeta), (1 - xi) * (1 - zeta),
             -(1 - xi) * (1 + zeta), -(1 + xi) * (1 + zeta),
             (1 + xi) * (1 + zeta), (1 - xi) * (1 + zeta)],
            # dN/dzeta for each node
            [-(1 - xi) * (1 - eta), -(1 + xi) * (1 - eta),
             -(1 + xi) * (1 + eta), -(1 - xi) * (1 + eta),
             (1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
             (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)],
        ])
    
    def _interpolate_in_element(self, elem_idx, xi, eta, zeta):
        """Interpolate value at isoparametric coordinates (xi, eta, zeta)"""
        # Get shape functions
        N = self._shape_functions(xi, eta, zeta)
        
        # Get nodal values for this element
        elem_nodes = self.elements[elem_idx]
        elem_values = self.values[elem_nodes]
        
        # Interpolate
        return float(N @ elem_values)


class BilinearQuadInterpolator:
    """
    Bilinear interpolator on quadrilateral mesh elements for 2D (r,z) Sneddon data.
    Matches Mathematica's Interpolation[..., InterpolationOrder->1] behavior.
    """
    def __init__(self, nodes_rz, values, name="2D"):
        """
        Args:
            nodes_rz: (N, 2) array of (r, z) coordinates
            values: (N,) array of values to interpolate
            name: Name for debugging
        """
        self.nodes = np.array(nodes_rz)
        self.values = np.array(values)
        self.name = name
        
        # Build structured grid if possible
        self._build_structured_grid()
    
    def _build_structured_grid(self):
        """Try to detect structured grid for fast interpolation"""
        # Extract unique r and z values
        r_vals = np.unique(self.nodes[:, 0])
        z_vals = np.unique(self.nodes[:, 1])
        
        # Check if it's a structured grid
        if len(r_vals) * len(z_vals) == len(self.nodes):
            # Likely structured - build fast lookup
            self.r_grid = r_vals
            self.z_grid = z_vals
            self.nr = len(r_vals)
            self.nz = len(z_vals)
            
            # Reshape values into 2D grid
            self.grid_values = np.full((self.nz, self.nr), np.nan)
            for i, (r, z) in enumerate(self.nodes):
                ir = np.searchsorted(r_vals, r)
                iz = np.searchsorted(z_vals, z)
                self.grid_values[iz, ir] = self.values[i]
            
            self.is_structured = True
        else:
            self.is_structured = False
    
    def __call__(self, r, z):
        """Interpolate at point (r, z)"""
        if self.is_structured:
            return self._interpolate_structured(r, z)
        else:
            return self._interpolate_unstructured(r, z)
    
    def _interpolate_structured(self, r, z):
        """Fast interpolation on structured grid"""
        # Find cell containing point
        if r < self.r_grid[0] or r > self.r_grid[-1] or \
           z < self.z_grid[0] or z > self.z_grid[-1]:
            return 0.0  # Outside domain
        
        # Find indices
        ir = np.searchsorted(self.r_grid, r) - 1
        iz = np.searchsorted(self.z_grid, z) - 1
        
        # Clamp to valid range
        ir = max(0, min(ir, self.nr - 2))
        iz = max(0, min(iz, self.nz - 2))
        
        # Bilinear interpolation
        r0, r1 = self.r_grid[ir], self.r_grid[ir + 1]
        z0, z1 = self.z_grid[iz], self.z_grid[iz + 1]
        
        # Normalized coordinates
        xi = (r - r0) / (r1 - r0) if r1 > r0 else 0.0
        eta = (z - z0) / (z1 - z0) if z1 > z0 else 0.0
        
        # Get corner values
        v00 = self.grid_values[iz, ir]
        v10 = self.grid_values[iz, ir + 1]
        v11 = self.grid_values[iz + 1, ir + 1]
        v01 = self.grid_values[iz + 1, ir]
        
        # Bilinear formula
        value = (1 - xi) * (1 - eta) * v00 + \
                xi * (1 - eta) * v10 + \
                xi * eta * v11 + \
                (1 - xi) * eta * v01
        
        return float(value)
    
    def _interpolate_unstructured(self, r, z):
        """Fallback for unstructured data"""
        # Use simple nearest neighbor for now
        dists = np.sqrt((self.nodes[:, 0] - r)**2 + (self.nodes[:, 1] - z)**2)
        idx = np.argmin(dists)
        if dists[idx] < 0.1:  # Within reasonable distance
            return float(self.values[idx])
        return 0.0


class TrilinearHexInterpolator:
    """
    Trilinear interpolator on hexahedral mesh elements for 3D displacement data.
    Matches Mathematica's Interpolation[..., InterpolationOrder->1] for 3D.
    """
    def __init__(self, nodes, elements, values, name="3D"):
        """
        Args:
            nodes: (N, 3) array of (x, y, z) coordinates
            elements: (M, 8) array of element connectivity (0-based)
            values: (N,) array of nodal values
            name: Name for debugging
        """
        self.nodes = np.array(nodes)
        self.values = np.array(values)
        self.name = name
        
        # Ensure elements are 0-based
        self.elements = np.array(elements, dtype=int)
        if len(elements) > 0 and np.min(self.elements) == 1:
            self.elements = self.elements - 1
        
        # Build bounding boxes for elements
        self._build_element_bbox()
    
    def _build_element_bbox(self):
        """Pre-compute bounding boxes for all elements"""
        nelem = len(self.elements)
        self.bbox_min = np.zeros((nelem, 3))
        self.bbox_max = np.zeros((nelem, 3))
        
        for e in range(nelem):
            elem_nodes = self.nodes[self.elements[e]]
            self.bbox_min[e] = np.min(elem_nodes, axis=0)
            self.bbox_max[e] = np.max(elem_nodes, axis=0)
    
    def _find_containing_element(self, x, y, z, tol=1e-6):
        """Find which element contains point (x, y, z)"""
        bbox_tol = tol * 10
        
        # Quick filtering using bounding boxes
        candidates = np.where(
            (x >= self.bbox_min[:, 0] - bbox_tol) & (x <= self.bbox_max[:, 0] + bbox_tol) &
            (y >= self.bbox_min[:, 1] - bbox_tol) & (y <= self.bbox_max[:, 1] + bbox_tol) &
            (z >= self.bbox_min[:, 2] - bbox_tol) & (z <= self.bbox_max[:, 2] + bbox_tol)
        )[0]
        
        # Check each candidate
        for e in candidates:
            elem_coords = self.nodes[self.elements[e]]
            xi, eta, zeta, success = self._point_to_parent(x, y, z, elem_coords, tol)
            if success:
                return e, xi, eta, zeta
        
        return None, None, None, None
    
    def _point_to_parent(self, x, y, z, elem_coords, tol=1e-6):
        """
        Convert physical point (x,y,z) to parent coordinates (xi, eta, zeta).
        Returns: (xi, eta, zeta, success)
        """
        # Check if element is axis-aligned box (fast path)
        x_coords = elem_coords[:, 0]
        y_coords = elem_coords[:, 1]
        z_coords = elem_coords[:, 2]
        
        rect_tol = tol * 100
        
        # Check if it's an axis-aligned hexahedron
        x_vals = np.unique(np.round(x_coords / rect_tol) * rect_tol)
        y_vals = np.unique(np.round(y_coords / rect_tol) * rect_tol)
        z_vals = np.unique(np.round(z_coords / rect_tol) * rect_tol)
        
        if len(x_vals) == 2 and len(y_vals) == 2 and len(z_vals) == 2:
            # Axis-aligned box
            x_min, x_max = np.min(x_coords), np.max(x_coords)
            y_min, y_max = np.min(y_coords), np.max(y_coords)
            z_min, z_max = np.min(z_coords), np.max(z_coords)
            
            x_center = 0.5 * (x_min + x_max)
            y_center = 0.5 * (y_min + y_max)
            z_center = 0.5 * (z_min + z_max)
            
            x_half = 0.5 * (x_max - x_min)
            y_half = 0.5 * (y_max - y_min)
            z_half = 0.5 * (z_max - z_min)
            
            if x_half > 1e-12 and y_half > 1e-12 and z_half > 1e-12:
                xi = (x - x_center) / x_half
                eta = (y - y_center) / y_half
                zeta = (z - z_center) / z_half
                
                if abs(xi) <= 1.0 + tol*10 and abs(eta) <= 1.0 + tol*10 and abs(zeta) <= 1.0 + tol*10:
                    return xi, eta, zeta, True
        
        # General case: Newton's method (skip for performance)
        return None, None, None, False
    
    def __call__(self, x, y, z):
        """Interpolate at point (x, y, z)"""
        e, xi, eta, zeta = self._find_containing_element(x, y, z)
        
        if e is None:
            return 0.0  # Outside all elements
        
        # Trilinear interpolation using H8 shape functions
        N = np.array([
            0.125 * (1 - xi) * (1 - eta) * (1 - zeta),
            0.125 * (1 + xi) * (1 - eta) * (1 - zeta),
            0.125 * (1 + xi) * (1 + eta) * (1 - zeta),
            0.125 * (1 - xi) * (1 + eta) * (1 - zeta),
            0.125 * (1 - xi) * (1 - eta) * (1 + zeta),
            0.125 * (1 + xi) * (1 - eta) * (1 + zeta),
            0.125 * (1 + xi) * (1 + eta) * (1 + zeta),
            0.125 * (1 - xi) * (1 + eta) * (1 + zeta)
        ])
        
        elem_values = self.values[self.elements[e]]
        return float(N @ elem_values)


class SneddonSolution:
    """Sneddon analytical solution for penny-shaped crack"""
    
    def __init__(self, mat_file='sneddon_python.mat', p0=1.0, c=1.0, EE=100.0, nu=0.3):
        """
        Initialize Sneddon solution
        
        Args:
            mat_file: Path to Sneddon data file
            p0: Applied stress
            c: Crack radius
            EE: Young's modulus
            nu: Poisson's ratio
        """
        self.p0 = p0
        self.c = c
        self.EE = EE
        self.nu = nu
        
        # Load Sneddon data
        print(f"Loading Sneddon solution from: {mat_file}")
        mat_data = sio.loadmat(mat_file)
        
        # Extract position and solution arrays
        # posA is already (N, 2) array of [r, z] pairs
        self.pos_2d = mat_data['posA']  # Shape: (N, 2)
        self.SA = mat_data['SA']  # Solution array: [ur1, ur2, uz1, uz2], Shape: (N, 4)
        
        print(f"  pos shape: {self.pos_2d.shape}")
        print(f"  SA shape: {self.SA.shape}")
        
        # Create bilinear interpolators (matching Mathematica's InterpolationOrder->1)
        print("  Building BilinearQuadInterpolator for Sneddon solution...")
        self.ur1_interp = BilinearQuadInterpolator(self.pos_2d, self.SA[:, 0], name="ur1")
        self.ur2_interp = BilinearQuadInterpolator(self.pos_2d, self.SA[:, 1], name="ur2")
        self.uz1_interp = BilinearQuadInterpolator(self.pos_2d, self.SA[:, 2], name="uz1")
        self.uz2_interp = BilinearQuadInterpolator(self.pos_2d, self.SA[:, 3], name="uz2")
    
    def evaluate(self, r, z):
        """
        Evaluate Sneddon solution at (r, z) in cylindrical coordinates
        
        Args:
            r: Radial coordinate
            z: Axial coordinate
            
        Returns:
            (ur, uz): Radial and axial displacements
        """
        # Far-field solution
        ur0 = -(self.nu * self.p0 / self.EE) * r
        uz0 = (self.p0 / self.EE) * z
        
        # Interpolate Sneddon correction terms using bilinear interpolation
        ur1 = self.ur1_interp(r, z)
        ur2 = self.ur2_interp(r, z)
        uz1 = self.uz1_interp(r, z)
        uz2 = self.uz2_interp(r, z)
        
        # Total displacement
        if r == 0:
            ur = 0.0
        else:
            ur = (2 * self.p0 * self.c * (1 + self.nu) / (np.pi * self.EE)) * \
                 ((1 - 2*self.nu) * ur1 - ur2)
        
        if r >= self.c and z == 0:
            uz = 0.0
        else:
            uz = -(4 * self.p0 * self.c * (1 - self.nu**2) / (np.pi * self.EE)) * \
                 (uz1 + uz2 / (2 * (1 - self.nu)))
        
        return ur0 + ur, uz0 + uz
    
    def evaluate_cartesian(self, x, y, z):
        """
        Evaluate Sneddon solution at (x, y, z) in Cartesian coordinates
        
        Args:
            x, y, z: Cartesian coordinates
            
        Returns:
            (ux, uy, uz): Displacements in Cartesian coordinates
        """
        r = np.sqrt(x**2 + y**2)
        
        if r == 0:
            theta = np.pi / 2
        else:
            theta = np.arctan2(y, x)
        
        ur, uz = self.evaluate(r, z)
        
        # Convert to Cartesian
        ux = ur * np.cos(theta)
        uy = ur * np.sin(theta)
        
        return ux, uy, uz


class BackgroundMesh:
    """Background mesh for integration"""
    
    def __init__(self, hB, max_x, max_y, max_z):
        """
        Create background mesh for integration
        
        Args:
            hB: Background mesh size
            max_x, max_y, max_z: Domain dimensions
        """
        self.hB = hB
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z
        
        # Number of elements in each direction
        self.nGxy = int(np.ceil(max_x / hB))
        self.nGz = int(np.ceil(max_z / hB))
        
        # Adjust mesh size for exact fit
        self.hGxy = max_x / self.nGxy
        self.hGz = max_z / self.nGz
        
        # Create nodal coordinates
        nodeGxy = np.linspace(0, max_x, self.nGxy + 1)
        nodeGz = np.linspace(0, max_z, self.nGz + 1)
        
        # Create 3D mesh nodes (z, y, x order as in Mathematica)
        self.nodes = []
        for iz in range(len(nodeGz)):
            for iy in range(len(nodeGxy)):
                for ix in range(len(nodeGxy)):
                    x = nodeGxy[ix]
                    y = nodeGxy[iy]
                    z = nodeGz[iz]
                    # Adjust small values to zero
                    if abs(x) < 1e-9 * max_x:
                        x = 0.0
                    if abs(y) < 1e-9 * max_x:
                        y = 0.0
                    if abs(z) < 1e-9 * max_x:
                        z = 0.0
                    self.nodes.append([x, y, z])
        
        self.nodes = np.array(self.nodes)
        
        # Create element connectivity
        self.elements = self._create_elements()
        
        print(f"Background mesh created:")
        print(f"  hB = {hB:.6f}, hGxy = {self.hGxy:.6f}, hGz = {self.hGz:.6f}")
        print(f"  nGxy = {self.nGxy}, nGz = {self.nGz}")
        print(f"  Total nodes: {len(self.nodes)}")
        print(f"  Total elements: {len(self.elements)}")
    
    def _create_elements(self):
        """Create hexahedral element connectivity"""
        nGxy = self.nGxy
        nGz = self.nGz
        
        elements = []
        for iz in range(nGz):
            for iy in range(nGxy):
                for ix in range(nGxy):
                    # Node indices (1-based in Mathematica, 0-based here)
                    n1 = iz * (nGxy + 1)**2 + iy * (nGxy + 1) + ix
                    n2 = iz * (nGxy + 1)**2 + iy * (nGxy + 1) + (ix + 1)
                    n3 = iz * (nGxy + 1)**2 + (iy + 1) * (nGxy + 1) + (ix + 1)
                    n4 = iz * (nGxy + 1)**2 + (iy + 1) * (nGxy + 1) + ix
                    n5 = (iz + 1) * (nGxy + 1)**2 + iy * (nGxy + 1) + ix
                    n6 = (iz + 1) * (nGxy + 1)**2 + iy * (nGxy + 1) + (ix + 1)
                    n7 = (iz + 1) * (nGxy + 1)**2 + (iy + 1) * (nGxy + 1) + (ix + 1)
                    n8 = (iz + 1) * (nGxy + 1)**2 + (iy + 1) * (nGxy + 1) + ix
                    
                    elements.append([n1, n2, n3, n4, n5, n6, n7, n8])
        
        return np.array(elements)


class GaussQuadrature:
    """Gauss quadrature for integration"""
    
    def __init__(self, order=4):
        """
        Initialize Gauss quadrature
        
        Args:
            order: Quadrature order (number of points per dimension)
        """
        self.order = order
        
        # Gauss points and weights - matching Mathematica intpco and wlist
        # These are taken from the Mathematica code's intpco[[order]] and wlist[[order]]
        gauss_data = {
            2: {
                'points': [-0.5773502691896258, 0.5773502691896258],
                'weights': [1.0, 1.0]
            },
            3: {
                'points': [-0.7745966692414834, 0.0, 0.7745966692414834],
                'weights': [0.5555555555555556, 0.8888888888888888, 0.5555555555555556]
            },
            4: {
                'points': [-0.8302961484013275, -0.40957436820775056, 
                          0.40957436820775056, 0.8302961484013275],
                'weights': [0.3478548451374538, 0.6521451548625461,
                           0.6521451548625461, 0.3478548451374538]
            }
        }
        
        if order not in gauss_data:
            raise ValueError(f"Order {order} not supported. Use 2 or 4.")
        
        self.points_1d = np.array(gauss_data[order]['points'])
        self.weights_1d = np.array(gauss_data[order]['weights'])
        
        # Create 3D quadrature points and weights
        self.points_3d = []
        self.weights_3d = []
        
        for i in range(order):
            for j in range(order):
                for k in range(order):
                    xi = self.points_1d[k]
                    eta = self.points_1d[j]
                    zeta = self.points_1d[i]
                    w = self.weights_1d[k] * self.weights_1d[j] * self.weights_1d[i]
                    self.points_3d.append([xi, eta, zeta])
                    self.weights_3d.append(w)
        
        self.points_3d = np.array(self.points_3d)
        self.weights_3d = np.array(self.weights_3d)


class HexahedronElement:
    """Hexahedral element shape functions"""
    
    @staticmethod
    def shape_functions(xi, eta, zeta):
        """
        Evaluate shape functions at (xi, eta, zeta)
        
        Returns:
            N: Shape function values (8,)
        """
        N = 0.125 * np.array([
            (1 - xi) * (1 - eta) * (1 - zeta),
            (1 + xi) * (1 - eta) * (1 - zeta),
            (1 + xi) * (1 + eta) * (1 - zeta),
            (1 - xi) * (1 + eta) * (1 - zeta),
            (1 - xi) * (1 - eta) * (1 + zeta),
            (1 + xi) * (1 - eta) * (1 + zeta),
            (1 + xi) * (1 + eta) * (1 + zeta),
            (1 - xi) * (1 + eta) * (1 + zeta)
        ])
        return N
    
    @staticmethod
    def shape_derivatives(xi, eta, zeta):
        """
        Evaluate shape function derivatives at (xi, eta, zeta)
        
        Returns:
            dN: Shape function derivatives (3, 8)
        """
        dN = 0.125 * np.array([
            # dN/dxi
            [-(1 - eta) * (1 - zeta), (1 - eta) * (1 - zeta),
             (1 + eta) * (1 - zeta), -(1 + eta) * (1 - zeta),
             -(1 - eta) * (1 + zeta), (1 - eta) * (1 + zeta),
             (1 + eta) * (1 + zeta), -(1 + eta) * (1 + zeta)],
            # dN/deta
            [-(1 - xi) * (1 - zeta), -(1 + xi) * (1 - zeta),
             (1 + xi) * (1 - zeta), (1 - xi) * (1 - zeta),
             -(1 - xi) * (1 + zeta), -(1 + xi) * (1 + zeta),
             (1 + xi) * (1 + zeta), (1 - xi) * (1 + zeta)],
            # dN/dzeta
            [-(1 - xi) * (1 - eta), -(1 + xi) * (1 - eta),
             -(1 + xi) * (1 + eta), -(1 - xi) * (1 + eta),
             (1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
             (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)]
        ])
        return dN
    
    @staticmethod
    def jacobian(dN, node_coords):
        """
        Compute Jacobian matrix
        
        Args:
            dN: Shape function derivatives (3, 8)
            node_coords: Element node coordinates (8, 3)
            
        Returns:
            J: Jacobian matrix (3, 3)
        """
        J = dN @ node_coords
        return J


class L2NormCalculator:
    """Calculate L2 norm between numerical and analytical solutions"""
    
    def __init__(self, result_folder, sneddon_file='sneddon_python.mat',
                 p0=1.0, c=1.0, EE=100.0, nu=0.3):
        """
        Initialize L2 norm calculator
        
        Args:
            result_folder: Path to result folder (str or Path)
            sneddon_file: Path to Sneddon data file
            p0, c, EE, nu: Material parameters
        """
        self.result_folder = Path(result_folder) if isinstance(result_folder, str) else result_folder
        
        # Load configuration
        config_file = self.result_folder / "run_config.json"
        with open(config_file, 'r') as f:
            self.config = json.load(f)
        
        self.hL = self.config['hL']
        
        # Initialize Sneddon solution
        self.sneddon = SneddonSolution(sneddon_file, p0, c, EE, nu)
        
        # Find step directory
        step_dirs = list(self.result_folder.glob("step*"))
        if not step_dirs:
            raise ValueError(f"No step directory found in {self.result_folder}")
        self.step_dir = step_dirs[0]
        
        print(f"\nProcessing: {self.result_folder.name}")
        print(f"  hL = {self.hL:.8f}")
        print(f"  Step directory: {self.step_dir.name}")
        
        # Load mesh and displacement data
        self._load_data()
        
        # Create background mesh
        self._create_background_mesh()
        
        # Create interpolators
        self._create_interpolators()
    
    def _load_data(self):
        """Load node coordinates and displacement data"""
        # Load global nodes and displacements
        node_g_file = self.step_dir / "node.g.dat"
        u_g_file = self.step_dir / "log" / "u.g.dat"
        
        self.node_g = np.loadtxt(node_g_file, skiprows=1, usecols=(1, 2, 3))
        self.u_g = np.loadtxt(u_g_file, skiprows=1, usecols=(1, 2, 3))
        
        # Reconstruct global mesh elements from structured grid
        # Global nodes are arranged in (x, y, z) order
        self.elem_g = self._reconstruct_structured_elements(self.node_g)
        
        # Load local nodes and displacements
        node_l_file = self.step_dir / "node.l.dat"
        u_gl_l_file = self.step_dir / "log" / "u_gl.l.dat"
        elem_l_file = self.step_dir / "elem.l.dat"
        
        self.node_l = np.loadtxt(node_l_file, skiprows=1, usecols=(1, 2, 3))
        self.u_gl_l = np.loadtxt(u_gl_l_file, skiprows=1, usecols=(1, 2, 3))
        
        # Load local elements (already H8 format)
        elem_l_full = np.loadtxt(elem_l_file, skiprows=1, dtype=int)
        self.elem_l = elem_l_full[:, 1:9] - 1  # Skip elem ID, convert to 0-based
        
        print(f"  Global nodes: {len(self.node_g)}, elements: {len(self.elem_g)}")
        print(f"  Local nodes: {len(self.node_l)}, elements: {len(self.elem_l)}")
        print(f"  Total DOF: {(len(self.node_g) + len(self.node_l)) * 3}")
    
    def _reconstruct_structured_elements(self, nodes):
        """
        Reconstruct H8 elements from structured node array.
        Nodes are arranged in (x, y, z) lexicographic order.
        
        Returns:
            elements: (num_elem, 8) array of element connectivity (0-based)
        """
        # Extract unique coordinates in each direction
        x_coords = np.unique(np.round(nodes[:, 0], 10))
        y_coords = np.unique(np.round(nodes[:, 1], 10))
        z_coords = np.unique(np.round(nodes[:, 2], 10))
        
        nx, ny, nz = len(x_coords), len(y_coords), len(z_coords)
        
        print(f"  Detected structured grid: nx={nx}, ny={ny}, nz={nz}")
        
        # Verify node count
        expected_nodes = nx * ny * nz
        if expected_nodes != len(nodes):
            print(f"  WARNING: Expected {expected_nodes} nodes, got {len(nodes)}")
        
        # Generate H8 elements
        # Node ordering: first x, then y, then z
        # Node index: i + j*nx + k*nx*ny
        elements = []
        
        for k in range(nz - 1):
            for j in range(ny - 1):
                for i in range(nx - 1):
                    # 8 corner nodes of hexahedron (H8 ordering)
                    n0 = i + j*nx + k*nx*ny
                    n1 = (i+1) + j*nx + k*nx*ny
                    n2 = (i+1) + (j+1)*nx + k*nx*ny
                    n3 = i + (j+1)*nx + k*nx*ny
                    n4 = i + j*nx + (k+1)*nx*ny
                    n5 = (i+1) + j*nx + (k+1)*nx*ny
                    n6 = (i+1) + (j+1)*nx + (k+1)*nx*ny
                    n7 = i + (j+1)*nx + (k+1)*nx*ny
                    
                    elements.append([n0, n1, n2, n3, n4, n5, n6, n7])
        
        return np.array(elements, dtype=int)
    
    def _create_background_mesh(self):
        """Create background mesh for integration"""
        max_x = np.max(self.node_g[:, 0])
        max_y = np.max(self.node_g[:, 1])
        max_z = np.max(self.node_g[:, 2])
        
        hB = self.hL * 1.5
        
        self.bg_mesh = BackgroundMesh(hB, max_x, max_y, max_z)
    
    def _create_interpolators(self):
        """Create displacement interpolators using ElementMeshInterpolation"""
        print("Creating displacement interpolators...")
        print("  Using ElementMeshInterpolator (FEM-based) - mimics Mathematica's ElementMeshInterpolation")
        
        # Global displacement interpolators using actual FEM mesh
        self.ugx = ElementMeshInterpolator(self.node_g, self.elem_g, self.u_g[:, 0], name="Global-X")
        self.ugy = ElementMeshInterpolator(self.node_g, self.elem_g, self.u_g[:, 1], name="Global-Y")
        self.ugz = ElementMeshInterpolator(self.node_g, self.elem_g, self.u_g[:, 2], name="Global-Z")
        
        # Local displacement interpolators using actual FEM mesh
        self.uglx = ElementMeshInterpolator(self.node_l, self.elem_l, self.u_gl_l[:, 0], name="Local-X")
        self.ugly = ElementMeshInterpolator(self.node_l, self.elem_l, self.u_gl_l[:, 1], name="Local-Y")
        self.uglz = ElementMeshInterpolator(self.node_l, self.elem_l, self.u_gl_l[:, 2], name="Local-Z")
    
    def _is_in_local_region(self, x, y, z):
        """
        Check if point is in local region
        
        Local region: 0.751 <= r <= 1.249, 0 <= z <= 0.25
        """
        r = np.sqrt(x**2 + y**2)
        return (0.751 <= r <= 1.249) and (0.0 <= z <= 0.25)
    
    def _get_numerical_displacement(self, x, y, z):
        """Get numerical displacement at (x, y, z)"""
        if self._is_in_local_region(x, y, z):
            # Use local interpolator
            ux = self.uglx(x, y, z)
            uy = self.ugly(x, y, z)
            uz = self.uglz(x, y, z)
        else:
            # Use global interpolator
            ux = self.ugx(x, y, z)
            uy = self.ugy(x, y, z)
            uz = self.ugz(x, y, z)
        
        return np.array([ux, uy, uz])
    
    def _process_element_batch(self, args):
        """
        Process a batch of elements for parallel computation.
        
        Args:
            args: Tuple of (element_indices, quad_points_3d, quad_weights_3d)
            
        Returns:
            Tuple of (integral_error, integral_exact) for this batch
        """
        elem_indices, quad_points_3d, quad_weights_3d = args
        
        batch_error = 0.0
        batch_exact = 0.0
        
        for ie in elem_indices:
            elem_nodes = self.bg_mesh.elements[ie]
            node_coords = self.bg_mesh.nodes[elem_nodes]
            
            # Process all quadrature points for this element
            for (xi, eta, zeta), w in zip(quad_points_3d, quad_weights_3d):
                # Evaluate shape functions
                N = HexahedronElement.shape_functions(xi, eta, zeta)
                dN = HexahedronElement.shape_derivatives(xi, eta, zeta)
                
                # Compute Jacobian
                J = HexahedronElement.jacobian(dN, node_coords)
                det_J = np.abs(np.linalg.det(J))  # Use absolute value for integration
                
                # Physical coordinates of quadrature point
                xyz = N @ node_coords
                x, y, z = xyz
                
                # Get numerical displacement
                u_num = self._get_numerical_displacement(x, y, z)
                
                # Get analytical displacement
                u_ana = np.array(self.sneddon.evaluate_cartesian(x, y, z))
                
                # Compute error
                error = u_num - u_ana
                
                # Add to integrals
                batch_error += w * np.dot(error, error) * det_J
                batch_exact += w * np.dot(u_ana, u_ana) * det_J
        
        return batch_error, batch_exact
    
    def calculate(self, quadrature_order=4, n_processes=None):
        """
        Calculate L2 norm using Gauss quadrature with parallel processing.
        
        Args:
            quadrature_order: Quadrature order (2 or 4)
            n_processes: Number of parallel processes (default: CPU count)
            
        Returns:
            Dictionary with results
        """
        if n_processes is None:
            n_processes = max(1, cpu_count() - 1)  # Leave one CPU free
        
        print(f"\nCalculating L2 norm with {quadrature_order}-point Gauss quadrature...")
        print(f"  Using {n_processes} parallel processes")
        
        start_time = time.time()
        
        # Initialize quadrature
        quad = GaussQuadrature(quadrature_order)
        
        # Adaptive batch sizing strategy
        n_elements = len(self.bg_mesh.elements)
        
        # Step 1: Shuffle element indices to balance load
        # Background elements are spatially ordered - elements near local region are MUCH slower
        # Random shuffling ensures each batch has a mix of fast/slow elements
        print(f"  Shuffling {n_elements} elements for load balancing...")
        element_indices = np.arange(n_elements)
        np.random.seed(42)  # Reproducible shuffling
        np.random.shuffle(element_indices)
        
        # Step 2: Measure processing time for a small sample
        print(f"  Calibrating batch size...")
        calibration_size = min(200, n_elements // 20)
        calibration_batch = (element_indices[:calibration_size].tolist(), quad.points_3d, quad.weights_3d)
        
        cal_start = time.time()
        _ = self._process_element_batch(calibration_batch)
        cal_time = time.time() - cal_start
        
        # Estimate time per element
        time_per_element = cal_time / calibration_size if calibration_size > 0 else 0.001
        
        # Step 3: Determine optimal batch size
        # Target: each batch takes 2-5 seconds to balance computation vs communication
        # Target: each process gets 8-15 batches for good load balancing
        target_batch_time = 3.0  # seconds
        optimal_batch_size = max(50, int(target_batch_time / time_per_element))
        
        # Ensure each process gets enough batches (at least 8)
        min_batches_per_process = 8
        max_batch_size = max(50, n_elements // (n_processes * min_batches_per_process))
        batch_size = min(optimal_batch_size, max_batch_size)
        
        # Ensure batch_size is reasonable
        batch_size = max(50, min(batch_size, n_elements // 4))
        
        num_batches = int(np.ceil(n_elements / batch_size))
        batches_per_process = num_batches / n_processes
        estimated_batch_time = batch_size * time_per_element
        
        print(f"  Calibration: {calibration_size} elements in {cal_time:.2f}s ({time_per_element*1000:.2f}ms/elem)")
        print(f"  Adaptive batch_size: {batch_size} elements/batch")
        print(f"  Total batches: {num_batches} (~{batches_per_process:.1f} per process)")
        print(f"  Estimated batch time: {estimated_batch_time:.2f}s")
        
        # Step 4: Create batches from shuffled indices
        batches = []
        for i in range(0, n_elements, batch_size):
            batch_elem_indices = element_indices[i:min(i + batch_size, n_elements)].tolist()
            batches.append((batch_elem_indices, quad.points_3d, quad.weights_3d))
        
        print(f"  Processing {n_elements} elements in {len(batches)} batches...")
        
        # Process batches in parallel with pre-initialized workers
        if n_processes > 1:
            # Initialize worker processes with interpolators (done once, not per batch!)
            with Pool(processes=n_processes, 
                     initializer=_init_worker,
                     initargs=(self.ugx, self.ugy, self.ugz,
                              self.uglx, self.ugly, self.uglz,
                              self.bg_mesh, self.sneddon)) as pool:
                # Use global worker function instead of instance method
                results = list(pool.imap_unordered(_worker_process_batch, batches, chunksize=1))
                print(f"  Progress: 100.0% ({len(batches)}/{len(batches)} batches)")
        else:
            # Serial processing fallback
            results = []
            for i, batch in enumerate(batches):
                results.append(self._process_element_batch(batch))
                if (i + 1) % max(1, len(batches) // 20) == 0 or (i + 1) == len(batches):
                    progress = 100 * (i + 1) / len(batches)
                    elapsed = time.time() - start_time
                    rate = (i + 1) / elapsed if elapsed > 0 else 0
                    eta = (len(batches) - i - 1) / rate if rate > 0 else 0
                    print(f"  Progress: {progress:.1f}% ({i+1}/{len(batches)} batches, {elapsed:.1f}s, ETA: {eta:.0f}s)")
        
        # Sum up results from all batches
        integral_error = sum(r[0] for r in results)
        integral_exact = sum(r[1] for r in results)
        
        # Compute relative L2 norm
        relative_L2 = np.sqrt(integral_error / integral_exact)
        
        elapsed_time = time.time() - start_time
        
        print(f"\n  ✓ Integration complete in {elapsed_time:.1f}s")
        print(f"  Total computation time: {elapsed_time:.2f}s")
        print(f"  ∫|u_num - u_ana|² dV = {integral_error:.6e}")
        print(f"  ∫|u_ana|² dV = {integral_exact:.6e}")
        print(f"  Relative L2 norm = {relative_L2:.6e}")
        
        # Print cache statistics for local interpolators (most benefit from caching)
        print(f"\n  Cache performance (Local mesh):")
        for interp_name, interp in [('X', self.uglx), ('Y', self.ugly), ('Z', self.uglz)]:
            total_queries = interp._cache_hits + interp._cache_misses
            if total_queries > 0:
                hit_rate = 100.0 * interp._cache_hits / total_queries
                print(f"    {interp_name}: {hit_rate:.1f}% hit rate ({interp._cache_hits}/{total_queries} hits)")
        
        # Also show global mesh cache performance
        print(f"\n  Cache performance (Global mesh):")
        for interp_name, interp in [('X', self.ugx), ('Y', self.ugy), ('Z', self.ugz)]:
            total_queries = interp._cache_hits + interp._cache_misses
            if total_queries > 0:
                hit_rate = 100.0 * interp._cache_hits / total_queries
                print(f"    {interp_name}: {hit_rate:.1f}% hit rate ({interp._cache_hits}/{total_queries} hits)")
        
        return {
            'hL': self.hL,
            'hG': self.config['hG'],
            'rGL': self.config['rGL'],
            'dof': (len(self.node_g) + len(self.node_l)) * 3,
            'integral_error': integral_error,
            'integral_exact': integral_exact,
            'relative_L2_norm': relative_L2,
            'computation_time': elapsed_time
        }


def process_all_results(base_dir='results/verification_5_2', rGL=2,
                       sneddon_file='sneddon_python.mat',
                       output_file=None):
    """
    Process all result folders for a given rGL value
    
    Args:
        base_dir: Base results directory
        rGL: rGL value to process
        sneddon_file: Path to Sneddon data file
        output_file: Output CSV file (optional)
    """
    base_path = Path(base_dir)
    rGL_folder = base_path / f"rGL{rGL}_0.25"
    
    if not rGL_folder.exists():
        print(f"Error: Folder not found: {rGL_folder}")
        return
    
    # Find all result folders and sort by hL value (descending: coarse to fine mesh)
    result_folders = [d for d in rGL_folder.iterdir() if d.is_dir() and d.name.startswith('hL_')]
    
    # Extract hL value and sort from large to small (coarse to fine)
    def extract_hL(folder):
        try:
            # Extract hL value from folder name like "hL_0.009615_hG_0.019231"
            hL_str = folder.name.split('_')[1]
            return float(hL_str)
        except:
            return 0.0
    
    result_folders = sorted(result_folders, key=extract_hL, reverse=True)
    
    print(f"\n{'='*70}")
    print(f"Processing rGL = {rGL} (coarse to fine mesh)")
    print(f"{'='*70}")
    print(f"Found {len(result_folders)} result folders\n")
    
    results = []
    
    for i, folder in enumerate(result_folders):
        print(f"\n--- Processing {i+1}/{len(result_folders)}: {folder.name} ---")
        folder_start = time.time()
        
        try:
            calc = L2NormCalculator(folder, sneddon_file=sneddon_file)
            result = calc.calculate(quadrature_order=4)
            results.append(result)
            folder_time = time.time() - folder_start
            print(f"\n>>> Folder {folder.name} completed in {folder_time:.1f}s <<<")
        except Exception as e:
            print(f"Error processing {folder.name}: {e}")
            continue
    
    # Save results
    if output_file is None:
        output_file = base_path / f"L2_norm_rGL{rGL}.csv"
    
    # Write CSV
    import csv
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['rGL', 'hL', 'hG', 'dof', 'relative_L2_norm'])
        writer.writeheader()
        for r in results:
            writer.writerow({
                'rGL': r['rGL'],
                'hL': r['hL'],
                'hG': r['hG'],
                'dof': r['dof'],
                'relative_L2_norm': r['relative_L2_norm']
            })
    
    print(f"\n{'='*70}")
    print(f"Results saved to: {output_file}")
    print(f"{'='*70}")
    
    return results


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Calculate L2 norm for S-IGA verification results"
    )
    parser.add_argument(
        '--rGL',
        type=int,
        default=2,
        help='rGL value to process (default: 2)'
    )
    parser.add_argument(
        '--base-dir',
        type=str,
        default='results/verification_5_2',
        help='Base results directory (default: results/verification_5_2)'
    )
    parser.add_argument(
        '--sneddon-file',
        type=str,
        default='sneddon_python.mat',
        help='Sneddon data file (default: sneddon_python.mat)'
    )
    parser.add_argument(
        '--output',
        type=str,
        help='Output CSV file (default: auto-generated)'
    )
    
    args = parser.parse_args()
    
    process_all_results(
        base_dir=args.base_dir,
        rGL=args.rGL,
        sneddon_file=args.sneddon_file,
        output_file=args.output
    )


if __name__ == "__main__":
    main()
