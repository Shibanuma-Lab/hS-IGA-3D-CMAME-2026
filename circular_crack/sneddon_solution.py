"""
Sneddon analytical solution for penny-shaped crack under uniform tension
Based on the Mathematica SneddonApp.mx module

This module provides TWO implementations:
1. sneddon_displacement_elliptic() - Fast approximation using elliptic integrals
2. sneddon_displacement_interpolated() - Accurate Bessel integral + interpolation

For production use in static mode, use the interpolated version (like Mathematica).

Reference:
Sneddon, I. N. (1946). "The distribution of stress in the neighbourhood 
of a crack in an elastic solid". Proceedings of the Royal Society A.
"""

import numpy as np
from scipy import special
import os

# Global interpolator (loaded on first use)
_sneddon_interpolator = None
_sneddon_c = None


def _load_interpolator():
    """Load precomputed Sneddon interpolation data (lazy loading)"""
    global _sneddon_interpolator, _sneddon_c
    
    if _sneddon_interpolator is not None:
        return _sneddon_interpolator, _sneddon_c
    
    # Get the directory where this file is located
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Priority 1: Try Mathematica .mat file (most accurate)
    mat_file = os.path.join(current_dir, 'sneddon_SA.mat')
    if os.path.exists(mat_file):
        from sneddon_precompute import load_interpolation_data
        _sneddon_interpolator, _sneddon_c = load_interpolation_data(mat_file)
        return _sneddon_interpolator, _sneddon_c
    
    # Priority 2: Try Python .npz file
    data_file = os.path.join(current_dir, 'sneddon_interpolation.npz')
    if not os.path.exists(data_file):
        # Fall back to test data
        data_file = os.path.join(current_dir, 'sneddon_interpolation_test.npz')
        if not os.path.exists(data_file):
            raise FileNotFoundError(
                f"Sneddon interpolation data not found!\n"
                f"Please either:\n"
                f"  1. Export from Mathematica: sneddon_SA.mat\n"
                f"  2. Run: python sneddon_precompute.py\n"
                f"Expected locations:\n"
                f"  {current_dir}/sneddon_SA.mat (preferred)\n"
                f"  {current_dir}/sneddon_interpolation.npz"
            )
    
    from sneddon_precompute import load_interpolation_data
    _sneddon_interpolator, _sneddon_c = load_interpolation_data(data_file)
    
    return _sneddon_interpolator, _sneddon_c


def sneddon_displacement_interpolated(sigma_app, a, E, nu, point):
    """
    Calculate displacement using precomputed Bessel integrals + interpolation
    
    This is the ACCURATE method matching Mathematica's SneddonApp.
    Based on numerical integration of Bessel functions, precomputed on a grid.
    
    Parameters:
    -----------
    sigma_app : float
        Applied stress [Pa]
    a : float
        Crack radius [m]
    E : float
        Young's modulus [Pa]
    nu : float
        Poisson's ratio [-]
    point : array-like [r, z]
        Point coordinates in cylindrical system
        r: radial distance from crack center [m]
        z: distance from crack plane [m]
    
    Returns:
    --------
    [u_r, u_z] : array
        Radial and axial displacements [m]
    """
    r, z = point
    
    # Load interpolator
    interpolator, c_ref = _load_interpolator()
    
    # Background displacement (far-field uniform stress)
    u_r0 = -(nu * sigma_app / E) * r
    u_z0 = (sigma_app / E) * z
    
    # Get interpolated integral values [ur1, ur2, uz1, uz2]
    integrals = interpolator([r, z])[0]  # Returns shape (4,)
    ur1, ur2, uz1, uz2 = integrals
    
    # Calculate crack-induced displacement
    # Based on Mathematica SneddonApp:
    # ur = (2*p0*c)/(π*E) * ((1-2ν)*ur1 - ur2)
    # uz = -(4*p0*c*(1-ν²))/(π*E) * (uz1 + uz2/(2*(1-ν)))
    #
    # IMPORTANT: Both formulas contain (1+ν) factor:
    #   - ur has explicit (1+ν) factor (see note below)
    #   - uz has (1-ν²) = (1-ν)(1+ν), so (1+ν) is already included
    #
    # NOTE: The (1+ν) factor in ur is required to match Mathematica's SneddonApp output
    # when converting from cylindrical to Cartesian coordinates.
    
    if abs(r) < 1e-15:
        u_r_crack = 0.0
    else:
        u_r_crack = (2 * sigma_app * a * (1 + nu)) / (np.pi * E) * ((1 - 2*nu) * ur1 - ur2)
    
    if r >= a and abs(z) < 1e-15:
        u_z_crack = 0.0
    else:
        u_z_crack = -(4 * sigma_app * a * (1 - nu**2)) / (np.pi * E) * \
                    (uz1 + uz2 / (2 * (1 - nu)))
    
    # Total displacement = background + crack-induced
    u_r = u_r0 + u_r_crack
    u_z = u_z0 + u_z_crack
    
    return np.array([u_r, u_z])


def sneddon_displacement_elliptic(sigma_app, a, E, nu, point):
    """
    Calculate displacement using elliptic integral approximation (FAST but less accurate)
    
    This is a simplified formula using complete elliptic integrals.
    Faster than Bessel integrals but may have accuracy differences.
    
    Parameters:
    -----------
    sigma_app : float
        Applied stress [Pa]
    a : float
        Crack radius [m]
    E : float
        Young's modulus [Pa]
    nu : float
        Poisson's ratio [-]
    point : array-like [r, z]
        Point coordinates in cylindrical system
        r: radial distance from crack center [m]
        z: distance from crack plane [m]
    
    Returns:
    --------
    [u_r, u_z] : array
        Radial and axial displacements [m]
    """
    r, z = point
    
    # Material parameter
    # For plane strain: mu = E / (2 * (1 + nu))
    # Compliance: (1 - nu^2) / E for plane strain
    
    # Distance from point to crack center
    R = np.sqrt(r**2 + z**2)
    
    # Special cases
    if R < 1e-15:
        # At origin
        return np.array([0.0, 0.0])
    
    if abs(z) < 1e-15 and r < a:
        # On crack surface (inside crack)
        u_r = 0.0
        # Opening displacement (COD)
        u_z = (2 * sigma_app * (1 - nu**2) / E) * np.sqrt(a**2 - r**2)
        return np.array([u_r, u_z])
    
    # General case: outside crack
    # Sneddon's solution for penny-shaped crack
    
    # Auxiliary parameters
    lambda_ = np.sqrt((R - z)**2 / 4.0 + r**2)
    mu_ = np.sqrt((R + z)**2 / 4.0 + r**2)
    
    # Complete elliptic integrals
    # k^2 = 4*lambda*mu / (lambda + mu)^2
    # For scipy: ellipk(m) where m = k^2
    
    if lambda_ + mu_ < 1e-15:
        k_squared = 0.0
    else:
        k_squared = 4 * lambda_ * mu_ / (lambda_ + mu_)**2
    
    # Ensure k_squared is in valid range [0, 1]
    k_squared = np.clip(k_squared, 0.0, 1.0)
    
    # Complete elliptic integrals of first and second kind
    K = special.ellipk(k_squared)  # K(k^2)
    E_ellip = special.ellipe(k_squared)  # E(k^2)
    
    # Displacement components from Sneddon's solution
    # Prefactor
    C = (2 * sigma_app * a * (1 - nu**2)) / (np.pi * E)
    
    if r < 1e-15:
        # On z-axis (r = 0)
        u_r = 0.0
        # Axial displacement
        if abs(z) > a:
            # Outside crack region
            u_z = C * np.arcsin(a / R)
        else:
            # Very close to crack (shouldn't happen if z != 0)
            u_z = C * np.pi / 2.0
    else:
        # General point (r > 0)
        # Radial displacement
        u_r = C * (r / (lambda_ + mu_)) * ((2 - k_squared) * K - 2 * E_ellip)
        
        # Axial displacement
        u_z = C * (z / (lambda_ + mu_)) * ((2 - k_squared) * K - 2 * E_ellip)
    
    return np.array([u_r, u_z])


def sneddon_displacement(sigma_app, a, E, nu, point, use_interpolation=True):
    """
    Calculate displacement at a point due to penny-shaped crack
    
    Wrapper function that selects between interpolated (accurate) or elliptic (fast) method.
    
    Parameters:
    -----------
    sigma_app : float
        Applied stress [Pa]
    a : float
        Crack radius [m]
    E : float
        Young's modulus [Pa]
    nu : float
        Poisson's ratio [-]
    point : array-like [r, z]
        Point coordinates in cylindrical system
        r: radial distance from crack center [m]
        z: distance from crack plane [m]
    use_interpolation : bool
        If True, use accurate Bessel integral + interpolation (like Mathematica)
        If False, use fast elliptic integral approximation
    
    Returns:
    --------
    [u_r, u_z] : array
        Radial and axial displacements [m]
    """
    if use_interpolation:
        return sneddon_displacement_interpolated(sigma_app, a, E, nu, point)
    else:
        return sneddon_displacement_elliptic(sigma_app, a, E, nu, point)


def sneddon_displacement_cartesian(sigma_app, a, E, nu, point, use_interpolation=True):
    """
    Calculate displacement in Cartesian coordinates
    
    Parameters:
    -----------
    sigma_app : float
        Applied stress [Pa]
    a : float
        Crack radius [m]
    E : float
        Young's modulus [Pa]
    nu : float
        Poisson's ratio [-]
    point : array-like [x, y, z]
        Point coordinates in Cartesian system [m]
    use_interpolation : bool
        If True, use accurate Bessel integral + interpolation (like Mathematica)
        If False, use fast elliptic integral approximation
    
    Returns:
    --------
    [u_x, u_y, u_z] : array
        Displacements in x, y, z directions [m]
    """
    x, y, z = point
    
    # Convert to cylindrical coordinates
    r = np.sqrt(x**2 + y**2)
    
    # Calculate theta for coordinate transformation
    if abs(x) < 1e-15:
        theta = np.pi / 2.0
    else:
        theta = np.arctan(y / x)
    
    # Get displacements in cylindrical coordinates
    u_r, u_z = sneddon_displacement(sigma_app, a, E, nu, [r, z], use_interpolation)
    
    # Convert radial displacement to Cartesian using theta
    # This matches Mathematica's getbc: 
    # {{i, 1, disp[[1]]*Cos[θ]}, {i, 2, disp[[1]]*Sin[θ]}, {i, 3, disp[[2]]}}
    if r < 1e-15:
        # On z-axis
        u_x = 0.0
        u_y = 0.0
    else:
        # General point
        u_x = u_r * np.cos(theta)
        u_y = u_r * np.sin(theta)
    
    return np.array([u_x, u_y, u_z])


def get_boundary_displacement(node_coords, sigma_app, a, E, nu, use_interpolation=True):
    """
    Get boundary displacements for a list of nodes
    
    Parameters:
    -----------
    node_coords : array-like, shape (N, 3)
        Node coordinates [x, y, z] in meters
    sigma_app : float
        Applied stress [Pa]
    a : float
        Crack radius [m]
    E : float
        Young's modulus [Pa]
    nu : float
        Poisson's ratio [-]
    use_interpolation : bool
        If True, use accurate Bessel integral + interpolation (like Mathematica)
        If False, use fast elliptic integral approximation
    
    Returns:
    --------
    displacements : array, shape (N, 3)
        Displacements [u_x, u_y, u_z] for each node [m]
    """
    N = len(node_coords)
    displacements = np.zeros((N, 3))
    
    for i, coords in enumerate(node_coords):
        displacements[i] = sneddon_displacement_cartesian(
            sigma_app, a, E, nu, coords, use_interpolation
        )
    
    return displacements


# For testing
if __name__ == "__main__":
    # Test parameters
    sigma_app = 1.0e6  # 1 MPa
    a = 1.0e-3  # 1 mm crack radius
    E = 3.2e9  # 3.2 GPa (PMMA)
    nu = 0.35  # Poisson's ratio
    
    # Test point outside crack
    point = np.array([2.0e-3, 2.0e-3, 0.5e-3])
    disp = sneddon_displacement_cartesian(sigma_app, a, E, nu, point)
    print(f"Point: {point}")
    print(f"Displacement: {disp}")
    print(f"  u_x = {disp[0]*1e9:.3f} nm")
    print(f"  u_y = {disp[1]*1e9:.3f} nm")
    print(f"  u_z = {disp[2]*1e9:.3f} nm")
    
    # Test point on crack surface
    point_crack = np.array([0.5e-3, 0.0, 0.0])
    disp_crack = sneddon_displacement_cartesian(sigma_app, a, E, nu, point_crack)
    print(f"\nPoint on crack: {point_crack}")
    print(f"Displacement: {disp_crack}")
    print(f"  COD = {disp_crack[2]*1e6:.3f} μm")
