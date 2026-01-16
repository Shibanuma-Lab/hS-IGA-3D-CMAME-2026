"""
Precompute Sneddon solution using Bessel integral formulation
Based on Mathematica SneddonApp.mx implementation

This script generates interpolation data for fast Sneddon solution evaluation.
Run once to generate 'sneddon_interpolation.npz', then use in sneddon_solution.py
"""

import numpy as np
from scipy import integrate, special
from scipy.interpolate import RegularGridInterpolator
import pickle
import os
from utils.logger import logger


class BoundaryCheckInterpolator:
    """
    Wrapper for RegularGridInterpolator that warns when points are outside domain
    """
    def __init__(self, interpolator, r_bounds, z_bounds, source='MAT file'):
        self.interpolator = interpolator
        self.r_min, self.r_max = r_bounds
        self.z_min, self.z_max = z_bounds
        self.source = source
        self._warning_issued = False
        
    def __call__(self, points):
        """Call interpolator with boundary checking"""
        # Convert to array if needed
        points = np.asarray(points)
        
        # Check if point is outside domain
        if points.ndim == 1:
            r, z = points
            if (r < self.r_min or r > self.r_max or z < self.z_min or z > self.z_max):
                if not self._warning_issued:
                    logger.warning(f"⚠ INTERPOLATION WARNING: Point outside {self.source} domain!")
                    logger.warning(f"  Domain: r ∈ [{self.r_min}, {self.r_max}], z ∈ [{self.z_min}, {self.z_max}]")
                    logger.warning(f"  Point: r={r:.4f}, z={z:.4f}")
                    logger.warning(f"  Using extrapolation (reduced accuracy)")
                    logger.warning(f"  Consider regenerating MAT file with larger domain")
                    self._warning_issued = True
        else:
            # Check multiple points
            rs = points[:, 0]
            zs = points[:, 1]
            out_of_bounds = ((rs < self.r_min) | (rs > self.r_max) | 
                           (zs < self.z_min) | (zs > self.z_max))
            if np.any(out_of_bounds) and not self._warning_issued:
                n_out = np.sum(out_of_bounds)
                r_max_actual = np.max(rs)
                z_max_actual = np.max(zs)
                logger.warning(f"⚠ INTERPOLATION WARNING: {n_out} points outside {self.source} domain!")
                logger.warning(f"  Domain: r ∈ [{self.r_min}, {self.r_max}], z ∈ [{self.z_min}, {self.z_max}]")
                logger.warning(f"  Actual range: r ∈ [0, {r_max_actual:.4f}], z ∈ [0, {z_max_actual:.4f}]")
                logger.warning(f"  Using extrapolation for out-of-bounds points (reduced accuracy)")
                logger.warning(f"  Consider regenerating MAT file with larger domain")
                self._warning_issued = True
                
        # Call the actual interpolator
        return self.interpolator(points)


def CS(eta):
    """
    Helper function: CS(η) = Cos(η)/η - Sin(η)/η²
    """
    if abs(eta) < 1e-10:
        # Taylor expansion near η=0: CS(η) ≈ 0
        return 0.0
    return np.cos(eta) / eta - np.sin(eta) / (eta**2)


def sneddon_integrals(r, z, c=1.0):
    """
    Compute Sneddon integrals using Bessel functions
    
    Based on Mathematica Sneddon0[{r, z}]:
    - ur1 = ∫ CS(η) E_ζ(η) J₁(ρη) dη
    - ur2 = ∫ ζη CS(η) E_ζ(η) J₁(ρη) dη
    - uz1 = ∫ CS(η) E_ζ(η) J₀(ρη) dη
    - uz2 = ∫ ζη CS(η) E_ζ(η) J₀(ρη) dη
    
    Where:
    - ρ = r/c (normalized radius)
    - ζ = z/c (normalized height)
    - E_ζ(η) = exp(-ζη)
    - J₀, J₁ are Bessel functions of first kind
    
    Parameters:
    -----------
    r : float
        Radial coordinate [m]
    z : float
        Axial coordinate [m]
    c : float
        Crack radius [m] (for normalization)
    
    Returns:
    --------
    [ur1, ur2, uz1, uz2] : list of floats
        Four integral values (dimensionless)
    """
    rho = r / c
    zeta = z / c
    
    # Define integrands
    def integrand_ur1(eta):
        if abs(eta) < 1e-10:
            return 0.0
        return CS(eta) * np.exp(-zeta * eta) * special.jv(1, rho * eta)
    
    def integrand_ur2(eta):
        if abs(eta) < 1e-10:
            return 0.0
        return zeta * eta * CS(eta) * np.exp(-zeta * eta) * special.jv(1, rho * eta)
    
    def integrand_uz1(eta):
        if abs(eta) < 1e-10:
            return 0.0
        return CS(eta) * np.exp(-zeta * eta) * special.jv(0, rho * eta)
    
    def integrand_uz2(eta):
        if abs(eta) < 1e-10:
            return 0.0
        return zeta * eta * CS(eta) * np.exp(-zeta * eta) * special.jv(0, rho * eta)
    
    # Integration limits
    # Mathematica uses {η, 0, Infinity}, we truncate at large value
    eta_max = max(50.0, 10.0 / max(zeta, 0.01))  # Adaptive upper limit
    
    # Compute integrals
    # ur1
    if abs(r) < 1e-10:
        ur1 = 0.0
    else:
        ur1, _ = integrate.quad(integrand_ur1, 0, eta_max, 
                                limit=100, epsabs=1e-6, epsrel=1e-6)
    
    # ur2
    if abs(r) < 1e-10:
        ur2 = 0.0
    else:
        ur2, _ = integrate.quad(integrand_ur2, 0, eta_max,
                                limit=100, epsabs=1e-6, epsrel=1e-6)
    
    # uz1
    if r >= c and abs(z) < 1e-10:
        uz1 = 0.0
    else:
        uz1, _ = integrate.quad(integrand_uz1, 0, eta_max,
                                limit=100, epsabs=1e-6, epsrel=1e-6)
    
    # uz2
    if r >= c and abs(z) < 1e-10:
        uz2 = 0.0
    else:
        uz2, _ = integrate.quad(integrand_uz2, 0, eta_max,
                                limit=100, epsabs=1e-6, epsrel=1e-6)
    
    return [ur1, ur2, uz1, uz2]


def precompute_sneddon_data(WG=3.0, HG=1.0, nW=600, nH=200, c=1.0):
    """
    Precompute Sneddon integrals on a regular grid
    
    Parameters:
    -----------
    WG : float
        Domain width in r-direction [m]
    HG : float
        Domain height in z-direction [m]
    nW : int
        Number of grid points in r-direction (default: 600, matching Mathematica)
    nH : int
        Number of grid points in z-direction (default: 200, matching Mathematica)
    c : float
        Crack radius [m]
    
    Returns:
    --------
    interpolator : RegularGridInterpolator
        Interpolation function for [ur1, ur2, uz1, uz2]
    """
    logger.info(f"Precomputing Sneddon integrals on {nW}x{nH} grid...")
    logger.info(f"  Domain: r ∈ [0, {WG}], z ∈ [0, {HG}]")
    logger.info(f"  Crack radius c = {c}")
    
    # Create grid
    r_grid = np.linspace(0, WG, nW)
    z_grid = np.linspace(0, HG, nH)
    
    # Storage for integral values
    data = np.zeros((nW, nH, 4))  # 4 integrals per point
    
    # Compute integrals at each grid point
    total_points = nW * nH
    computed = 0
    
    for i, r in enumerate(r_grid):
        for j, z in enumerate(z_grid):
            data[i, j, :] = sneddon_integrals(r, z, c)
            
            computed += 1
            if computed % 10000 == 0:
                progress = 100.0 * computed / total_points
                logger.info(f"  Progress: {progress:.1f}% ({computed}/{total_points})")
    
    logger.info("Precomputation complete!")
    
    # Create base interpolator
    base_interpolator = RegularGridInterpolator(
        (r_grid, z_grid), 
        data,
        method='linear',  # Linear interpolation (InterpolationOrder->1 in Mathematica)
        bounds_error=False,
        fill_value=None  # Extrapolate if needed
    )
    
    # Wrap with boundary checking
    interpolator = BoundaryCheckInterpolator(
        base_interpolator,
        r_bounds=(0, WG),
        z_bounds=(0, HG),
        source='precomputed data'
    )
    
    return interpolator, r_grid, z_grid, data


def save_interpolation_data(filename='sneddon_interpolation.npz', 
                            WG=3.0, HG=1.0, nW=600, nH=200, c=1.0):
    """
    Precompute and save Sneddon interpolation data
    
    Parameters:
    -----------
    filename : str
        Output filename (.npz format)
    WG, HG, nW, nH, c : see precompute_sneddon_data()
    """
    interpolator, r_grid, z_grid, data = precompute_sneddon_data(WG, HG, nW, nH, c)
    
    # Save to file
    np.savez_compressed(
        filename,
        r_grid=r_grid,
        z_grid=z_grid,
        data=data,
        c=c,
        WG=WG,
        HG=HG
    )
    
    logger.info(f"Interpolation data saved to {filename}")
    file_size_mb = os.path.getsize(filename) / 1024 / 1024
    logger.info(f"  File size: {file_size_mb:.2f} MB")
    
    return interpolator


def load_interpolation_data(filename='sneddon_interpolation.npz'):
    """
    Load precomputed Sneddon interpolation data
    
    Can load from:
    - .npz format (Python generated)
    - .mat format (Mathematica exported)
    
    Returns:
    --------
    interpolator : RegularGridInterpolator
    c : float (crack radius used for normalization)
    """
    if filename.endswith('.mat'):
        # Load from Mathematica .mat file
        return load_from_mathematica_mat(filename)
    else:
        # Load from Python .npz file
        data_file = np.load(filename)
        r_grid = data_file['r_grid']
        z_grid = data_file['z_grid']
        data = data_file['data']
        c = float(data_file['c'])
        WG = float(data_file['WG'])
        HG = float(data_file['HG'])
        
        # Create base interpolator
        base_interpolator = RegularGridInterpolator(
            (r_grid, z_grid),
            data,
            method='linear',
            bounds_error=False,
            fill_value=None
        )
        
        # Wrap with boundary checking
        interpolator = BoundaryCheckInterpolator(
            base_interpolator,
            r_bounds=(0, WG),
            z_bounds=(0, HG),
            source=f'NPZ file ({filename})'
        )
        
        logger.info(f"Loaded interpolation data from {filename}")
        logger.info(f"  Grid: {len(r_grid)}x{len(z_grid)}, c={c}")
        logger.info(f"  Domain: r ∈ [0, {WG}], z ∈ [0, {HG}]")
        
        return interpolator, c


def load_from_mathematica_mat(filename='sneddon_SA.mat'):
    """
    Load Sneddon interpolation data from .mat file
    (Compatible with both Mathematica and Python-generated .mat files)
    
    The .mat file should contain:
    - posA: Nx2 array of (z, r) coordinates
    - SA: Nx4 array of [ur1, ur2, uz1, uz2] values
    - c: crack radius
    - WG, HG: domain size
    - nW, nH: grid dimensions
    
    Returns:
    --------
    interpolator : RegularGridInterpolator
    c : float (crack radius used for normalization)
    """
    try:
        from scipy.io import loadmat
    except ImportError:
        logger.error("scipy.io.loadmat not available. Install scipy to load .mat files")
        raise
    
    logger.info(f"Loading Sneddon interpolation data from {filename}")
    
    # Load .mat file
    mat_data = loadmat(filename)
    
    # Extract data (Mathematica uses 1-based indexing, MATLAB format stores in specific way)
    # NOTE: posA from Mathematica Outer[{#2, #1}&, posAz, posAr] is actually [r, z]!
    posA = mat_data['posA']  # Nx2: [r, z] for each point (NOT [z, r]!)
    SA = mat_data['SA']      # Nx4: [ur1, ur2, uz1, uz2] for each point
    c = float(mat_data['c'][0, 0])
    WG = float(mat_data['WG'][0, 0])
    HG = float(mat_data['HG'][0, 0])
    nW = int(mat_data['nW'][0, 0])
    nH = int(mat_data['nH'][0, 0])
    
    actual_points = len(posA)
    expected_points = nW * nH
    
    logger.info(f"  Declared grid: {nW}x{nH} = {expected_points}, c={c}")
    logger.info(f"  Actual points: {actual_points}")
    logger.info(f"  Domain: r ∈ [0, {WG}], z ∈ [0, {HG}]")
    
    # Auto-detect actual grid dimensions (Mathematica Range issue: generates nW+1 points)
    if actual_points != expected_points:
        # Try to infer actual dimensions from data
        # posA is [r, z], so column 0 is r, column 1 is z
        unique_r = np.unique(posA[:, 0])
        unique_z = np.unique(posA[:, 1])
        nW_actual = len(unique_r)
        nH_actual = len(unique_z)
        
        logger.warning(f"  Grid mismatch detected!")
        logger.warning(f"  Declared: {nW}x{nH} = {expected_points}")
        logger.warning(f"  Actual: {nW_actual}x{nH_actual} = {nW_actual*nH_actual}")
        
        if nW_actual * nH_actual == actual_points:
            logger.info(f"  ✓ Auto-corrected to {nW_actual}x{nH_actual}")
            nW = nW_actual
            nH = nH_actual
        else:
            raise ValueError(
                f"Cannot determine grid dimensions!\n"
                f"Declared: {nW}x{nH} = {expected_points}\n"
                f"Actual points: {actual_points}\n"
                f"Detected unique: {nW_actual}x{nH_actual} = {nW_actual*nH_actual}"
            )
    
    # Reshape data to regular grid
    # posA from Mathematica is flattened: Flatten[Outer[{#2, #1}&, posAz, posAr], 1]
    # Outer[{#2, #1}&, posAz, posAr] means: {r from posAr, z from posAz}
    # Flatten with level 1 means: for each z in posAz, for each r in posAr
    # So the order is: (r0,z0), (r1,z0), ..., (rN,z0), (r0,z1), (r1,z1), ...
    r_grid = np.linspace(0, WG, nW)
    z_grid = np.linspace(0, HG, nH)
    
    # Reshape SA: (nH*nW, 4) -> (nH, nW, 4) -> (nW, nH, 4)
    # Data is ordered as: for each z, all r values
    # So reshape to (nH, nW, 4) first, then transpose to (nW, nH, 4)
    data = SA.reshape(nH, nW, 4).transpose(1, 0, 2)
    
    logger.info(f"  Reshaped data to ({nW}, {nH}, 4)")
    
    # Create base interpolator
    base_interpolator = RegularGridInterpolator(
        (r_grid, z_grid),
        data,
        method='linear',
        bounds_error=False,
        fill_value=None
    )
    
    # Detect data source from filename for accurate reporting
    source = 'Python' if 'python' in filename.lower() else 'Mathematica'
    
    # Wrap with boundary checking
    interpolator = BoundaryCheckInterpolator(
        base_interpolator,
        r_bounds=(0, WG),
        z_bounds=(0, HG),
        source=f'{source} MAT file (domain: r≤{WG:.1f}, z≤{HG:.1f})'
    )
    
    logger.info(f"✓ {source} interpolation data loaded successfully from .mat file")
    
    return interpolator, c


if __name__ == "__main__":
    import sys
    
    # Check if running with reduced resolution (for testing)
    if len(sys.argv) > 1 and sys.argv[1] == '--test':
        logger.info("Running in TEST mode (low resolution)")
        save_interpolation_data(
            filename='sneddon_interpolation_test.npz',
            WG=3.0, HG=1.0, nW=120, nH=40, c=1.0
        )
    else:
        logger.info("Running in FULL mode (high resolution)")
        logger.info("This will take about 30-45 minutes...")
        logger.info("For quick test, run with --test flag")
        
        save_interpolation_data(
            filename='sneddon_interpolation.npz',
            WG=3.0, HG=1.0, nW=600, nH=200, c=1.0
        )
    
    logger.info("Done!")
