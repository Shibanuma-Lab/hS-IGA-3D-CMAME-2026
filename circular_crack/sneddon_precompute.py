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
    
    # Create interpolator
    interpolator = RegularGridInterpolator(
        (r_grid, z_grid), 
        data,
        method='linear',  # Linear interpolation (InterpolationOrder->1 in Mathematica)
        bounds_error=False,
        fill_value=None  # Extrapolate if needed
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
    
    Returns:
    --------
    interpolator : RegularGridInterpolator
    c : float (crack radius used for normalization)
    """
    # Load data
    data_file = np.load(filename)
    r_grid = data_file['r_grid']
    z_grid = data_file['z_grid']
    data = data_file['data']
    c = float(data_file['c'])
    
    # Create interpolator
    interpolator = RegularGridInterpolator(
        (r_grid, z_grid),
        data,
        method='linear',
        bounds_error=False,
        fill_value=None
    )
    
    logger.info(f"Loaded interpolation data from {filename}")
    logger.info(f"  Grid: {len(r_grid)}x{len(z_grid)}, c={c}")
    
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
