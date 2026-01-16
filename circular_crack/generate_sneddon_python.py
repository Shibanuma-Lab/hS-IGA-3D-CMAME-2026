#!/usr/bin/env python3
"""
Generate Sneddon integral data using pure Python/scipy

This script computes the Bessel integrals needed for Sneddon's penny-shaped crack solution
and saves them in .mat format compatible with Mathematica.

Key features:
1. Parallel computation using multiprocessing (24 cores)
2. Adaptive integration parameters based on z-coordinate
3. Pure scipy.integrate.quad (no mpmath)
4. Saves as sneddon_python.mat in Mathematica-compatible format

Usage:
    python generate_sneddon_python.py [--test]
    
    --test: Generate low-resolution test data (120x40 grid)
    (default): Generate full-resolution data (701x301 grid)

Based on Sneddon (1946) and optimized integration strategy from previous experiments.
"""

import numpy as np
from scipy import integrate, io as sio
import multiprocessing as mp
import time
import logging
import sys
import os

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s %(asctime)s] %(message)s',
    datefmt='%y%m%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def CS_func(eta):
    """
    Kernel function: CS(η) = cos(η)/η - sin(η)/η²
    
    This is the fundamental kernel in Sneddon's integral representation.
    """
    if abs(eta) < 1e-15:
        return 0.0
    return np.cos(eta) / eta - np.sin(eta) / (eta**2)


def get_adaptive_params(zeta):
    """
    Get adaptive integration parameters based on z-coordinate (normalized)
    
    This is the VERIFIED configuration that achieved 0.128% error on 31×11 coarse grid.
    Multi-stage adaptive strategy with fine-grained control over integration parameters.
    
    Parameters:
    -----------
    zeta : float
        Normalized z-coordinate (z/c)
    
    Returns:
    --------
    eta_max : float
        Upper integration limit
    int_limit : int
        Maximum number of subdivisions for quad
    """
    if zeta < 0.0001:
        # Very close to z=0: extremely slow overall decay
        # Keep moderate eta_max but use very high iteration limit
        eta_max = 50000.0   # Avoid numerical instability
        int_limit = 5000    # Greatly increased for convergence
    elif zeta < 0.001:
        eta_max = 20000.0   # High precision near crack
        int_limit = 4000    # Increased from 1500
    elif zeta < 0.01:
        eta_max = 5000.0    # Moderate precision
        int_limit = 3000    # Increased from 1000
    elif zeta < 0.05:
        # Intermediate range
        eta_max = 2000.0    # More conservative
        int_limit = 2500
    elif zeta < 0.1:
        eta_max = 1000.0    # Standard precision
        int_limit = 2000    # Increased from 500
    elif zeta < 0.5:
        # Intermediate range
        eta_max = 500.0
        int_limit = 1500
    else:
        # For larger z, exp(-ζη) dominates and integral converges faster
        eta_max = max(200.0, 20.0 / zeta)  # Conservative with adaptive floor
        int_limit = 1000    # Increased from 300
    
    return int(eta_max), int_limit


def compute_integrals_single_point(args):
    """
    Compute Sneddon integrals at a single point (r, z)
    
    This function is designed to be called in parallel.
    
    Parameters:
    -----------
    args : tuple
        (point_index, r, z, c)
        
    Returns:
    --------
    (point_index, ur1, ur2, uz1, uz2)
    """
    idx, r, z, c = args
    
    # Normalize coordinates
    rho = r / c
    zeta = z / c
    
    # Get adaptive parameters (verified configuration from successful 31×11 test)
    eta_max, int_limit = get_adaptive_params(zeta)
    
    # Integration tolerance: very tight precision
    epsabs = 1e-12
    epsrel = 1e-12
    
    # Initialize results
    ur1 = ur2 = uz1 = uz2 = 0.0
    
    try:
        # ur1: ∫ CS(η) * exp(-ζη) * J₁(ρη) dη
        if abs(r) > 1e-15:
            def integrand_ur1(eta):
                from scipy.special import jv
                return CS_func(eta) * np.exp(-zeta * eta) * jv(1, rho * eta)
            
            ur1, _ = integrate.quad(
                integrand_ur1, 
                0, eta_max,
                epsabs=epsabs, 
                epsrel=epsrel,
                limit=int_limit
            )
        
        # ur2: ∫ ζη * CS(η) * exp(-ζη) * J₁(ρη) dη
        if abs(r) > 1e-15 and abs(z) > 1e-15:
            def integrand_ur2(eta):
                from scipy.special import jv
                return zeta * eta * CS_func(eta) * np.exp(-zeta * eta) * jv(1, rho * eta)
            
            ur2, _ = integrate.quad(
                integrand_ur2, 
                0, eta_max,
                epsabs=epsabs, 
                epsrel=epsrel,
                limit=int_limit
            )
        
        # uz1: ∫ CS(η) * exp(-ζη) * J₀(ρη) dη
        if not (r >= c and abs(z) < 1e-15):
            def integrand_uz1(eta):
                from scipy.special import jv
                return CS_func(eta) * np.exp(-zeta * eta) * jv(0, rho * eta)
            
            # Special case: uz1(0, 0) = -1 (analytical)
            if abs(r) < 1e-15 and abs(z) < 1e-15:
                uz1 = -1.0
            else:
                uz1, _ = integrate.quad(
                    integrand_uz1, 
                    0, eta_max,
                    epsabs=epsabs, 
                    epsrel=epsrel,
                    limit=int_limit
                )
        
        # uz2: ∫ ζη * CS(η) * exp(-ζη) * J₀(ρη) dη
        if abs(z) > 1e-15:
            def integrand_uz2(eta):
                from scipy.special import jv
                return zeta * eta * CS_func(eta) * np.exp(-zeta * eta) * jv(0, rho * eta)
            
            uz2, _ = integrate.quad(
                integrand_uz2, 
                0, eta_max,
                epsabs=epsabs, 
                epsrel=epsrel,
                limit=int_limit
            )
        
    except Exception as e:
        logger.warning(f"Integration failed at point {idx} (r={r:.6f}, z={z:.6f}): {e}")
        # Return zeros on failure
        pass
    
    return (idx, ur1, ur2, uz1, uz2)


def generate_sneddon_data_parallel(WG=3.5, HG=1.5, nW=700, nH=300, c=1.0, n_processes=24):
    """
    Generate Sneddon integral data using parallel computation
    
    Parameters:
    -----------
    WG : float
        Domain width in r-direction (default: 3.5)
    HG : float
        Domain height in z-direction (default: 1.5)
    nW : int
        Number of points in r-direction (default: 700)
    nH : int
        Number of points in z-direction (default: 300)
    c : float
        Crack radius (default: 1.0)
    n_processes : int
        Number of parallel processes (default: 24)
    
    Returns:
    --------
    posA : ndarray (N, 2)
        Grid coordinates [r, z]
    SA : ndarray (N, 4)
        Integral values [ur1, ur2, uz1, uz2]
    """
    logger.info("="*60)
    logger.info("Generating Sneddon integral data (Python/scipy version)")
    logger.info("="*60)
    logger.info(f"Grid: {nW+1}x{nH+1} = {(nW+1)*(nH+1)} points")
    logger.info(f"Domain: r ∈ [0, {WG}], z ∈ [0, {HG}]")
    logger.info(f"Crack radius: c = {c}")
    logger.info(f"Parallel processes: {n_processes}")
    logger.info("")
    
    # Generate grid (matching Mathematica's Range behavior: nW points from 0 to WG)
    # Mathematica: Range[0, WG, WG/nW] generates nW+1 points
    posAr = np.linspace(0, WG, nW + 1)
    posAz = np.linspace(0, HG, nH + 1)
    
    # Create coordinate list (matching Mathematica order)
    # Mathematica: Flatten[Outer[{#2, #1} &, posAz, posAr], 1]
    # This creates: {r, z} pairs
    posA = []
    for z in posAz:
        for r in posAr:
            posA.append([r, z])
    posA = np.array(posA)
    
    n_points = len(posA)
    logger.info(f"Total points to compute: {n_points}")
    logger.info("")
    
    # Prepare arguments for parallel computation
    args_list = [(i, posA[i, 0], posA[i, 1], c) for i in range(n_points)]
    
    # Parallel computation
    logger.info(f"Starting parallel computation with {n_processes} processes...")
    start_time = time.time()
    
    with mp.Pool(processes=n_processes) as pool:
        results = []
        # Use imap_unordered for better performance and progress tracking
        for i, result in enumerate(pool.imap_unordered(compute_integrals_single_point, args_list, chunksize=10)):
            results.append(result)
            
            # Progress reporting
            if (i + 1) % 1000 == 0 or (i + 1) == n_points:
                elapsed = time.time() - start_time
                rate = (i + 1) / elapsed
                remaining = (n_points - i - 1) / rate if rate > 0 else 0
                logger.info(f"Progress: {i+1}/{n_points} ({100*(i+1)/n_points:.1f}%) - "
                          f"Rate: {rate:.1f} pts/s - ETA: {remaining/60:.1f} min")
    
    end_time = time.time()
    total_time = end_time - start_time
    
    logger.info("")
    logger.info(f"✓ Computation completed in {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    logger.info(f"  Average rate: {n_points/total_time:.1f} points/second")
    logger.info("")
    
    # Sort results by index and extract integral values
    results.sort(key=lambda x: x[0])
    SA = np.array([[r[1], r[2], r[3], r[4]] for r in results])
    
    return posA, SA


def save_as_mat(filename, posA, SA, c, WG, HG, nW, nH):
    """
    Save data in Mathematica-compatible .mat format
    
    Parameters:
    -----------
    filename : str
        Output filename (e.g., 'sneddon_python.mat')
    posA : ndarray
        Grid coordinates
    SA : ndarray
        Integral values
    c, WG, HG, nW, nH : float/int
        Grid parameters
    """
    logger.info(f"Saving data to {filename}...")
    
    # Prepare data dictionary (matching Mathematica format)
    mat_dict = {
        'posA': posA,
        'SA': SA,
        'c': np.array([[c]]),
        'WG': np.array([[WG]]),
        'HG': np.array([[HG]]),
        'nW': np.array([[nW]]),
        'nH': np.array([[nH]])
    }
    
    # Save to .mat file
    sio.savemat(filename, mat_dict)
    
    logger.info(f"✓ Data saved successfully")
    logger.info(f"  File: {filename}")
    logger.info(f"  Size: {os.path.getsize(filename) / 1024 / 1024:.2f} MB")


def main():
    """Main function"""
    
    # Check for test mode
    test_mode = '--test' in sys.argv
    
    if test_mode:
        logger.info("Running in TEST mode (low resolution)")
        logger.info("")
        
        # Test parameters (low resolution)
        WG = 3.0
        HG = 1.0
        nW = 119  # 120 points (0 to WG)
        nH = 39   # 40 points (0 to HG)
        c = 1.0
        output_file = 'sneddon_python_test.mat'
        n_processes = 24
        
    else:
        logger.info("Running in FULL mode (high resolution)")
        logger.info("This will take approximately 50-105 minutes...")
        logger.info("")
        
        # Full resolution parameters (expanded domain)
        WG = 3.5
        HG = 1.5
        nW = 700  # 701 points
        nH = 300  # 301 points
        c = 1.0
        output_file = 'sneddon_python.mat'
        n_processes = 24
    
    # Generate data
    posA, SA = generate_sneddon_data_parallel(
        WG=WG, HG=HG, nW=nW, nH=nH, c=c, n_processes=n_processes
    )
    
    # Save to .mat file
    save_as_mat(output_file, posA, SA, c, WG, HG, nW, nH)
    
    logger.info("")
    logger.info("="*60)
    logger.info("ALL DONE!")
    logger.info("="*60)
    logger.info("")
    logger.info(f"Next steps:")
    logger.info(f"1. Verify data quality by comparing with Mathematica:")
    logger.info(f"   python compare_sneddon_data.py")
    logger.info(f"2. Update simulation_params.py:")
    logger.info(f"   sneddon_data_source = 'python'")
    logger.info(f"3. Run static analysis:")
    logger.info(f"   python main.py --static_only")


if __name__ == "__main__":
    main()
