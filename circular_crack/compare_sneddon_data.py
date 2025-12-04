#!/usr/bin/env python3
"""
Compare Python-generated Sneddon data with Mathematica reference

This script loads both sneddon_python.mat and sneddon_SA.mat,
compares the integral values, and reports statistics.

Usage:
    python compare_sneddon_data.py
"""

import numpy as np
from scipy.io import loadmat
import sys
import os

def load_and_compare():
    """Load both data files and compare"""
    
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Load Mathematica data
    ma_file = os.path.join(current_dir, 'sneddon_SA.mat')
    if not os.path.exists(ma_file):
        print(f"ERROR: Mathematica data not found: {ma_file}")
        return False
    
    print("Loading Mathematica data...")
    ma_data = loadmat(ma_file)
    ma_posA = ma_data['posA']
    ma_SA = ma_data['SA']
    ma_c = float(ma_data['c'][0, 0])
    
    print(f"  Mathematica: {len(ma_posA)} points, c={ma_c}")
    
    # Load Python data
    py_file = os.path.join(current_dir, 'sneddon_python.mat')
    if not os.path.exists(py_file):
        print(f"ERROR: Python data not found: {py_file}")
        print(f"Please run: python generate_sneddon_python.py")
        return False
    
    print("Loading Python data...")
    py_data = loadmat(py_file)
    py_posA = py_data['posA']
    py_SA = py_data['SA']
    py_c = float(py_data['c'][0, 0])
    
    print(f"  Python:      {len(py_posA)} points, c={py_c}")
    print()
    
    # Check dimensions match
    if len(ma_posA) != len(py_posA):
        print(f"ERROR: Different number of points!")
        print(f"  Mathematica: {len(ma_posA)}")
        print(f"  Python:      {len(py_posA)}")
        return False
    
    # Check coordinates match
    coord_diff = np.max(np.abs(ma_posA - py_posA))
    if coord_diff > 1e-10:
        print(f"WARNING: Coordinate mismatch detected!")
        print(f"  Max difference: {coord_diff}")
    else:
        print("✓ Coordinates match perfectly")
    
    print()
    print("="*70)
    print("COMPARING INTEGRAL VALUES")
    print("="*70)
    print()
    
    # Compare each integral
    integral_names = ['ur1', 'ur2', 'uz1', 'uz2']
    
    for i, name in enumerate(integral_names):
        ma_vals = ma_SA[:, i]
        py_vals = py_SA[:, i]
        
        # Calculate differences
        abs_diff = np.abs(ma_vals - py_vals)
        rel_diff = np.abs(abs_diff / (np.abs(ma_vals) + 1e-15))
        
        # Statistics
        max_abs_diff = np.max(abs_diff)
        max_rel_diff = np.max(rel_diff)
        mean_abs_diff = np.mean(abs_diff)
        mean_rel_diff = np.mean(rel_diff)
        
        # Find worst point
        worst_idx = np.argmax(abs_diff)
        worst_r = ma_posA[worst_idx, 0]
        worst_z = ma_posA[worst_idx, 1]
        worst_ma = ma_vals[worst_idx]
        worst_py = py_vals[worst_idx]
        
        print(f"{name}:")
        print(f"  Max absolute diff: {max_abs_diff:.6e}")
        print(f"  Max relative diff: {max_rel_diff*100:.4f}%")
        print(f"  Mean absolute diff: {mean_abs_diff:.6e}")
        print(f"  Mean relative diff: {mean_rel_diff*100:.4f}%")
        print(f"  Worst point: (r={worst_r:.4f}, z={worst_z:.4f})")
        print(f"    Mathematica: {worst_ma:.10e}")
        print(f"    Python:      {worst_py:.10e}")
        print()
    
    # Overall assessment
    print("="*70)
    print("OVERALL ASSESSMENT")
    print("="*70)
    print()
    
    all_abs_diff = np.abs(ma_SA - py_SA)
    all_rel_diff = np.abs(all_abs_diff / (np.abs(ma_SA) + 1e-15))
    
    overall_max_abs = np.max(all_abs_diff)
    overall_max_rel = np.max(all_rel_diff)
    
    print(f"Maximum absolute difference: {overall_max_abs:.6e}")
    print(f"Maximum relative difference: {overall_max_rel*100:.4f}%")
    print()
    
    # Quality check
    if overall_max_rel < 0.01:  # < 1%
        print("✓ EXCELLENT: Differences < 1%")
        quality = "EXCELLENT"
    elif overall_max_rel < 0.05:  # < 5%
        print("✓ GOOD: Differences < 5%")
        quality = "GOOD"
    elif overall_max_rel < 0.10:  # < 10%
        print("⚠ ACCEPTABLE: Differences < 10%")
        quality = "ACCEPTABLE"
    else:
        print("✗ POOR: Differences > 10%")
        quality = "POOR"
    
    print()
    
    # Recommendations
    print("="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    print()
    
    if quality in ["EXCELLENT", "GOOD"]:
        print("✓ Python-generated data is suitable for production use")
        print()
        print("To use Python data:")
        print("  1. Edit const/simulation_params.py")
        print("  2. Set: sneddon_data_source = 'python'")
        print("  3. Run: python main.py --static_only")
    elif quality == "ACCEPTABLE":
        print("⚠ Python data has acceptable accuracy but consider:")
        print("  - Increasing integration limits (eta_max)")
        print("  - Using tighter tolerances (epsabs, epsrel)")
        print("  - Checking points with large z-coordinates")
    else:
        print("✗ Python data needs improvement:")
        print("  - Review integration parameters")
        print("  - Check for numerical instabilities")
        print("  - Consider using Mathematica data for now")
    
    print()
    
    return quality in ["EXCELLENT", "GOOD", "ACCEPTABLE"]


if __name__ == "__main__":
    success = load_and_compare()
    sys.exit(0 if success else 1)
