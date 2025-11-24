#!/usr/bin/env python3
"""
Quick test of Sneddon solution implementation
Compare elliptic vs Bessel integral methods
"""

import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(__file__))

print("=" * 60)
print("Testing Sneddon Solution Implementation")
print("=" * 60)

# First, generate a small test interpolation dataset
print("\n[1/3] Generating test interpolation data (this may take a few minutes)...")
print("      Grid: 60x20 (very coarse for quick testing)")

from sneddon_precompute import save_interpolation_data

try:
    save_interpolation_data(
        filename='sneddon_interpolation_test.npz',
        WG=3.0, HG=1.0, 
        nW=60, nH=20,  # Very coarse grid for quick test
        c=1.0
    )
    print("✓ Test interpolation data generated")
except Exception as e:
    print(f"✗ Failed to generate interpolation data: {e}")
    sys.exit(1)

# Now test the solution
print("\n[2/3] Testing solution accuracy...")

from sneddon_solution import sneddon_displacement_cartesian

# Test parameters (matching your setup)
sigma_app = 1.0e6  # 1 MPa
a = 1.0  # 1 m (normalized)
E = 3.2e9  # 3.2 GPa
nu = 0.35

# Test points
test_points = [
    ([0.5, 0.0, 0.0], "On crack surface (r=0.5a)"),
    ([0.0, 0.0, 0.5], "On z-axis above crack"),
    ([1.5, 0.0, 0.0], "Outside crack (r=1.5a, z=0)"),
    ([1.0, 1.0, 0.5], "General point (x=y=a, z=0.5a)"),
]

print("\nComparing methods:")
print("  Method A: Elliptic integral (fast, approximate)")
print("  Method B: Bessel integral + interpolation (accurate, like Mathematica)")
print()

for point, description in test_points:
    print(f"\n{description}")
    print(f"  Point: x={point[0]:.2f}, y={point[1]:.2f}, z={point[2]:.2f}")
    
    # Method A: Elliptic integral
    try:
        disp_elliptic = sneddon_displacement_cartesian(
            sigma_app, a, E, nu, point, use_interpolation=False
        )
        print(f"  Method A: u_x={disp_elliptic[0]*1e6:.4f} μm, "
              f"u_y={disp_elliptic[1]*1e6:.4f} μm, "
              f"u_z={disp_elliptic[2]*1e6:.4f} μm")
    except Exception as e:
        print(f"  Method A: FAILED - {e}")
        disp_elliptic = None
    
    # Method B: Bessel integral + interpolation
    try:
        disp_bessel = sneddon_displacement_cartesian(
            sigma_app, a, E, nu, point, use_interpolation=True
        )
        print(f"  Method B: u_x={disp_bessel[0]*1e6:.4f} μm, "
              f"u_y={disp_bessel[1]*1e6:.4f} μm, "
              f"u_z={disp_bessel[2]*1e6:.4f} μm")
    except Exception as e:
        print(f"  Method B: FAILED - {e}")
        disp_bessel = None
    
    # Compare
    if disp_elliptic is not None and disp_bessel is not None:
        diff = np.abs(disp_bessel - disp_elliptic)
        rel_error = diff / (np.abs(disp_bessel) + 1e-15) * 100
        print(f"  Difference: Δu_x={diff[0]*1e6:.4f} μm ({rel_error[0]:.1f}%), "
              f"Δu_y={diff[1]*1e6:.4f} μm ({rel_error[1]:.1f}%), "
              f"Δu_z={diff[2]*1e6:.4f} μm ({rel_error[2]:.1f}%)")

# Test COD (Crack Opening Displacement) - analytical formula
print("\n[3/3] Checking COD (Crack Opening Displacement)...")
print("      Analytical formula: COD(r) = (2σ(1-ν²)/E)√(a²-r²)")

r_test = 0.0  # Center of crack
point_crack = [r_test, 0.0, 0.0]

disp_numerical = sneddon_displacement_cartesian(
    sigma_app, a, E, nu, point_crack, use_interpolation=True
)

# Analytical COD at crack center
COD_analytical = (2 * sigma_app * (1 - nu**2) / E) * np.sqrt(a**2 - r_test**2)

print(f"\n  Point: r={r_test:.2f}, z=0 (crack center)")
print(f"  Numerical COD:  {disp_numerical[2]*1e6:.4f} μm")
print(f"  Analytical COD: {COD_analytical*1e6:.4f} μm")
print(f"  Error: {abs(disp_numerical[2] - COD_analytical)*1e6:.4f} μm "
      f"({abs(disp_numerical[2] - COD_analytical)/COD_analytical*100:.2f}%)")

print("\n" + "=" * 60)
print("Test Complete!")
print("=" * 60)
print("\nNext steps:")
print("  1. For production use, generate full resolution data:")
print("     ./generate_sneddon_data.sh")
print("     (Choose option 2 for full 1200x400 grid)")
print()
print("  2. Run static analysis:")
print("     python3 main.py --static_only")
print()
