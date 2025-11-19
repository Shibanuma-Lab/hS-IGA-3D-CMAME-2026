#!/usr/bin/env python3
"""
Quick Start Script for Circular Crack Simulation

This script provides an easy way to run simulations with different configurations.
"""

import os
import sys
import argparse

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from const import simulation_params as sp
from utils.logger import logger


def setup_directories():
    """Create necessary directories"""
    dirs = ['inputfiles', 'logs', sp.OUTPUT_FOLDER]
    for d in dirs:
        os.makedirs(d, exist_ok=True)
        logger.info(f"Created directory: {d}")


def run_test():
    """Run test without solver"""
    logger.info("Running test mode (mesh generation only)")
    from test_mesh_generation import main as test_main
    return test_main()


def run_single_step(step):
    """Run a single step"""
    from main import makemeshs
    from linux_command import run
    
    logger.info(f"Running single step: {step}")
    REstart = 0 if step == sp.step_start else 1
    l, g = makemeshs(step, REstart)
    success = run(step)
    
    return success


def run_full_simulation():
    """Run full simulation from step_start to step_end"""
    from main import main
    logger.info(f"Running full simulation: steps {sp.step_start} to {sp.step_end}")
    main()


def display_config():
    """Display current configuration"""
    from const import const_global_mesh as cgm, const_local_mesh as clm, material_property as mp
    
    print("\n" + "="*70)
    print("CURRENT CONFIGURATION")
    print("="*70)
    
    print("\n[Simulation Parameters]")
    print(f"  Crack radius:        {sp.c*1e3:.3f} mm")
    print(f"  Thickness:           {sp.thi:.3f} m")
    print(f"  Crack velocity:      {sp.V:.1f} m/s")
    print(f"  Step range:          {sp.step_start} to {sp.step_end}")
    print(f"  Total steps:         {sp.stepall}")
    
    print("\n[Global Mesh]")
    print(f"  Refinement levels:   {cgm.m0}")
    print(f"  Elements in X:       {cgm.nx1}")
    print(f"  Elements in Y:       {cgm.ny1}")
    print(f"  Global/local ratio:  {cgm.rGL}")
    
    print("\n[Local Mesh]")
    print(f"  Element size:        {clm.hL*1e6:.2f} μm")
    print(f"  Crack elements:      {clm.aL}")
    print(f"  Ligament elements:   {clm.lL}")
    print(f"  Thickness elements:  {clm.HL}")
    
    print("\n[Material Properties]")
    print(f"  Young's modulus:     {mp.EE/1e9:.2f} GPa")
    print(f"  Poisson's ratio:     {mp.Nu}")
    print(f"  Density:             {mp.Rho:.1f} kg/m³")
    print(f"  Applied stress:      {mp.SigmaInfinity/1e6:.1f} MPa")
    
    print("\n[Time Integration (HHT)]")
    print(f"  Alpha:               {sp.Alpha}")
    print(f"  Beta:                {sp.Beta}")
    print(f"  Gamma:               {sp.Gamma}")
    
    print("\n[Solver Settings]")
    print(f"  OpenMP threads:      {sp.OPENMP}")
    print(f"  Integration points:  {sp.ngp}")
    print(f"  Refinement level:    {sp.nrefLlist}")
    
    print("\n" + "="*70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Circular Crack Simulation - Quick Start',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --test              Run test mode (no solver)
  %(prog)s --config            Display configuration
  %(prog)s --step 0            Run only step 0
  %(prog)s --run               Run full simulation
  %(prog)s --setup             Setup directories only
        """
    )
    
    parser.add_argument('--test', action='store_true',
                        help='Run test mode (mesh generation only)')
    parser.add_argument('--config', action='store_true',
                        help='Display current configuration')
    parser.add_argument('--step', type=int, metavar='N',
                        help='Run single step N')
    parser.add_argument('--run', action='store_true',
                        help='Run full simulation')
    parser.add_argument('--setup', action='store_true',
                        help='Setup directories only')
    
    args = parser.parse_args()
    
    # If no arguments, show help
    if len(sys.argv) == 1:
        parser.print_help()
        print("\n")
        display_config()
        return 0
    
    # Setup directories
    if args.setup or args.run or args.step is not None:
        setup_directories()
    
    # Execute requested action
    if args.config:
        display_config()
    
    if args.test:
        return run_test()
    
    if args.step is not None:
        success = run_single_step(args.step)
        return 0 if success else 1
    
    if args.run:
        run_full_simulation()
        return 0
    
    if args.setup:
        logger.info("Setup complete")
        return 0
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
