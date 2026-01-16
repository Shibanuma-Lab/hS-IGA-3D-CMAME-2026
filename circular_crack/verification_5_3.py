#!/usr/bin/env python3
"""
Verification script for Section 5.3: Verification of the proposed strategy 
using local mesh coarsened in the crack front direction

Based on the paper's Figure 16 and 17 methodology

This script studies the effect of d_theta (crack front direction length) on
the solution accuracy. Unlike Section 5.2 where d_theta was determined by the
isotropic condition (2*sin(d_theta/2) = hL), here we independently vary d_theta.

Paper parameters (Section 5.3):
- Local mesh dimensions: WL=0.5, aL=0.25, lL=0.25, HL=0.25
- Fixed rGL = 4
- Two mesh configurations:
  a) hL = 1/48  
  b) hL = 1/96
- d_theta values to test: 1°, 3°, 6°, 10°, 15°

Total: 5 d_theta values × 2 mesh configs = 10 cases
"""

import os
import sys
import subprocess
import shutil
import numpy as np
import json
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

class Verification53:
    def __init__(self, domain_scale=0.25, skip_existing=True):
        # Paper parameters for local mesh (in normalized units)
        self.WL = domain_scale * 2
        self.domain_scale = domain_scale
        self.aL = domain_scale  # Crack radius (0.25 in paper)
        self.lL = domain_scale  # Ligament length (0.25 in paper)
        self.HL = domain_scale  # Height (0.25 in paper)
        
        # Fixed rGL value for section 5.3
        self.rGL = 4
        
        # Fixed rBL value
        self.rBL = 0.25
        
        # Two mesh configurations
        self.hL_values = [1.0/48.0, 1.0/96.0]
        
        # d_theta values to test (degrees)
        # Order: from large to small (coarse to fine in crack front direction)
        # Larger d_theta = fewer elements in crack front = smaller DOF
        # Note: Maximum d_theta limited to 9° due to IGA global mesh constraint
        # (larger d_theta causes local mesh to span multiple IGA elements, reducing accuracy)
        # All values are divisors of 90° to ensure integer number of elements in quarter circle
        # 90°/9° = 10, 90°/6° = 15, 90°/3° = 30, 90°/2° = 45, 90°/1° = 90
        self.d_theta_values = [9.0, 6.0, 3.0, 2.0, 1.0]
        
        # Skip existing results
        self.skip_existing = skip_existing
        
        # Base directory for results
        self.base_results_dir = "results/verification_5_3"
        
        # Store all run configurations
        self.run_history = []
        
    def calculate_local_elements(self, hL):
        """
        Calculate number of local elements based on paper dimensions
        
        Paper dimensions (normalized):
        - WL = 0.5 (crack width in local mesh)
        - aL = 0.25 (crack radius in local mesh)
        - lL = 0.25 (ligament length)
        - HL = 0.25 (height)
        
        Element counts:
        - aL_elements = aL / hL
        - lL_elements = lL / hL
        - HL_elements = HL / hL
        
        Args:
            hL: Local element size
            
        Returns:
            Dictionary with element counts
        """
        aL_elements = int(round(self.aL / hL))
        lL_elements = int(round(self.lL / hL))
        HL_elements = int(round(self.HL / hL))
        
        return {
            'aL': aL_elements,
            'lL': lL_elements,
            'HL': HL_elements
        }
    
    def estimate_DOF(self, hL):
        """
        Estimate degrees of freedom for a given configuration
        
        Args:
            hL: Local element size
            
        Returns:
            Estimated DOF count
        """
        # Local mesh elements
        local_elems = self.calculate_local_elements(hL)
        aL_el = local_elems['aL']
        lL_el = local_elems['lL']
        HL_el = local_elems['HL']
        
        # Estimate local mesh DOF (approximate)
        local_DOF = aL_el * lL_el * HL_el * 8
        
        # Global mesh size
        hG = hL * self.rGL
        
        # Global domain (normalized): WidthG=2.0, HeightG=1.0
        WidthG = 2.0
        HeightG = 1.0
        
        # Global mesh control points (rough estimate)
        mu_G = 0.99 ** 0.5
        nPtsX = int(np.ceil(WidthG * mu_G / hG)) + 2
        nPtsY = nPtsX
        nPtsZ = int(np.ceil(HeightG * mu_G / hG)) + 2
        
        global_DOF = nPtsX * nPtsY * nPtsZ * 3
        
        # Total DOF
        total_DOF = local_DOF + global_DOF
        
        return total_DOF
    
    def update_const_files(self, hL, d_theta, local_elems):
        """
        Update const files with current parameters
        
        Args:
            hL: Local element size
            d_theta: Angular resolution (degrees)
            local_elems: Dictionary with element counts
        """
        import re
        
        # Update const_local_mesh.py
        const_local_mesh_path = "const/const_local_mesh.py"
        
        with open(const_local_mesh_path, 'r') as f:
            lines = f.readlines()
        
        new_lines = []
        for line in lines:
            # Update static mode values using regex
            if re.match(r'^hL_static\s*=', line):
                new_lines.append(f"hL_static = {hL}  # Normalized element size (dimensionless)\n")
            elif re.match(r'^aL_static\s*=', line):
                new_lines.append(f"aL_static = {local_elems['aL']}\n")
            elif re.match(r'^lL_static\s*=', line):
                new_lines.append(f"lL_static = {local_elems['lL']}\n")
            elif re.match(r'^HL_static\s*=', line):
                new_lines.append(f"HL_static = {local_elems['HL']}\n")
            elif re.match(r'^d_theta\s*=', line):
                new_lines.append(f"d_theta = {d_theta:.6f}  # Angular resolution [degrees]\n")
            else:
                new_lines.append(line)
        
        with open(const_local_mesh_path, 'w') as f:
            f.writelines(new_lines)
        
        # Update const_global_mesh.py
        const_global_mesh_path = "const/const_global_mesh.py"
        
        with open(const_global_mesh_path, 'r') as f:
            lines = f.readlines()
        
        new_lines = []
        for line in lines:
            if re.match(r'^rGL\s*=', line):
                new_lines.append(f"rGL = {self.rGL}  # Ratio of global to local element size (hG/hL)\n")
            elif re.match(r'^rBL\s*=', line):
                new_lines.append(f"rBL = {self.rBL}  # Ratio of background to local element size (hB/hL)\n")
            else:
                new_lines.append(line)
        
        with open(const_global_mesh_path, 'w') as f:
            f.writelines(new_lines)
    
    def create_result_folder(self, hL, d_theta):
        """
        Create result folder with naming convention: 
        verification_5_3/hL_{hL}/dtheta_{d_theta}
        
        Args:
            hL: Local element size
            d_theta: Angular resolution (degrees)
            
        Returns:
            Tuple of (path to result folder, whether it already existed)
        """
        # Create hL folder
        hL_folder = os.path.join(self.base_results_dir, f"hL_{hL:.8f}")
        os.makedirs(hL_folder, exist_ok=True)
        
        # Create specific d_theta folder
        run_folder_name = f"dtheta_{d_theta:.1f}"
        run_folder = os.path.join(hL_folder, run_folder_name)
        
        # Check if folder already exists and has results
        already_exists = self._check_result_exists(run_folder, hL)
        
        os.makedirs(run_folder, exist_ok=True)
        
        return run_folder, already_exists
    
    def _check_result_exists(self, run_folder, hL):
        """
        Check if a result folder already contains valid simulation results
        
        Args:
            run_folder: Path to the result folder
            hL: Local element size (used to determine step number)
            
        Returns:
            True if valid results exist, False otherwise
        """
        if not os.path.exists(run_folder):
            return False
        
        # Calculate expected step number
        c = 1.0
        step = int(round(c / hL))
        step_str = f"{step:05d}"
        
        # Check if result directory exists
        result_dir = os.path.join(run_folder, f"step{step_str}")
        if not os.path.exists(result_dir):
            return False
        
        # Check if essential files exist
        essential_files = [
            os.path.join(result_dir, "node.g.dat"),
            os.path.join(result_dir, "node.l.dat"),
            os.path.join(result_dir, "log", "u.g.dat"),
            os.path.join(result_dir, "log", "u_gl.l.dat")
        ]
        
        for f in essential_files:
            if not os.path.exists(f):
                return False
        
        return True
    
    def save_run_config(self, run_folder, hL, hG, d_theta, local_elems, est_DOF):
        """
        Save run configuration to JSON file
        
        Args:
            run_folder: Path to result folder
            hL: Local element size
            hG: Global element size
            d_theta: Angular resolution (degrees)
            local_elems: Dictionary with element counts
            est_DOF: Estimated degrees of freedom
        """
        config = {
            'timestamp': datetime.now().isoformat(),
            'section': '5.3',
            'description': 'Verification using local mesh coarsened in crack front direction',
            'domain_scale': self.domain_scale,
            'hL': hL,
            'hG': hG,
            'rGL': self.rGL,
            'rBL': self.rBL,
            'd_theta': d_theta,
            'local_elements': local_elems,
            'estimated_DOF': est_DOF,
            'local_dimensions': {
                'WL': self.WL,
                'aL': self.aL,
                'lL': self.lL,
                'HL': self.HL
            }
        }
        
        config_file = os.path.join(run_folder, "run_config.json")
        with open(config_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        # Also add to history
        self.run_history.append(config)
    
    def run_simulation(self, run_folder, hL):
        """
        Run the simulation
        
        Args:
            run_folder: Path to result folder
            hL: Local element size
            
        Returns:
            True if successful, False otherwise
        """
        print(f"\n  Running simulation...")
        
        # Calculate step number
        c = 1.0
        step = int(round(c / hL))
        
        # Run main.py with static_only mode
        cmd = ['python3', 'main.py', '--static_only']
        
        try:
            # Run simulation
            print(f"  Command: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            
            print(f"  Simulation completed successfully")
            
            # Copy results to result folder
            step_str = f"{step:05d}"
            inputfiles_dir = f"inputfiles/step{step_str}"
            results_dir = f"results/step{step_str}"
            
            if os.path.exists(inputfiles_dir):
                dest_inputfiles = os.path.join(run_folder, f"step{step_str}")
                if os.path.exists(dest_inputfiles):
                    shutil.rmtree(dest_inputfiles)
                shutil.copytree(inputfiles_dir, dest_inputfiles)
                print(f"  Copied inputfiles to {dest_inputfiles}")
            
            if os.path.exists(results_dir):
                dest_results = os.path.join(run_folder, f"step{step_str}")
                if not os.path.exists(dest_results):
                    os.makedirs(dest_results)
                
                # Copy results
                for item in os.listdir(results_dir):
                    src = os.path.join(results_dir, item)
                    dst = os.path.join(dest_results, item)
                    if os.path.isdir(src):
                        if os.path.exists(dst):
                            shutil.rmtree(dst)
                        shutil.copytree(src, dst)
                    else:
                        shutil.copy2(src, dst)
                print(f"  Copied results to {dest_results}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            print(f"  ERROR: Simulation failed!")
            print(f"  Return code: {e.returncode}")
            print(f"  stdout: {e.stdout}")
            print(f"  stderr: {e.stderr}")
            return False
        except Exception as e:
            print(f"  ERROR: {str(e)}")
            return False
    
    def run_single_case(self, hL, d_theta):
        """
        Run a single case with given hL and d_theta
        
        Args:
            hL: Local element size
            d_theta: Angular resolution (degrees)
            
        Returns:
            True if successful, False otherwise
        """
        # Calculate parameters
        hG = hL * self.rGL
        local_elems = self.calculate_local_elements(hL)
        est_DOF = self.estimate_DOF(hL)
        
        print(f"\n{'='*80}")
        print(f"Running case: hL={hL:.8f} (1/{1/hL:.0f}), d_theta={d_theta}°")
        print(f"{'='*80}")
        print(f"  rGL = {self.rGL}")
        print(f"  hL = {hL:.8f}")
        print(f"  hG = {hG:.8f}")
        print(f"  d_theta = {d_theta}°")
        print(f"  Local elements: aL={local_elems['aL']}, lL={local_elems['lL']}, HL={local_elems['HL']}")
        print(f"  Estimated DOF: {est_DOF:,}")
        
        # Create result folder
        run_folder, already_exists = self.create_result_folder(hL, d_theta)
        print(f"  Result folder: {run_folder}")
        
        if already_exists and self.skip_existing:
            print(f"  SKIPPED: Results already exist")
            # Still save config for record keeping
            self.save_run_config(run_folder, hL, hG, d_theta, local_elems, est_DOF)
            return True
        
        # Update const files
        print(f"  Updating const files...")
        self.update_const_files(hL, d_theta, local_elems)
        
        # Save configuration
        self.save_run_config(run_folder, hL, hG, d_theta, local_elems, est_DOF)
        
        # Run simulation
        success = self.run_simulation(run_folder, hL)
        
        if success:
            print(f"\n  ✓ Case completed successfully")
        else:
            print(f"\n  ✗ Case failed")
        
        return success
    
    def run_all_cases(self):
        """
        Run all cases for section 5.3
        
        Total cases: 5 d_theta values × 2 hL values = 10 cases
        """
        print(f"\n{'='*80}")
        print(f"Section 5.3 Verification")
        print(f"Study of d_theta effect on solution accuracy")
        print(f"{'='*80}")
        print(f"\nConfiguration:")
        print(f"  rGL = {self.rGL}")
        print(f"  rBL = {self.rBL}")
        print(f"  hL values: {self.hL_values}")
        print(f"  d_theta values: {self.d_theta_values}°")
        print(f"  Total cases: {len(self.hL_values) * len(self.d_theta_values)}")
        print(f"  Skip existing: {self.skip_existing}")
        print(f"  Results directory: {self.base_results_dir}")
        
        # Create base results directory
        os.makedirs(self.base_results_dir, exist_ok=True)
        
        # Track statistics
        total_cases = len(self.hL_values) * len(self.d_theta_values)
        completed = 0
        failed = 0
        skipped = 0
        
        # Run all cases
        case_num = 0
        for hL in self.hL_values:
            for d_theta in self.d_theta_values:
                case_num += 1
                print(f"\n{'='*80}")
                print(f"Case {case_num}/{total_cases}")
                print(f"{'='*80}")
                
                # Check if already exists before running
                run_folder, already_exists = self.create_result_folder(hL, d_theta)
                
                if already_exists and self.skip_existing:
                    skipped += 1
                    print(f"SKIPPED: hL={hL:.8f}, d_theta={d_theta}° (results exist)")
                    continue
                
                success = self.run_single_case(hL, d_theta)
                
                if success:
                    completed += 1
                else:
                    failed += 1
        
        # Print summary
        print(f"\n{'='*80}")
        print(f"Section 5.3 Verification Summary")
        print(f"{'='*80}")
        print(f"Total cases: {total_cases}")
        print(f"Completed: {completed}")
        print(f"Skipped: {skipped}")
        print(f"Failed: {failed}")
        
        # Save run history
        history_file = os.path.join(self.base_results_dir, "run_history.json")
        with open(history_file, 'w') as f:
            json.dump({
                'timestamp': datetime.now().isoformat(),
                'section': '5.3',
                'total_cases': total_cases,
                'completed': completed,
                'skipped': skipped,
                'failed': failed,
                'runs': self.run_history
            }, f, indent=2)
        print(f"\nRun history saved to: {history_file}")


def main():
    """Main entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Run Section 5.3 verification cases"
    )
    parser.add_argument(
        '--domain-scale',
        type=float,
        default=0.25,
        help='Domain scale factor (default: 0.25)'
    )
    parser.add_argument(
        '--no-skip',
        action='store_true',
        help='Re-run cases even if results exist'
    )
    parser.add_argument(
        '--hL',
        type=str,
        choices=['48', '96', 'all'],
        default='all',
        help='Which hL config to run: 48 (1/48), 96 (1/96), or all (default: all)'
    )
    parser.add_argument(
        '--dtheta',
        type=float,
        help='Run only specific d_theta value (degrees)'
    )
    
    args = parser.parse_args()
    
    # Create verification object
    verif = Verification53(
        domain_scale=args.domain_scale,
        skip_existing=not args.no_skip
    )
    
    # Filter hL values if specified
    if args.hL != 'all':
        if args.hL == '48':
            verif.hL_values = [1.0/48.0]
        elif args.hL == '96':
            verif.hL_values = [1.0/96.0]
    
    # Filter d_theta values if specified
    if args.dtheta is not None:
        if args.dtheta in verif.d_theta_values:
            verif.d_theta_values = [args.dtheta]
        else:
            print(f"ERROR: d_theta={args.dtheta} not in valid values: {verif.d_theta_values}")
            return 1
    
    # Run all cases
    verif.run_all_cases()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
