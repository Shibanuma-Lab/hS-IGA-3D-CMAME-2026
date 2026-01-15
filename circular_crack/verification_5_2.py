#!/usr/bin/env python3
"""
Verification script for Section 5.2: Verification using isotropic local mesh
Based on the paper's Figure 11 methodology

This script systematically varies hL for different rGL values to study convergence
with increasing DOF. Results are saved in separate folders for each rGL value.

Paper parameters (Section 5.2):
- Local mesh dimensions: WL=0.5, aL=0.25, lL=0.25, HL=0.25
- Isotropic condition: crack front direction element size = crack plane element size
  This means: 2*sin(d_theta/2) = hL
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

class Verification52:
    def __init__(self, domain_scale=0.25, skip_existing=True):
        # Paper parameters for local mesh (in normalized units)
        self.WL = domain_scale * 2
        # Scale factor for local domain size (0.25 in paper, 0.5 for larger domain)
        self.domain_scale = domain_scale
        self.aL = domain_scale  # Crack radius (0.25 in paper)
        self.lL = domain_scale  # Ligament length (0.25 in paper)
        self.HL = domain_scale  # Height (0.25 in paper)
        
        # rGL values to test
        self.rGL_values = [2, 4, 8]
        
        # DOF limit
        self.max_DOF = 2e6
        
        # Skip existing results
        self.skip_existing = skip_existing
        
        # Base directory for results
        self.base_results_dir = "results/verification_5_2"
        
        # Store all run configurations
        self.run_history = []
        
    def calculate_d_theta(self, hL):
        """
        Calculate d_theta to satisfy isotropic condition
        2*sin(d_theta/2) = hL
        d_theta = 2 * arcsin(hL/2)
        
        Args:
            hL: Local element size
        
        Returns:
            d_theta in degrees
        """
        if hL > 2.0:
            raise ValueError(f"hL={hL} is too large, must be <= 2.0 for isotropic condition")
        
        d_theta_rad = 2.0 * np.arcsin(hL / 2.0)
        d_theta_deg = np.degrees(d_theta_rad)
        return d_theta_deg
    
    def calculate_local_elements(self, hL):
        """
        Calculate number of local elements based on paper dimensions
        
        Paper dimensions (normalized):
        - WL = 0.5 (crack width in local mesh)
        - aL = 0.25 (crack radius in local mesh)
        - lL = 0.25 (ligament length)
        - HL = 0.25 (height)
        
        Element counts:
        - aL_elements = aL / hL (exact division since hL = 0.25/n)
        - lL_elements = lL / hL
        - HL_elements = HL / hL
        
        Args:
            hL: Local element size (must be 0.25/n for integer n)
            
        Returns:
            Dictionary with element counts
        """
        # Since hL = 1/(4*n), we have 0.25/hL = n (exact integer)
        # Use round() to handle any floating point errors
        aL_elements = int(round(self.aL / hL))
        lL_elements = int(round(self.lL / hL))
        HL_elements = int(round(self.HL / hL))
        
        return {
            'aL': aL_elements,
            'lL': lL_elements,
            'HL': HL_elements
        }
    
    def estimate_DOF(self, hL, rGL):
        """
        Estimate degrees of freedom for a given configuration
        
        This is a rough estimation based on:
        - Local mesh elements
        - Global mesh elements (depends on rGL)
        - B-spline degree and control points
        
        Args:
            hL: Local element size
            rGL: Global to local element size ratio
            
        Returns:
            Estimated DOF count
        """
        # Local mesh elements
        local_elems = self.calculate_local_elements(hL)
        aL_el = local_elems['aL']
        lL_el = local_elems['lL']
        HL_el = local_elems['HL']
        
        # Estimate local mesh DOF (approximate)
        # This is rough - actual DOF depends on mesh topology
        local_DOF = aL_el * lL_el * HL_el * 8  # Rough estimate
        
        # Global mesh size
        hG = hL * rGL
        
        # Global domain (normalized): WidthG=2.0, HeightG=1.0
        WidthG = 2.0
        HeightG = 1.0
        
        # Global mesh control points (rough estimate)
        mu_G = 0.99 ** 0.5  # From const_global_mesh.py
        nPtsX = int(np.ceil(WidthG * mu_G / hG)) + 2  # +p for B-spline
        nPtsY = nPtsX
        nPtsZ = int(np.ceil(HeightG * mu_G / hG)) + 2
        
        global_DOF = nPtsX * nPtsY * nPtsZ * 3  # 3 DOF per node
        
        # Total DOF (rough estimate - not exact)
        total_DOF = local_DOF + global_DOF
        
        return total_DOF
    
    def generate_hL_sequence(self, rGL):
        """
        Generate sequence of hL values from coarse to fine
        Stop when estimated DOF exceeds max_DOF
        
        Strategy:
        - hL must satisfy: hL = domain_scale / n, where n is an integer
        - For domain_scale=0.5: hL = 0.5 / n = 1 / (2*n)
        - This ensures aL, lL, HL are exactly integers
        - Start with small n (coarse mesh), increase n (finer mesh)
        - Stop when DOF > max_DOF
        
        Args:
            rGL: Global to local element size ratio
            
        Returns:
            List of hL values
        """
        hL_values = []
        
        # hL = domain_scale / n, where n is the number of elements
        # For domain_scale=0.5: hL = 0.5 / n = 1 / (2*n)
        # Start with coarse mesh: n = 3 gives hL = 0.5/3 ≈ 0.1667
        n_start = 3
        
        # Maximum n to try (safety limit)
        # n = 50 gives hL = 1/200 = 0.005
        n_max = 50
        
        for n in range(n_start, n_max + 1):
            # Calculate hL = domain_scale / n to ensure exact representation
            # For domain_scale=0.5: hL = 0.5 / n = 1 / (2*n)
            # For domain_scale=0.25: hL = 0.25 / n = 1 / (4*n)
            hL = self.domain_scale / float(n)
            
            # Check if hL is valid for isotropic condition (hL <= 2.0)
            if hL > 2.0:
                continue
                
            # Estimate DOF
            est_DOF = self.estimate_DOF(hL, rGL)
            
            if est_DOF > self.max_DOF:
                break
            
            hL_values.append(hL)
        
        return hL_values
    
    def update_const_files(self, hL, rGL, d_theta, local_elems):
        """
        Update const files with current parameters
        
        Args:
            hL: Local element size
            rGL: Global to local element size ratio
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
                new_lines.append(f"rGL = {rGL}  # Ratio of global to local element size (hG/hL)\n")
            else:
                new_lines.append(line)
        
        with open(const_global_mesh_path, 'w') as f:
            f.writelines(new_lines)
    
    def create_result_folder(self, rGL, hL, hG):
        """
        Create result folder with naming convention: rGL{rGL}_{domain_scale}/hL_{hL}_hG_{hG}
        
        Args:
            rGL: Global to local element size ratio
            hL: Local element size
            hG: Global element size
            
        Returns:
            Tuple of (path to result folder, whether it already existed)
        """
        # Create rGL folder with domain_scale identifier
        rGL_folder = os.path.join(self.base_results_dir, f"rGL{rGL}_{self.domain_scale}")
        os.makedirs(rGL_folder, exist_ok=True)
        
        # Create specific run folder
        run_folder_name = f"hL_{hL:.6f}_hG_{hG:.6f}"
        run_folder = os.path.join(rGL_folder, run_folder_name)
        
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
        
        # Check log directory
        log_dir = os.path.join(result_dir, "log")
        if not os.path.exists(log_dir):
            return False
        
        # Check visual directory
        visual_dir = os.path.join(result_dir, "visual")
        if not os.path.exists(visual_dir):
            return False
        
        # Check for essential result files in log directory
        essential_log_files = [
            os.path.join(log_dir, "u.g.dat"),
            os.path.join(log_dir, "u.l.dat"),
            os.path.join(log_dir, "u_gl.l.dat")
        ]
        
        for file_path in essential_log_files:
            if not os.path.exists(file_path):
                return False
            # Check if file is not empty
            if os.path.getsize(file_path) == 0:
                return False
        
        # Check for at least one .vtu file in visual directory
        vtu_files = [f for f in os.listdir(visual_dir) if f.endswith('.vtu')]
        if len(vtu_files) == 0:
            return False
        
        return True
    
    def run_simulation(self, hL, rGL, force_run=False):
        """
        Run a single simulation with given parameters
        
        Args:
            hL: Local element size
            rGL: Global to local element size ratio
            force_run: If True, run even if results already exist
            
        Returns:
            Dictionary with run information, or None if skipped
        """
        # Calculate parameters
        d_theta = self.calculate_d_theta(hL)
        local_elems = self.calculate_local_elements(hL)
        hG = hL * rGL
        
        # Check if results already exist
        run_folder, already_exists = self.create_result_folder(rGL, hL, hG)
        
        if already_exists and self.skip_existing and not force_run:
            print(f"\n{'='*70}")
            print(f"SKIPPING (already exists): rGL={rGL}, hL={hL:.6f}, hG={hG:.6f}")
            print(f"{'='*70}")
            print(f"  Result folder: {run_folder}")
            print(f"  Use --force to rerun existing simulations")
            
            # Return existing configuration info
            run_info = {
                'rGL': rGL,
                'hL': hL,
                'hG': hG,
                'd_theta': d_theta,
                'aL': local_elems['aL'],
                'lL': local_elems['lL'],
                'HL': local_elems['HL'],
                'step': int(round(1.0 / hL)),
                'estimated_DOF': self.estimate_DOF(hL, rGL),
                'result_folder': run_folder,
                'skipped': True
            }
            return run_info
        
        print(f"\n{'='*70}")
        print(f"Running simulation: rGL={rGL}, hL={hL:.6f}, hG={hG:.6f}")
        print(f"{'='*70}")
        print(f"Local elements: aL={local_elems['aL']}, lL={local_elems['lL']}, HL={local_elems['HL']}")
        print(f"d_theta: {d_theta:.6f} degrees")
        print(f"Estimated DOF: {self.estimate_DOF(hL, rGL):.2e}")
        
        # Update configuration files
        self.update_const_files(hL, rGL, d_theta, local_elems)
        
        # Run main.py in static mode
        cmd = ["python3", "main.py", "--static_only"]
        
        print(f"\nExecuting: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"ERROR: Simulation failed!")
            print(f"STDOUT:\n{result.stdout}")
            print(f"STDERR:\n{result.stderr}")
            return None
        
        print(f"Simulation completed successfully")
        
        # Determine step number (should be auto-calculated by main.py)
        # step = round(c / hL) where c=1.0 in static mode
        c = 1.0
        step = int(round(c / hL))
        step_str = f"{step:05d}"
        
        # Copy result files
        source_result_dir = f"results/step{step_str}"
        if os.path.exists(source_result_dir):
            # Copy entire result directory
            dest_result_dir = os.path.join(run_folder, f"step{step_str}")
            if os.path.exists(dest_result_dir):
                shutil.rmtree(dest_result_dir)
            shutil.copytree(source_result_dir, dest_result_dir)
            print(f"Results copied to: {run_folder}")
        else:
            print(f"WARNING: Result directory not found: {source_result_dir}")
        
        # Copy input files as well
        source_input_dir = f"inputfiles/step{step_str}"
        if os.path.exists(source_input_dir):
            dest_input_dir = os.path.join(run_folder, f"inputfiles_step{step_str}")
            if os.path.exists(dest_input_dir):
                shutil.rmtree(dest_input_dir)
            shutil.copytree(source_input_dir, dest_input_dir)
        
        # Save run configuration
        run_info = {
            'rGL': rGL,
            'hL': hL,
            'hG': hG,
            'd_theta': d_theta,
            'aL': local_elems['aL'],
            'lL': local_elems['lL'],
            'HL': local_elems['HL'],
            'step': step,
            'estimated_DOF': self.estimate_DOF(hL, rGL),
            'timestamp': datetime.now().isoformat(),
            'result_folder': run_folder,
            'skipped': False
        }
        
        # Save configuration to JSON
        config_file = os.path.join(run_folder, "run_config.json")
        with open(config_file, 'w') as f:
            json.dump(run_info, f, indent=2)
        
        return run_info
    
    def run_all(self):
        """
        Run all simulations for all rGL values
        """
        print("="*70)
        print("Verification 5.2: Isotropic Local Mesh Convergence Study")
        print("="*70)
        print(f"Base results directory: {self.base_results_dir}")
        print(f"rGL values to test: {self.rGL_values}")
        print(f"Maximum DOF: {self.max_DOF:.2e}")
        print(f"\nLocal domain parameters:")
        print(f"  Domain scale: {self.domain_scale}")
        print(f"  WL = {self.WL}")
        print(f"  aL = {self.aL}")
        print(f"  lL = {self.lL}")
        print(f"  HL = {self.HL}")
        print(f"  Isotropic condition: 2*sin(d_theta/2) = hL")
        
        # Create base results directory
        os.makedirs(self.base_results_dir, exist_ok=True)
        
        # Run for each rGL value
        for rGL in self.rGL_values:
            print(f"\n\n{'#'*70}")
            print(f"# Processing rGL = {rGL}")
            print(f"{'#'*70}")
            
            # Generate hL sequence
            hL_values = self.generate_hL_sequence(rGL)
            print(f"\nGenerated {len(hL_values)} hL values for rGL={rGL}:")
            for i, hL in enumerate(hL_values):
                est_DOF = self.estimate_DOF(hL, rGL)
                print(f"  {i+1}. hL={hL:.6f}, hG={hL*rGL:.6f}, est. DOF={est_DOF:.2e}")
            
            # Run simulations
            for i, hL in enumerate(hL_values):
                print(f"\n--- Simulation {i+1}/{len(hL_values)} for rGL={rGL} ---")
                run_info = self.run_simulation(hL, rGL)
                
                if run_info:
                    self.run_history.append(run_info)
        
        # Save complete run history
        history_file = os.path.join(self.base_results_dir, "run_history.json")
        with open(history_file, 'w') as f:
            json.dump(self.run_history, f, indent=2)
        
        print(f"\n\n{'='*70}")
        print("All simulations completed!")
        print(f"{'='*70}")
        
        # Count skipped vs executed
        total_runs = len(self.run_history)
        skipped_runs = sum(1 for r in self.run_history if r.get('skipped', False))
        executed_runs = total_runs - skipped_runs
        
        print(f"Total configurations: {total_runs}")
        print(f"  Executed: {executed_runs}")
        print(f"  Skipped (already exist): {skipped_runs}")
        print(f"Results saved in: {self.base_results_dir}")
        print(f"Run history saved in: {history_file}")
        
        # Print summary
        print(f"\nSummary by rGL:")
        for rGL in self.rGL_values:
            rGL_runs = [r for r in self.run_history if r['rGL'] == rGL]
            if rGL_runs:
                rGL_skipped = sum(1 for r in rGL_runs if r.get('skipped', False))
                rGL_executed = len(rGL_runs) - rGL_skipped
                print(f"  rGL={rGL}: {len(rGL_runs)} total ({rGL_executed} executed, {rGL_skipped} skipped)")
                print(f"    hL range: {min(r['hL'] for r in rGL_runs):.6f} to {max(r['hL'] for r in rGL_runs):.6f}")
                print(f"    DOF range: {min(r['estimated_DOF'] for r in rGL_runs):.2e} to {max(r['estimated_DOF'] for r in rGL_runs):.2e}")


def main():
    """Main entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Verification script for Section 5.2 of the paper",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--rGL",
        type=int,
        nargs='+',
        default=[2, 4, 8],
        help="List of rGL values to test (default: 2 4 8)"
    )
    parser.add_argument(
        "--max_DOF",
        type=float,
        default=2e6,
        help="Maximum DOF limit (default: 2e6)"
    )
    parser.add_argument(
        "--domain_scale",
        type=float,
        default=0.25,
        help="Local domain scale factor (aL=lL=HL=domain_scale, default: 0.25)"
    )
    parser.add_argument(
        "--no-skip",
        action="store_true",
        help="Do not skip existing results, rerun all simulations"
    )
    
    args = parser.parse_args()
    
    # Create verification object with specified domain scale
    verif = Verification52(
        domain_scale=args.domain_scale,
        skip_existing=not args.no_skip
    )
    verif.rGL_values = args.rGL
    verif.max_DOF = args.max_DOF
    
    # Run all simulations
    verif.run_all()


if __name__ == "__main__":
    main()
