#!/usr/bin/env python3
"""
Verification script for Section 5.3: Verification of the proposed strategy 
using local mesh coarsened in the crack front direction

Based on the paper's Figure 16 and 17 methodology

This script studies the effect of the crack-front-direction element size on
the solution accuracy. Section 5.3 varies the target ratio

    crack-front-direction element size / crack-plane element size = x

using the same theta rule as Section 5.2:

    2*sin(d_theta/2) = x*hL

Paper parameters (Section 5.3):
- Local mesh dimensions: WL=0.5, aL=0.25, lL=0.25, HL=0.25
- Fixed nominal rGL = 4
- Two mesh configurations:
  a) hL = 1/48
  b) hL = 1/88
- Global control point overrides:
  a) (27, 27, 15)
  b) (47, 47, 24)
- Target crack-front/crack-plane element-size ratios to test: 1.5, 2, 2.5, 3, 4

Total: 5 target ratios × 2 mesh configs = 10 cases
"""

import math
import os
import sys
import subprocess
import shutil
import json
from datetime import datetime

# Add script directory to path so "const" imports work from any launch cwd.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from const import const_global_mesh as cgm


class Verification53:
    def __init__(self, domain_scale=0.25, skip_existing=True):
        self.script_dir = os.path.dirname(os.path.abspath(__file__))

        # Paper parameters for local mesh (in normalized units)
        self.WL = domain_scale * 2
        self.domain_scale = domain_scale
        self.aL = domain_scale  # Crack radius (0.25 in paper)
        self.lL = domain_scale  # Ligament length (0.25 in paper)
        self.HL = domain_scale  # Height (0.25 in paper)
        
        # Fixed nominal rGL value for section 5.3. The actual hG is taken from
        # the explicit static global nPts overrides below.
        self.rGL = 4
        
        # Fixed rBL value
        self.rBL = 0.25
        
        # Static crack/domain constants used by main.py in --static_only mode.
        self.static_c = 1.0
        self.static_widthG = 2.0
        self.static_heightG = 1.0

        # Two mesh configurations
        self.hL_values = [1.0/48.0, 1.0/88.0]

        # Target ratios x for 2*sin(d_theta/2) = x*hL.
        self.target_ratio_values = [1.5, 2.0, 2.5, 3.0, 4.0]

        # Section 5.3 special global mesh settings (control points).
        # The user's requested pairs are interpreted as:
        #   hL=1/48  -> (nPtsX, nPtsZ) = (27, 15)
        #   hL=1/88  -> (nPtsX, nPtsZ) = (47, 24)
        # and we keep XY symmetry by setting nPtsY = nPtsX.
        self.global_npts_overrides = {
            48: (27, 27, 15),  # (nPtsX, nPtsY, nPtsZ)
            88: (47, 47, 24)
        }
        
        # Skip existing results
        self.skip_existing = skip_existing
        
        # Base directory for results
        self.base_results_dir = os.path.join(self.script_dir, "results/verification_5_3")
        
        # Store all run configurations
        self.run_history = []

    def get_global_npts_override(self, hL):
        """
        Get global control point counts (nPtsX, nPtsY, nPtsZ) for a given hL.
        Section 5.3 uses explicit overrides only.
        """
        step_count = int(round(1.0 / hL))
        if step_count not in self.global_npts_overrides:
            raise ValueError(
                f"No global_npts_overrides entry for hL={hL} (1/{step_count}). "
                f"Please add step_count={step_count} to self.global_npts_overrides."
            )
        return self.global_npts_overrides[step_count]

    def get_static_domain_lengths(self):
        """
        Return adjusted static global domain lengths used by main.py.

        In static mode, main.py sets WidthG=2.0 and HeightG=1.0, then applies
        cgm.mu_G to get Lx, Ly, and Lz.
        """
        return (
            self.static_widthG * cgm.mu_G,
            self.static_widthG * cgm.mu_G,
            self.static_heightG * cgm.mu_G,
        )

    def calculate_global_mesh_sizes(self, global_npts):
        """
        Calculate actual global element sizes from static nPts overrides.

        For open quadratic B-splines, the number of elements in X is nPtsX-p,
        so the X-direction element size is Lx/(nPtsX-p).
        """
        nPtsX, nPtsY, nPtsZ = global_npts
        Lx, Ly, Lz = self.get_static_domain_lengths()

        if nPtsX <= cgm.p or nPtsY <= cgm.q or nPtsZ <= cgm.r:
            raise ValueError(
                "Global nPts must be larger than the spline degree in every direction."
            )

        return {
            'hGx': Lx / (nPtsX - cgm.p),
            'hGy': Ly / (nPtsY - cgm.q),
            'hGz': Lz / (nPtsZ - cgm.r),
            'Lx': Lx,
            'Ly': Ly,
            'Lz': Lz,
        }

    def calculate_theta_from_ratio(self, target_ratio, hL):
        """
        Calculate d_theta from the Section 5.2-style crack-front ratio rule.

        2*sin(d_theta/2) = target_ratio*hL

        local_mesh.py converts the input d_theta to nLtheta by round(90/d_theta)
        and then spaces the actual local mesh uniformly over 90 degrees, so this
        method records both the exact input d_theta and the rounded actual value.
        """
        if target_ratio <= 0:
            raise ValueError(f"target_ratio must be positive, got {target_ratio}")

        target_arg = target_ratio * hL / 2.0
        if target_arg > 1.0:
            raise ValueError(
                f"target_ratio={target_ratio} is too large: "
                f"target_ratio*hL/2={target_arg:.6f} > 1."
            )

        target_theta_rad = 2.0 * math.asin(target_arg)
        target_theta_deg = math.degrees(target_theta_rad)
        if target_theta_deg <= 0.0:
            raise ValueError(
                f"target_ratio={target_ratio} produced non-positive target theta."
            )

        nLtheta = max(1, int(round(90.0 / target_theta_deg)))
        actual_d_theta = 90.0 / nLtheta
        theta_element_size = 2.0 * math.sin(math.radians(actual_d_theta) / 2.0)
        actual_ratio = theta_element_size / hL
        ratio_error = actual_ratio - target_ratio

        return {
            'target_ratio': target_ratio,
            'target_d_theta': target_theta_deg,
            'd_theta': target_theta_deg,
            'nLtheta': nLtheta,
            'actual_d_theta': actual_d_theta,
            'theta_element_size': theta_element_size,
            'actual_ratio': actual_ratio,
            'ratio_error': ratio_error,
            'relative_ratio_error': ratio_error / target_ratio,
        }

    def build_case_parameters(self, hL, target_ratio):
        """
        Build all derived parameters for a single Section 5.3 run.
        """
        local_elems = self.calculate_local_elements(hL)
        global_npts = self.get_global_npts_override(hL)
        global_mesh_sizes = self.calculate_global_mesh_sizes(global_npts)
        hG = global_mesh_sizes['hGx']
        theta_info = self.calculate_theta_from_ratio(target_ratio, hL)

        return {
            'hL': hL,
            'hG': hG,
            'actual_rGL': hG / hL,
            'target_ratio': target_ratio,
            'local_elems': local_elems,
            'global_npts': global_npts,
            'global_mesh_sizes': global_mesh_sizes,
            'estimated_DOF': self.estimate_DOF(hL, theta_info['nLtheta']),
            **theta_info,
        }
        
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
    
    def estimate_DOF(self, hL, nLtheta=None):
        """
        Estimate degrees of freedom for a given configuration
        
        Args:
            hL: Local element size
            nLtheta: Number of local elements in theta direction
            
        Returns:
            Estimated DOF count
        """
        # Local mesh elements
        local_elems = self.calculate_local_elements(hL)
        aL_el = local_elems['aL']
        lL_el = local_elems['lL']
        HL_el = local_elems['HL']
        
        # Local mesh nodes: (nLr+1) * (nLtheta+1) * (HL+1), 3 DOF per node.
        # If nLtheta is not provided, fall back to the old coarse estimate.
        if nLtheta is None:
            local_DOF = aL_el * lL_el * HL_el * 8
        else:
            nLr = aL_el + lL_el
            local_DOF = (nLr + 1) * (nLtheta + 1) * (HL_el + 1) * 3
        
        # Global mesh control points
        nPtsX, nPtsY, nPtsZ = self.get_global_npts_override(hL)
        
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
        const_local_mesh_path = os.path.join(self.script_dir, "const/const_local_mesh.py")
        
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
                new_lines.append(f"d_theta = {d_theta:.12f}  # Angular resolution [degrees]\n")
            else:
                new_lines.append(line)
        
        with open(const_local_mesh_path, 'w') as f:
            f.writelines(new_lines)
        
        # Update const_global_mesh.py
        const_global_mesh_path = os.path.join(self.script_dir, "const/const_global_mesh.py")
        
        with open(const_global_mesh_path, 'r') as f:
            lines = f.readlines()
        
        new_lines = []
        for line in lines:
            if re.match(r'^rGL\s*=', line):
                new_lines.append(f"rGL = {self.rGL}  # Nominal global/local element-size ratio\n")
            elif re.match(r'^rBL\s*=', line):
                new_lines.append(f"rBL = {self.rBL}  # Ratio of background to local element size (hB/hL)\n")
            else:
                new_lines.append(line)
        
        with open(const_global_mesh_path, 'w') as f:
            f.writelines(new_lines)
    
    @staticmethod
    def format_float_for_path(value, decimals=6):
        """
        Format a float compactly for result folder names.
        """
        text = f"{value:.{decimals}f}".rstrip('0').rstrip('.')
        return text if text else "0"

    def create_result_folder(self, hL, target_ratio, d_theta):
        """
        Create result folder with naming convention: 
        verification_5_3/hL_{hL}/ratio_{target_ratio}_dtheta_{d_theta}
        
        Args:
            hL: Local element size
            target_ratio: Requested crack-front/crack-plane element-size ratio
            d_theta: Angular resolution (degrees)
            
        Returns:
            Tuple of (path to result folder, whether it already existed)
        """
        # Create hL folder
        hL_folder = os.path.join(self.base_results_dir, f"hL_{hL:.8f}")
        os.makedirs(hL_folder, exist_ok=True)
        
        # Create specific target-ratio folder
        ratio_text = self.format_float_for_path(target_ratio, decimals=3)
        dtheta_text = self.format_float_for_path(d_theta, decimals=6)
        run_folder_name = f"ratio_{ratio_text}_dtheta_{dtheta_text}"
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
        step = int(round(self.static_c / hL))
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
    
    def save_run_config(self, run_folder, case, status):
        """
        Save run configuration to JSON file
        
        Args:
            run_folder: Path to result folder
            case: Dictionary with derived run parameters
            status: Case status string
        """
        global_npts = case['global_npts']
        global_mesh_sizes = case['global_mesh_sizes']

        config = {
            'timestamp': datetime.now().isoformat(),
            'section': '5.3',
            'status': status,
            'description': 'Verification using target crack-front/crack-plane element-size ratios',
            'theta_rule': '2*sin(d_theta/2) = target_ratio*hL',
            'domain_scale': self.domain_scale,
            'hL': case['hL'],
            'hG': case['hG'],
            'rGL': self.rGL,
            'actual_rGL': case['actual_rGL'],
            'rBL': self.rBL,
            'target_ratio': case['target_ratio'],
            'target_d_theta': case['target_d_theta'],
            'd_theta': case['d_theta'],
            'nLtheta': case['nLtheta'],
            'actual_d_theta': case['actual_d_theta'],
            'theta_element_size': case['theta_element_size'],
            'actual_ratio': case['actual_ratio'],
            'ratio_error': case['ratio_error'],
            'relative_ratio_error': case['relative_ratio_error'],
            'local_elements': case['local_elems'],
            'global_control_points': {
                'nPtsX': global_npts[0],
                'nPtsY': global_npts[1],
                'nPtsZ': global_npts[2]
            },
            'global_mesh_sizes': global_mesh_sizes,
            'estimated_DOF': case['estimated_DOF'],
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
    
    def run_simulation(self, run_folder, case):
        """
        Run the simulation
        
        Args:
            run_folder: Path to result folder
            case: Dictionary with derived run parameters
            
        Returns:
            True if successful, False otherwise
        """
        print(f"\n  Running simulation...")

        hL = case['hL']
        global_npts = case['global_npts']
        
        # Calculate step number
        step = int(round(self.static_c / hL))
        
        # Run main.py with static_only mode and section 5.3 nPts overrides
        cmd = [
            'python3', 'main.py', '--static_only',
            '--static_nptsx', str(global_npts[0]),
            '--static_nptsy', str(global_npts[1]),
            '--static_nptsz', str(global_npts[2])
        ]
        
        try:
            # Run simulation
            print(f"  Command: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                cwd=self.script_dir
            )
            
            print(f"  Simulation completed successfully")
            
            # Copy results to result folder
            step_str = f"{step:05d}"
            inputfiles_dir = os.path.join(self.script_dir, f"inputfiles/step{step_str}")
            results_dir = os.path.join(self.script_dir, f"results/step{step_str}")
            
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
    
    def run_single_case(self, hL, target_ratio):
        """
        Run a single case with given hL and target crack-front/crack-plane ratio
        
        Args:
            hL: Local element size
            target_ratio: Requested crack-front/crack-plane element-size ratio
            
        Returns:
            Case status: "completed", "skipped", or "failed"
        """
        # Calculate parameters
        case = self.build_case_parameters(hL, target_ratio)
        hG = case['hG']
        d_theta = case['d_theta']
        local_elems = case['local_elems']
        global_npts = case['global_npts']
        est_DOF = case['estimated_DOF']
        
        print(f"\n{'='*80}")
        print(
            f"Running case: hL={hL:.8f} (1/{1/hL:.0f}), "
            f"target ratio={target_ratio:.3f}, d_theta={d_theta:.6f}°"
        )
        print(f"{'='*80}")
        print(f"  nominal rGL = {self.rGL}")
        print(f"  actual hG/hL = {case['actual_rGL']:.6f}")
        print(f"  hL = {hL:.8f}")
        print(f"  hG = Lx/(nPtsX-{cgm.p}) = {hG:.8f}")
        print(f"  target theta/hL ratio = {target_ratio:.6f}")
        print(f"  actual theta/hL ratio = {case['actual_ratio']:.6f}")
        print(f"  theta element size = {case['theta_element_size']:.8f}")
        print(f"  d_theta input = {d_theta:.6f}°")
        print(f"  nLtheta = {case['nLtheta']}, actual d_theta = {case['actual_d_theta']:.6f}°")
        print(f"  Global control points: nPtsX={global_npts[0]}, nPtsY={global_npts[1]}, nPtsZ={global_npts[2]}")
        print(f"  Local elements: aL={local_elems['aL']}, lL={local_elems['lL']}, HL={local_elems['HL']}")
        print(f"  Estimated DOF: {est_DOF:,}")
        
        # Create result folder
        run_folder, already_exists = self.create_result_folder(hL, target_ratio, d_theta)
        print(f"  Result folder: {run_folder}")
        
        if already_exists and self.skip_existing:
            print(f"  SKIPPED: Results already exist")
            # Still save config for record keeping
            self.save_run_config(run_folder, case, status="skipped")
            return "skipped"
        
        # Update const files
        print(f"  Updating const files...")
        self.update_const_files(hL, d_theta, local_elems)
        
        # Run simulation
        success = self.run_simulation(run_folder, case)
        
        if success:
            print(f"\n  ✓ Case completed successfully")
            self.save_run_config(run_folder, case, status="completed")
            return "completed"
        else:
            print(f"\n  ✗ Case failed")
            self.save_run_config(run_folder, case, status="failed")
            return "failed"
    
    def run_all_cases(self):
        """
        Run all cases for section 5.3
        
        Total cases: target ratio values × hL values
        """
        print(f"\n{'='*80}")
        print(f"Section 5.3 Verification")
        print(f"Study of crack-front/crack-plane element-size ratio effect on solution accuracy")
        print(f"{'='*80}")
        print(f"\nConfiguration:")
        print(f"  nominal rGL = {self.rGL}")
        print(f"  rBL = {self.rBL}")
        print(f"  hL values: {self.hL_values}")
        print(f"  target theta/hL ratios: {self.target_ratio_values}")
        print(f"  Total cases: {len(self.hL_values) * len(self.target_ratio_values)}")
        print(f"  Skip existing: {self.skip_existing}")
        print(f"  Results directory: {self.base_results_dir}")
        
        # Create base results directory
        os.makedirs(self.base_results_dir, exist_ok=True)
        
        # Track statistics
        total_cases = len(self.hL_values) * len(self.target_ratio_values)
        completed = 0
        failed = 0
        skipped = 0
        
        # Run all cases
        case_num = 0
        for hL in self.hL_values:
            for target_ratio in self.target_ratio_values:
                case_num += 1
                print(f"\n{'='*80}")
                print(f"Case {case_num}/{total_cases}")
                print(f"{'='*80}")
                
                status = self.run_single_case(hL, target_ratio)

                if status == "completed":
                    completed += 1
                elif status == "skipped":
                    skipped += 1
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
        choices=['48', '88', 'all'],
        default='all',
        help='Which hL config to run: 48 (1/48), 88 (1/88), or all (default: all)'
    )
    parser.add_argument(
        '--ratios',
        '--ratio',
        dest='ratios',
        type=float,
        nargs='+',
        help='Run only specific target theta/hL ratio(s), e.g. --ratios 1.5 2.5'
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
        elif args.hL == '88':
            verif.hL_values = [1.0/88.0]
    
    # Filter target ratios if specified
    if args.ratios is not None:
        verif.target_ratio_values = args.ratios
    
    # Run all cases
    verif.run_all_cases()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
