#!/usr/bin/env python3
"""
Run the Section 5.2 plotting case with manually specified mesh density.

This case keeps the static global domain size used by main.py:

    WidthG x WidthG x HeightG = 2 x 2 x 1

and uses a larger local window:

    WL = 1.0, aL = 0.5, lL = 0.5, HL = 0.5, hL = 1/8

The global IGA mesh is specified directly by control point counts:

    nPtsX = 12, nPtsY = 12, nPtsZ = 7

With the current quadratic basis, this gives 10 x 10 x 5 IGA elements.
"""

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

# Add script directory to path so "const" imports work from any launch cwd.
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from const import const_global_mesh as cgm


class PlotCase52:
    def __init__(self, skip_existing=True, keep_const=False, meshonly=False):
        self.script_dir = SCRIPT_DIR
        self.skip_existing = skip_existing
        self.keep_const = keep_const
        self.meshonly = meshonly

        self.static_c = 1.0
        self.static_widthG = 2.0
        self.static_heightG = 1.0

        self.hL = 1.0 / 8.0
        self.local_lengths = {
            "WL": 1.0,
            "aL": 0.5,
            "lL": 0.5,
            "HL": 0.5,
        }
        self.global_npts = (12, 12, 7)
        self.expected_global_elements = (10, 10, 5)

        self.base_results_dir = self.script_dir / "results" / "verification_5_2" / "plot_case"
        self.const_local_mesh_path = self.script_dir / "const" / "const_local_mesh.py"

    def calculate_d_theta(self):
        """Return the input d_theta satisfying 2*sin(d_theta/2) = hL."""
        return math.degrees(2.0 * math.asin(self.hL / 2.0))

    def calculate_local_elements(self):
        elems = {}
        for key in ("aL", "lL", "HL"):
            count = self.local_lengths[key] / self.hL
            rounded = int(round(count))
            if not math.isclose(count, rounded, rel_tol=0.0, abs_tol=1.0e-12):
                raise ValueError(f"{key}={self.local_lengths[key]} is not divisible by hL={self.hL}")
            elems[key] = rounded

        radial_width = (elems["aL"] + elems["lL"]) * self.hL
        if not math.isclose(radial_width, self.local_lengths["WL"], rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError(
                f"WL mismatch: (aL+lL)*hL={radial_width} but WL={self.local_lengths['WL']}"
            )

        return elems

    def calculate_global_mesh_info(self):
        nPtsX, nPtsY, nPtsZ = self.global_npts
        elements = (nPtsX - cgm.p, nPtsY - cgm.q, nPtsZ - cgm.r)

        if elements != self.expected_global_elements:
            raise ValueError(
                f"Global element mismatch: got {elements}, expected {self.expected_global_elements}. "
                f"Check nPts and B-spline degrees."
            )

        Lx = self.static_widthG * cgm.mu_G
        Ly = self.static_widthG * cgm.mu_G
        Lz = self.static_heightG * cgm.mu_G

        return {
            "global_domain": {
                "WidthG": self.static_widthG,
                "HeightG": self.static_heightG,
                "adjusted_Lx": Lx,
                "adjusted_Ly": Ly,
                "adjusted_Lz": Lz,
                "mu_G": cgm.mu_G,
            },
            "basis_degree": {
                "p": cgm.p,
                "q": cgm.q,
                "r": cgm.r,
            },
            "control_points": {
                "nPtsX": nPtsX,
                "nPtsY": nPtsY,
                "nPtsZ": nPtsZ,
            },
            "elements": {
                "x": elements[0],
                "y": elements[1],
                "z": elements[2],
            },
            "element_sizes": {
                "hGx": Lx / elements[0],
                "hGy": Ly / elements[1],
                "hGz": Lz / elements[2],
            },
        }

    def build_case(self):
        d_theta = self.calculate_d_theta()
        nLtheta = int(round(90.0 / d_theta))
        actual_d_theta = 90.0 / nLtheta
        actual_theta_size = 2.0 * math.sin(math.radians(actual_d_theta) / 2.0)
        local_elems = self.calculate_local_elements()
        global_mesh = self.calculate_global_mesh_info()
        step = int(round(self.static_c / self.hL))

        nLr = local_elems["aL"] + local_elems["lL"]
        local_nodes = (nLr + 1) * (nLtheta + 1) * (local_elems["HL"] + 1)
        global_nodes = (
            self.global_npts[0]
            * self.global_npts[1]
            * self.global_npts[2]
        )

        return {
            "timestamp": datetime.now().isoformat(),
            "section": "5.2",
            "case_name": "plot_case_manual_global_10x10x5_local_hL_1_8",
            "description": "Plotting case based on verification 5.2 with manually specified mesh density.",
            "static_c": self.static_c,
            "hL": self.hL,
            "hL_text": "1/8",
            "step": step,
            "theta_rule": "2*sin(d_theta/2) = hL",
            "d_theta": d_theta,
            "nLtheta": nLtheta,
            "actual_d_theta": actual_d_theta,
            "actual_theta_element_size": actual_theta_size,
            "actual_theta_element_size_over_hL": actual_theta_size / self.hL,
            "local_dimensions": self.local_lengths,
            "local_elements": local_elems,
            "global_mesh": global_mesh,
            "estimated_nodes": {
                "global": global_nodes,
                "local": local_nodes,
                "total": global_nodes + local_nodes,
            },
            "estimated_DOF": 3 * (global_nodes + local_nodes),
            "meshonly": self.meshonly,
        }

    @staticmethod
    def format_float_for_path(value, decimals=8):
        text = f"{value:.{decimals}f}".rstrip("0").rstrip(".")
        return text if text else "0"

    def result_folder(self, case):
        npts = case["global_mesh"]["control_points"]
        elems = case["global_mesh"]["elements"]
        hL_text = self.format_float_for_path(case["hL"])
        folder = (
            f"hL_{hL_text}_globalElem_{elems['x']}_{elems['y']}_{elems['z']}_"
            f"nPts_{npts['nPtsX']}_{npts['nPtsY']}_{npts['nPtsZ']}"
        )
        return self.base_results_dir / folder

    def check_result_exists(self, run_folder, case):
        step_str = f"{case['step']:05d}"
        result_dir = run_folder / f"step{step_str}"
        if not result_dir.exists():
            return False

        required_files = [
            result_dir / "node.g.dat",
            result_dir / "node.l.dat",
            result_dir / "elem.g.dat",
            result_dir / "elem.l.dat",
        ]

        if not self.meshonly:
            required_files.extend(
                [
                    result_dir / "log" / "u.g.dat",
                    result_dir / "log" / "u.l.dat",
                    result_dir / "log" / "u_gl.l.dat",
                ]
            )

        return all(path.exists() and path.stat().st_size > 0 for path in required_files)

    def update_const_files(self, case):
        local_elems = case["local_elements"]
        replacements = {
            r"^hL_static\s*=": f"hL_static = {case['hL']}  # Normalized element size (dimensionless)\n",
            r"^aL_static\s*=": f"aL_static = {local_elems['aL']}\n",
            r"^lL_static\s*=": f"lL_static = {local_elems['lL']}\n",
            r"^HL_static\s*=": f"HL_static = {local_elems['HL']}\n",
            r"^d_theta\s*=": f"d_theta = {case['d_theta']:.12f}  # Angular resolution [degrees]\n",
        }

        lines = self.const_local_mesh_path.read_text().splitlines(keepends=True)
        new_lines = []
        for line in lines:
            for pattern, replacement in replacements.items():
                if re.match(pattern, line):
                    new_lines.append(replacement)
                    break
            else:
                new_lines.append(line)

        self.const_local_mesh_path.write_text("".join(new_lines))

    def copy_case_outputs(self, run_folder, case):
        step_str = f"{case['step']:05d}"
        source_result_dir = self.script_dir / "results" / f"step{step_str}"
        source_input_dir = self.script_dir / "inputfiles" / f"step{step_str}"
        dest_result_dir = run_folder / f"step{step_str}"
        dest_input_dir = run_folder / f"inputfiles_step{step_str}"

        if source_result_dir.exists():
            if dest_result_dir.exists():
                shutil.rmtree(dest_result_dir)
            shutil.copytree(source_result_dir, dest_result_dir)
            print(f"  Copied results to: {dest_result_dir}")
        elif self.meshonly and source_input_dir.exists():
            if dest_result_dir.exists():
                shutil.rmtree(dest_result_dir)
            shutil.copytree(source_input_dir, dest_result_dir)
            print(f"  Copied mesh-only files to: {dest_result_dir}")
        else:
            print(f"  WARNING: result directory not found: {source_result_dir}")

        if source_input_dir.exists():
            if dest_input_dir.exists():
                shutil.rmtree(dest_input_dir)
            shutil.copytree(source_input_dir, dest_input_dir)
            print(f"  Copied input files to: {dest_input_dir}")
        else:
            print(f"  WARNING: input directory not found: {source_input_dir}")

    def save_config(self, run_folder, case, status):
        config = dict(case)
        config["status"] = status
        config["result_folder"] = str(run_folder)
        config["completed_at"] = datetime.now().isoformat()

        run_folder.mkdir(parents=True, exist_ok=True)
        with (run_folder / "run_config.json").open("w") as f:
            json.dump(config, f, indent=2)

        history_file = self.base_results_dir / "run_history.json"
        if history_file.exists():
            with history_file.open("r") as f:
                history = json.load(f)
        else:
            history = []
        history.append(config)
        with history_file.open("w") as f:
            json.dump(history, f, indent=2)

    def print_case_summary(self, case, run_folder):
        npts = case["global_mesh"]["control_points"]
        elems = case["global_mesh"]["elements"]
        hgs = case["global_mesh"]["element_sizes"]
        local = case["local_elements"]

        print("=" * 78)
        print("Verification 5.2 plotting case")
        print("=" * 78)
        print(f"Result folder: {run_folder}")
        print(f"hL = {case['hL']:.8f} ({case['hL_text']}), step = {case['step']}")
        print(f"Local lengths: WL=1.0, aL=0.5, lL=0.5, HL=0.5")
        print(f"Local elements: aL={local['aL']}, lL={local['lL']}, HL={local['HL']}")
        print(f"d_theta = {case['d_theta']:.12f} deg from 2*sin(d_theta/2)=hL")
        print(f"nLtheta = {case['nLtheta']}, actual d_theta in mesh = {case['actual_d_theta']:.12f} deg")
        print(
            "Global elements/control points: "
            f"{elems['x']}x{elems['y']}x{elems['z']} / "
            f"{npts['nPtsX']}x{npts['nPtsY']}x{npts['nPtsZ']}"
        )
        print(f"Global hG: x={hgs['hGx']:.12f}, y={hgs['hGy']:.12f}, z={hgs['hGz']:.12f}")
        print(f"Estimated DOF: {case['estimated_DOF']}")
        print("=" * 78)

    def run(self, force=False):
        case = self.build_case()
        run_folder = self.result_folder(case)
        run_folder.mkdir(parents=True, exist_ok=True)
        self.print_case_summary(case, run_folder)

        if self.check_result_exists(run_folder, case) and self.skip_existing and not force:
            print("Existing completed result found; skipping. Use --force to rerun.")
            self.save_config(run_folder, case, status="skipped")
            return run_folder

        step_str = f"{case['step']:05d}"
        transient_dirs = [
            ("results", self.script_dir / "results" / f"step{step_str}"),
            ("inputfiles", self.script_dir / "inputfiles" / f"step{step_str}"),
        ]
        backup_dir = Path(tempfile.mkdtemp(prefix="verification_5_2_plot_case_"))
        original_const = self.const_local_mesh_path.read_text()

        try:
            for label, path in transient_dirs:
                if path.exists():
                    shutil.copytree(path, backup_dir / f"{label}_{path.name}")

            print("Updating static local mesh constants...")
            self.update_const_files(case)

            npts = case["global_mesh"]["control_points"]
            cmd = [
                "python3",
                "main.py",
                "--static_only",
                "--static_nptsx",
                str(npts["nPtsX"]),
                "--static_nptsy",
                str(npts["nPtsY"]),
                "--static_nptsz",
                str(npts["nPtsZ"]),
            ]
            if self.meshonly:
                cmd.append("--meshonly")

            print(f"Executing: {' '.join(cmd)}")
            sys.stdout.flush()
            subprocess.run(cmd, cwd=self.script_dir, check=True)

            self.copy_case_outputs(run_folder, case)
            self.save_config(run_folder, case, status="completed")
            print("Plotting case completed successfully.")
            return run_folder
        except subprocess.CalledProcessError as exc:
            self.save_config(run_folder, case, status="failed")
            raise SystemExit(f"Simulation failed with return code {exc.returncode}") from exc
        finally:
            if not self.keep_const:
                self.const_local_mesh_path.write_text(original_const)

            for label, path in transient_dirs:
                if path.exists():
                    shutil.rmtree(path)
                backup_path = backup_dir / f"{label}_{path.name}"
                if backup_path.exists():
                    shutil.copytree(backup_path, path)

            shutil.rmtree(backup_dir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(
        description="Run the manually specified Section 5.2 plotting case."
    )
    parser.add_argument("--force", action="store_true", help="Rerun even if archived results already exist.")
    parser.add_argument("--no-skip", action="store_true", help="Alias for --force.")
    parser.add_argument("--keep-const", action="store_true", help="Leave const_local_mesh.py set to this case after running.")
    parser.add_argument("--meshonly", action="store_true", help="Generate mesh/input files only; skip the solver.")
    args = parser.parse_args()

    runner = PlotCase52(
        skip_existing=not (args.force or args.no_skip),
        keep_const=args.keep_const,
        meshonly=args.meshonly,
    )
    run_folder = runner.run(force=args.force or args.no_skip)
    print(f"Archived case folder: {run_folder}")


if __name__ == "__main__":
    main()
