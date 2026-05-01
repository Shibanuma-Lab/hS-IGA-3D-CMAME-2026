#!/usr/bin/env python3
"""
Parameter sweep driver for dynamic circular-crack simulations.

The sweep follows the Section 6.2-style one-factor-at-a-time plan:
vary rGL, aL, lL, and HL around the baseline rGL=8 condition, run each
configuration for a fixed number of dynamic steps, and record an incremental
CSV summary for later inspection.
"""

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent

FIELDNAMES = [
    "idx",
    "group",
    "label",
    "v",
    "rGL",
    "aL",
    "lL",
    "HL",
    "d_theta",
    "nLtheta",
    "hG",
    "hL_theta_max",
    "actual_ratio",
    "folder",
    "status",
    "seconds",
    "message",
]

BASE_RGL = 8
TARGET_THETA_RATIO = 1.2
MIN_LL = 6
RGL_VALUES = [3, 4, 6, 8, 10]
A_FACTORS = [
    ("15/10", 15.0 / 10.0),
    ("20/10", 20.0 / 10.0),
    ("25/10", 25.0 / 10.0),
    ("30/10", 30.0 / 10.0),
    ("40/10", 40.0 / 10.0),
]
LL_FACTORS = [
    ("8/10", 8.0 / 10.0),
    ("10/10", 10.0 / 10.0),
    ("12/10", 12.0 / 10.0),
    ("16/10", 16.0 / 10.0),
    ("22/10", 22.0 / 10.0),
]
HL_FACTORS = [
    ("8/10", 8.0 / 10.0),
    ("10/10", 10.0 / 10.0),
    ("12/10", 12.0 / 10.0),
    ("14/10", 14.0 / 10.0),
    ("16/10", 16.0 / 10.0),
    ("18/10", 18.0 / 10.0),
    ("20/10", 20.0 / 10.0),
    ("22/10", 22.0 / 10.0),
    ("24/10", 24.0 / 10.0),
]


def ceil_scaled(factor, rGL):
    return int(math.ceil(factor * rGL))


def baseline_aL(rGL=BASE_RGL):
    return ceil_scaled(25.0 / 10.0, rGL)


def baseline_lL(rGL=BASE_RGL):
    return max(MIN_LL, ceil_scaled(12.0 / 10.0, rGL))


def sweep_lL(factor, rGL=BASE_RGL):
    return max(MIN_LL, ceil_scaled(factor, rGL))


def baseline_HL(rGL=BASE_RGL):
    return ceil_scaled(18.0 / 10.0, rGL)


def format_velocity(value):
    value_float = float(value)
    if value_float.is_integer():
        return str(int(value_float))
    return f"{value_float:g}".replace(".", "p")


def format_float_for_path(value, decimals=6):
    text = f"{value:.{decimals}f}".rstrip("0").rstrip(".")
    return text if text else "0"


def load_current_dynamic_constants():
    """
    Read hL and global-mesh constants from the current dynamic const files.
    """
    if str(SCRIPT_DIR) not in sys.path:
        sys.path.insert(0, str(SCRIPT_DIR))

    from const import const_local_mesh as clm
    from const import const_global_mesh as cgm
    from const import simulation_params as sp

    return {
        "hL": float(clm.hL),
        "widthG": float(sp.WidthG),
        "mu_G": float(cgm.mu_G),
        "p": int(cgm.p),
    }


def calculate_actual_hG(hL, rGL, mesh_constants):
    """
    Match simulation_params.py: nPtsX = ceil(Lx/(hL*rGL)) + p.
    The actual element size is therefore Lx/(nPtsX-p).
    """
    Lx = mesh_constants["widthG"] * mesh_constants["mu_G"]
    n_elem_x = int(math.ceil(Lx / (hL * rGL)))
    return Lx / n_elem_x


def calculate_theta_from_ratio(target_ratio, hG, hL, lL, theta_max_step):
    """
    Use the Section 5.3 rule with the dynamic outer radius:

        R = theta_max_step*hL + lL*hL
        hL_theta_max = 2*R*sin(d_theta/2)
        hL_theta_max / hG = target_ratio

    The local mesh generator expects d_theta to divide 90 degrees, so choose
    the nearest admissible value d_theta = 90/nLtheta.
    """
    if target_ratio <= 0:
        raise ValueError(f"target_ratio must be positive, got {target_ratio}")

    outer_radius = (theta_max_step + lL) * hL
    target_arg = target_ratio * hG / (2.0 * outer_radius)
    if target_arg > 1.0:
        raise ValueError(
            f"target_ratio={target_ratio} is too large: "
            f"target_ratio*hG/(2R)={target_arg:.6f} > 1."
        )

    target_theta_rad = 2.0 * math.asin(target_arg)
    target_theta_deg = math.degrees(target_theta_rad)
    if target_theta_deg <= 0.0:
        raise ValueError(
            f"target_ratio={target_ratio} produced non-positive target theta."
        )

    ideal_nLtheta = 90.0 / target_theta_deg
    center = max(1, int(round(ideal_nLtheta)))
    candidate_min = max(1, center - 4)
    candidate_max = center + 4

    best = None
    for nLtheta in range(candidate_min, candidate_max + 1):
        d_theta = 90.0 / nLtheta
        hL_theta_max = 2.0 * outer_radius * math.sin(math.radians(d_theta) / 2.0)
        actual_ratio = hL_theta_max / hG
        ratio_error = actual_ratio - target_ratio
        theta_error = d_theta - target_theta_deg
        score = (abs(ratio_error), abs(theta_error), nLtheta)

        if best is None or score < best["score"]:
            best = {
                "score": score,
                "target_d_theta": target_theta_deg,
                "d_theta": d_theta,
                "nLtheta": nLtheta,
                "outer_radius": outer_radius,
                "hL_theta_max": hL_theta_max,
                "actual_ratio": actual_ratio,
                "ratio_error": ratio_error,
                "relative_ratio_error": ratio_error / target_ratio,
            }

    best.pop("score")
    return best


@dataclass(frozen=True)
class SweepCase:
    idx: int
    group: str
    label: str
    v: float
    rGL: int
    aL: int
    lL: int
    HL: int
    d_theta: float
    nLtheta: int
    hG: float
    hL_theta_max: float
    actual_ratio: float
    target_ratio: float
    theta_max_step: int

    @property
    def velocity_label(self):
        return format_velocity(self.v)

    @property
    def folder_name(self):
        dtheta = format_float_for_path(self.d_theta)
        return (
            f"v{self.velocity_label}_rGL{self.rGL}_"
            f"aL{self.aL}_lL{self.lL}_HL{self.HL}_dtheta{dtheta}"
        )

    def csv_row(self, folder, status, seconds, message):
        return {
            "idx": self.idx,
            "group": self.group,
            "label": self.label,
            "v": self.velocity_label,
            "rGL": self.rGL,
            "aL": self.aL,
            "lL": self.lL,
            "HL": self.HL,
            "d_theta": f"{self.d_theta:.12g}",
            "nLtheta": self.nLtheta,
            "hG": f"{self.hG:.12g}",
            "hL_theta_max": f"{self.hL_theta_max:.12g}",
            "actual_ratio": f"{self.actual_ratio:.12g}",
            "folder": str(folder),
            "status": status,
            "seconds": f"{seconds:.3f}",
            "message": message,
        }

    def run_config(self, step_start, step_end, command, hL=None):
        return {
            "idx": self.idx,
            "group": self.group,
            "label": self.label,
            "velocity": self.v,
            "rGL": self.rGL,
            "aL": self.aL,
            "lL": self.lL,
            "HL": self.HL,
            "hL": hL,
            "target_ratio": self.target_ratio,
            "theta_max_step": self.theta_max_step,
            "d_theta": self.d_theta,
            "nLtheta": self.nLtheta,
            "hG": self.hG,
            "hL_theta_max": self.hL_theta_max,
            "actual_ratio": self.actual_ratio,
            "step_start": step_start,
            "step_end": step_end,
            "steps_run": step_end - step_start,
            "command": command,
            "timestamp": datetime.now().isoformat(),
        }


class ConstFileEditor:
    """
    Temporarily update dynamic-mode constants for a single subprocess run.
    """

    def __init__(self, script_dir):
        self.paths = {
            "local": Path(script_dir) / "const" / "const_local_mesh.py",
            "global": Path(script_dir) / "const" / "const_global_mesh.py",
            "simulation": Path(script_dir) / "const" / "simulation_params.py",
        }
        self.original_text = {
            name: path.read_text()
            for name, path in self.paths.items()
        }

    def restore(self):
        for name, path in self.paths.items():
            path.write_text(self.original_text[name])

    def apply_case(self, case):
        local_text = self.original_text["local"]
        local_text = self._replace_assignment(local_text, "aL", str(case.aL))
        local_text = self._replace_assignment(local_text, "lL", str(case.lL))
        local_text = self._replace_assignment(local_text, "HL", str(case.HL))
        local_text = self._replace_assignment(
            local_text,
            "d_theta",
            f"{case.d_theta:.12f}  # Angular resolution [degrees]",
        )
        self.paths["local"].write_text(local_text)

        global_text = self.original_text["global"]
        global_text = self._replace_assignment(
            global_text,
            "rGL",
            f"{case.rGL}  # Ratio of global to local element size (hG/hL)",
        )
        self.paths["global"].write_text(global_text)

        simulation_text = self.original_text["simulation"]
        simulation_text = self._replace_assignment(
            simulation_text,
            "V",
            f"{float(case.v)}  # Velocity [m/s]",
        )
        self.paths["simulation"].write_text(simulation_text)

    @staticmethod
    def _replace_assignment(text, name, replacement):
        pattern = re.compile(rf"^{re.escape(name)}\s*=.*$", re.MULTILINE)
        new_text, count = pattern.subn(f"{name} = {replacement}", text, count=1)
        if count != 1:
            raise ValueError(f"Could not find assignment for {name}")
        return new_text


class DynamicParamSweep:
    def __init__(
        self,
        velocities,
        output_dir,
        step_start=0,
        max_step=100,
        theta_max_step=None,
        target_theta_ratio=TARGET_THETA_RATIO,
        mesh_constants=None,
        dry_run=False,
        force=False,
        postprocess=True,
        postprocess_skip_dsif=False,
        only_baseline=False,
    ):
        self.velocities = [float(v) for v in velocities]
        self.output_dir = Path(output_dir).resolve()
        self.step_start = step_start
        self.max_step = max_step
        self.step_end = max_step + 1
        self.theta_max_step = theta_max_step if theta_max_step is not None else max_step
        self.target_theta_ratio = target_theta_ratio
        self.mesh_constants = mesh_constants or load_current_dynamic_constants()
        self.hL = self.mesh_constants["hL"]
        self.dry_run = dry_run
        self.force = force
        self.postprocess = postprocess
        self.postprocess_skip_dsif = postprocess_skip_dsif
        self.only_baseline = only_baseline

    def build_cases_for_velocity(self, velocity):
        cases = []

        if self.only_baseline:
            cases.append(
                (
                    "baseline",
                    "baseline",
                    BASE_RGL,
                    baseline_aL(),
                    baseline_lL(),
                    baseline_HL(),
                )
            )
        else:

            for rGL in RGL_VALUES:
                cases.append(
                    (
                        "rGL",
                        f"rGL={rGL}",
                        rGL,
                        baseline_aL(rGL),
                        baseline_lL(rGL),
                        baseline_HL(rGL),
                    )
                )

            for label, factor in A_FACTORS:
                cases.append(
                    (
                        "aL",
                        f"aL=ceil({label}*rGL)",
                        BASE_RGL,
                        ceil_scaled(factor, BASE_RGL),
                        baseline_lL(),
                        baseline_HL(),
                    )
                )

            for label, factor in LL_FACTORS:
                cases.append(
                    (
                        "lL",
                        f"lL=ceil({label}*rGL)",
                        BASE_RGL,
                        baseline_aL(),
                        sweep_lL(factor, BASE_RGL),
                        baseline_HL(),
                    )
                )

            for label, factor in HL_FACTORS:
                cases.append(
                    (
                        "HL",
                        f"HL=ceil({label}*rGL)",
                        BASE_RGL,
                        baseline_aL(),
                        baseline_lL(),
                        ceil_scaled(factor, BASE_RGL),
                    )
                )

        sweep_cases = []
        for i, (group, label, rGL, aL, lL, HL) in enumerate(cases, start=1):
            hG = calculate_actual_hG(self.hL, rGL, self.mesh_constants)
            theta_info = calculate_theta_from_ratio(
                self.target_theta_ratio,
                hG,
                self.hL,
                lL,
                self.theta_max_step,
            )
            sweep_cases.append(
                SweepCase(
                    idx=i,
                    group=group,
                    label=label,
                    v=velocity,
                    rGL=rGL,
                    aL=aL,
                    lL=lL,
                    HL=HL,
                    d_theta=theta_info["d_theta"],
                    nLtheta=theta_info["nLtheta"],
                    hG=hG,
                    hL_theta_max=theta_info["hL_theta_max"],
                    actual_ratio=theta_info["actual_ratio"],
                    target_ratio=self.target_theta_ratio,
                    theta_max_step=self.theta_max_step,
                )
            )

        return sweep_cases

    def run(self):
        self.output_dir.mkdir(parents=True, exist_ok=True)
        const_editor = None

        if not self.dry_run:
            self.cleanup_workdirs()
            const_editor = ConstFileEditor(SCRIPT_DIR)

        try:
            for velocity in self.velocities:
                self.run_velocity(velocity, const_editor)
        finally:
            if const_editor is not None:
                const_editor.restore()

    def run_velocity(self, velocity, const_editor):
        v_label = format_velocity(velocity)
        summary_path = self.output_dir / f"param_sweep_v{v_label}_summary.csv"
        cases = self.build_cases_for_velocity(velocity)

        with summary_path.open("w", newline="") as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=FIELDNAMES)
            writer.writeheader()

            print("=" * 80)
            print(f"Dynamic parameter sweep: v={v_label} m/s")
            print("=" * 80)
            print(f"Steps: {self.step_start}..{self.max_step} ({self.step_end - self.step_start} solver calls)")
            print(f"hL: {self.hL:.12g}")
            print(f"Target hL_theta_max/hG: {self.target_theta_ratio}")
            print(f"Theta max step: {self.theta_max_step}")
            print(f"Output directory: {self.output_dir}")
            print(f"Summary CSV: {summary_path}")
            print(f"Dry run: {self.dry_run}")
            print(f"Postprocess: {self.postprocess}")

            planned_folders = set()
            for case in cases:
                row = self.handle_case(case, const_editor, planned_folders)
                writer.writerow(row)
                csv_file.flush()
                planned_folders.add(row["folder"])

        print(f"\nSummary written to: {summary_path}")

    def handle_case(self, case, const_editor, planned_folders=None):
        folder = self.output_dir / f"v{case.velocity_label}" / case.folder_name
        folder_key = str(folder)

        if self.result_exists(folder):
            message = "Existing result folder contains the expected final-step outputs"
            print(f"[{case.idx:02d}] SKIP {case.group:>3} {case.label}: {folder}")
            if self.postprocess and not self.dry_run:
                try:
                    from postprocess_dynamic import postprocess_case

                    outputs = postprocess_case(
                        folder,
                        step=self.max_step,
                        skip_dsif=self.postprocess_skip_dsif,
                    )
                    message = f"{message}; postprocess={outputs['summary']}"
                    return case.csv_row(folder, "skipped_postprocessed", 0.0, message)
                except Exception as exc:
                    message = f"{message}; postprocess failed: {exc}"
                    return case.csv_row(folder, "skipped_postprocess_failed", 0.0, message)
            if self.postprocess and self.dry_run:
                message = f"{message}; dry-run does not execute postprocess"
            return case.csv_row(folder, "skipped", 0.0, message)

        if self.dry_run and planned_folders is not None and folder_key in planned_folders:
            message = "Duplicate of an earlier planned case; formal run will skip after the first result exists"
            print(
                f"[{case.idx:02d}] DUP  {case.group:>3} {case.label}: "
                f"same folder {folder}"
            )
            return case.csv_row(folder, "duplicate_plan", 0.0, message)

        if self.dry_run:
            message = "Dry run only; simulation was not executed"
            print(
                f"[{case.idx:02d}] DRY  {case.group:>3} {case.label}: "
                f"rGL={case.rGL}, aL={case.aL}, lL={case.lL}, HL={case.HL}, "
                f"d_theta={case.d_theta:.6f}, ratio={case.actual_ratio:.6f}"
            )
            return case.csv_row(folder, "would_run", 0.0, message)

        if const_editor is None:
            raise RuntimeError("Internal error: const editor missing for non-dry run")

        folder.mkdir(parents=True, exist_ok=True)
        command = [
            "python3",
            "main.py",
            "--step_start",
            str(self.step_start),
            "--step_end",
            str(self.step_end),
            "--delete",
        ]

        if self.force:
            self.clear_case_outputs(folder)

        print(
            f"[{case.idx:02d}] RUN  {case.group:>3} {case.label}: "
            f"v={case.velocity_label}, rGL={case.rGL}, "
            f"aL={case.aL}, lL={case.lL}, HL={case.HL}, "
            f"d_theta={case.d_theta:.6f}"
        )
        print(f"     Result folder: {folder}")

        start_time = time.time()
        log_path = folder / "param_sweep_run.log"

        try:
            const_editor.apply_case(case)
            self.clear_source_outputs()

            with log_path.open("w") as log_file:
                log_file.write(f"Command: {' '.join(command)}\n")
                log_file.write(f"Started: {datetime.now().isoformat()}\n\n")
                completed = subprocess.run(
                    command,
                    cwd=SCRIPT_DIR,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=False,
                )

            seconds = time.time() - start_time
            if completed.returncode != 0:
                message = f"main.py exited with code {completed.returncode}; see {log_path}"
                return case.csv_row(folder, "failed", seconds, message)

            missing = self.missing_source_results()
            if missing:
                message = f"Missing expected source results after run: {', '.join(missing[:5])}"
                if len(missing) > 5:
                    message += f" ... (+{len(missing) - 5} more)"
                return case.csv_row(folder, "failed", seconds, message)

            self.copy_outputs(folder)
            self.cleanup_workdirs()
            config_path = folder / "run_config.json"
            config_path.write_text(
                json.dumps(
                    case.run_config(self.step_start, self.step_end, command, hL=self.hL),
                    indent=2,
                )
            )

            if self.postprocess:
                try:
                    from postprocess_dynamic import postprocess_case

                    outputs = postprocess_case(
                        folder,
                        step=self.max_step,
                        skip_dsif=self.postprocess_skip_dsif,
                    )
                    message = f"Completed; postprocess={outputs['summary']}"
                    return case.csv_row(folder, "done_postprocessed", seconds, message)
                except Exception as exc:
                    message = f"Simulation completed, but postprocess failed: {exc}"
                    return case.csv_row(folder, "postprocess_failed", seconds, message)

            message = f"Completed and copied outputs; log={log_path.name}"
            return case.csv_row(folder, "done", seconds, message)

        except Exception as exc:
            seconds = time.time() - start_time
            return case.csv_row(folder, "failed", seconds, str(exc))

    def clear_case_outputs(self, folder):
        for step in range(self.step_start, self.step_end):
            step_dir = folder / f"step{step:05d}"
            if step_dir.exists():
                shutil.rmtree(step_dir)

        input_dir = folder / "inputfiles"
        if input_dir.exists():
            shutil.rmtree(input_dir)

    def clear_source_outputs(self):
        for step in range(self.step_start, self.step_end):
            source = SCRIPT_DIR / "results" / f"step{step:05d}"
            if source.exists():
                shutil.rmtree(source)

    def cleanup_workdirs(self):
        results_dir = SCRIPT_DIR / "results"
        if results_dir.exists():
            for step_dir in results_dir.glob("step*"):
                if step_dir.is_dir():
                    shutil.rmtree(step_dir)

        inputfiles_dir = SCRIPT_DIR / "inputfiles"
        if inputfiles_dir.exists():
            shutil.rmtree(inputfiles_dir)

    def copy_outputs(self, folder):
        folder.mkdir(parents=True, exist_ok=True)

        for step in range(self.step_start, self.step_end):
            step_name = f"step{step:05d}"
            source = SCRIPT_DIR / "results" / step_name
            destination = folder / step_name
            if destination.exists():
                shutil.rmtree(destination)
            shutil.copytree(source, destination)

        source_input = SCRIPT_DIR / "inputfiles"
        if source_input.exists():
            destination_input = folder / "inputfiles"
            if destination_input.exists():
                shutil.rmtree(destination_input)
            shutil.copytree(source_input, destination_input)

    def result_exists(self, folder):
        if self.force or not folder.exists():
            return False

        final_step = self.step_end - 1
        final_dir = folder / f"step{final_step:05d}"
        log_dir = final_dir / "log"
        visual_dir = final_dir / "visual"

        if not log_dir.is_dir() or not visual_dir.is_dir():
            return False

        essential_files = [
            log_dir / "u.g.dat",
            log_dir / "u.l.dat",
            log_dir / "u_gl.l.dat",
        ]
        if any((not path.exists()) or path.stat().st_size == 0 for path in essential_files):
            return False

        return any(path.suffix == ".vtu" for path in visual_dir.iterdir())

    def missing_source_results(self):
        missing = []
        for step in range(self.step_start, self.step_end):
            step_name = f"step{step:05d}"
            step_dir = SCRIPT_DIR / "results" / step_name
            log_dir = step_dir / "log"
            visual_dir = step_dir / "visual"

            if not log_dir.is_dir() or not visual_dir.is_dir():
                missing.append(step_name)
                continue

            essential_files = [
                log_dir / "u.g.dat",
                log_dir / "u.l.dat",
                log_dir / "u_gl.l.dat",
            ]
            if any((not path.exists()) or path.stat().st_size == 0 for path in essential_files):
                missing.append(step_name)
                continue

            if not any(path.suffix == ".vtu" for path in visual_dir.iterdir()):
                missing.append(step_name)

        return missing


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run dynamic circular-crack parameter sweeps for v=500/1000.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--velocities",
        type=float,
        nargs="+",
        default=[500.0, 1000.0],
        help="Velocity values to sweep",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=None,
        help=(
            "Legacy step-count mode. If set, max_step becomes "
            "step_start + steps - 1"
        ),
    )
    parser.add_argument(
        "--max-step",
        type=int,
        default=100,
        help="Last dynamic step to run/check; main.py receives step_end=max_step+1",
    )
    parser.add_argument(
        "--theta-max-step",
        type=int,
        default=None,
        help="Step used in R=(theta_max_step+lL)*hL for d_theta; defaults to max_step",
    )
    parser.add_argument(
        "--theta-ratio",
        type=float,
        default=TARGET_THETA_RATIO,
        help="Target hL_theta_max/hG ratio used to calculate d_theta",
    )
    parser.add_argument(
        "--step-start",
        type=int,
        default=0,
        help="First step passed to main.py",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR / "results" / "param_sweep_dynamic",
        help="Directory for copied result folders and summary CSV files",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write/check the sweep plan without running simulations or editing const files",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun even when a matching result folder already exists",
    )
    parser.add_argument(
        "--postprocess",
        action="store_true",
        dest="postprocess",
        help="Run local-stress and DSIF post-processing after each completed case (default)",
    )
    parser.add_argument(
        "--no-postprocess",
        action="store_false",
        dest="postprocess",
        help="Do not run post-processing after each completed case",
    )
    parser.add_argument(
        "--postprocess-skip-dsif",
        action="store_true",
        help="With --postprocess, only calculate local-stress outputs",
    )
    parser.add_argument(
        "--only-baseline",
        action="store_true",
        help="Run only the baseline case rGL=8, aL=20, lL=10, HL=15",
    )
    parser.set_defaults(postprocess=True)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.steps is not None and args.steps <= 0:
        raise ValueError("--steps must be positive")
    if args.theta_ratio <= 0.0:
        raise ValueError("--theta-ratio must be positive")

    max_step = args.max_step
    if args.steps is not None:
        max_step = args.step_start + args.steps - 1
    if max_step < args.step_start:
        raise ValueError("--max-step must be >= --step-start")
    if args.theta_max_step is not None and args.theta_max_step < 0:
        raise ValueError("--theta-max-step must be non-negative")

    sweep = DynamicParamSweep(
        velocities=args.velocities,
        output_dir=args.output_dir,
        step_start=args.step_start,
        max_step=max_step,
        theta_max_step=args.theta_max_step,
        target_theta_ratio=args.theta_ratio,
        dry_run=args.dry_run,
        force=args.force,
        postprocess=args.postprocess,
        postprocess_skip_dsif=args.postprocess_skip_dsif,
        only_baseline=args.only_baseline,
    )
    sweep.run()


if __name__ == "__main__":
    main()
