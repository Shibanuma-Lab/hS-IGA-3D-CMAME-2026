#!/usr/bin/env python3
"""
Crack-velocity sweep driver for dynamic circular-crack simulations.

This sweep keeps the parameter-sweep baseline mesh fixed:

    rGL = 8, aL = 20, lL = 12, HL = 15

and runs one dynamic case per crack velocity.  By default it plans
velocities 200..1500 m/s in increments of 100 m/s.  Each velocity must have
the matching FEM boundary-condition file:

    data/FEMdata/FEM_v_<velocity>_uva.mat

Post-processing reuses postprocess_dynamic.py, so it writes the same local
stress, DSIF, and normalized comparison CSV files as the parameter sweep.
"""

import argparse
import csv
import math
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from param_sweep_dynamic import (  # noqa: E402
    BASE_RGL,
    FIELDNAMES,
    ConstFileEditor,
    DynamicParamSweep,
    SweepCase,
    baseline_HL,
    baseline_aL,
    baseline_lL,
    calculate_actual_hG,
    calculate_uniform_theta,
    format_velocity,
)


DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "results" / "crack_velocity_sweep_dynamic"
DEFAULT_VELOCITY_START = 200
DEFAULT_VELOCITY_STOP = 1500
DEFAULT_VELOCITY_STEP = 100

BASELINE_SOURCE = (
    "param_sweep_dynamic baseline: "
    f"rGL={BASE_RGL}, aL={baseline_aL()}, "
    f"lL={baseline_lL()}, HL={baseline_HL()}"
)

VELOCITY_FIELDNAMES = FIELDNAMES + [
    "fem_mat",
    "fem_local_stress",
    "baseline_source",
]


def velocity_file_label(velocity):
    """
    FEMDataLoader casts velocity to int before building FEM_v_<v>_uva.mat.
    Keep this driver explicit and only accept integer velocity labels.
    """
    value = float(velocity)
    nearest = round(value)
    if not math.isclose(value, nearest, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError(
            f"Velocity {velocity!r} is not an integer. "
            "FEM reference files use integer names such as FEM_v_500_uva.mat."
        )
    return str(int(nearest))


class CrackVelocitySweep(DynamicParamSweep):
    """
    One-factor sweep over crack velocity using the dynamic baseline mesh.
    """

    def __init__(self, *args, allow_missing_fem_data=False, **kwargs):
        super().__init__(
            *args,
            only_baseline=True,
            selected_groups=None,
            **kwargs,
        )
        self.allow_missing_fem_data = allow_missing_fem_data

    def build_cases(self):
        cases = []
        for idx, velocity in enumerate(self.velocities, start=1):
            hG = calculate_actual_hG(self.hL, BASE_RGL, self.mesh_constants)
            theta_info = calculate_uniform_theta(
                self.hL,
                self.theta_reference_radius,
            )
            v_label = format_velocity(velocity)
            cases.append(
                SweepCase(
                    idx=idx,
                    group="velocity",
                    label=f"v={v_label}",
                    v=velocity,
                    rGL=BASE_RGL,
                    aL=baseline_aL(),
                    lL=baseline_lL(),
                    HL=baseline_HL(),
                    d_theta=theta_info["d_theta"],
                    nLtheta=theta_info["nLtheta"],
                    hG=hG,
                    hL_theta_max=theta_info["hL_theta_max"],
                    actual_ratio=theta_info["actual_ratio"],
                    theta_reference_radius=theta_info["reference_radius"],
                )
            )
        return cases

    def fem_mat_path(self, case):
        v_label = velocity_file_label(case.v)
        return SCRIPT_DIR / "data" / "FEMdata" / f"FEM_v_{v_label}_uva.mat"

    def fem_local_stress_path(self, case):
        v_label = velocity_file_label(case.v)
        return (
            SCRIPT_DIR
            / "data"
            / "FEM_local_stress"
            / f"FEM_v_{v_label}_last.xlsx"
        )

    def add_reference_columns(self, row, case):
        row["fem_mat"] = str(self.fem_mat_path(case))
        row["fem_local_stress"] = str(self.fem_local_stress_path(case))
        row["baseline_source"] = BASELINE_SOURCE
        return row

    def missing_references_before_run(self, case, folder):
        """
        Return (status, message) when references needed for the next action are
        missing.  Return (None, None) when the case can proceed.
        """
        fem_mat = self.fem_mat_path(case)
        fem_local_stress = self.fem_local_stress_path(case)

        result_exists = self.result_exists(folder, ignore_force=True)
        if self.postprocess_only:
            reruns_simulation = False
            needs_postprocess = result_exists and (
                self.force or not self.postprocess_exists(folder, step=self.max_step)
            )
        else:
            reruns_simulation = self.force or not result_exists
            needs_postprocess = self.postprocess and (
                reruns_simulation
                or not self.postprocess_exists(folder, step=self.max_step)
            )
        needs_dsif_reference = needs_postprocess and not self.postprocess_skip_dsif
        needs_fem_mat = reruns_simulation or needs_dsif_reference

        if (
            needs_fem_mat
            and not self.allow_missing_fem_data
            and not fem_mat.exists()
        ):
            return (
                "missing_fem_data",
                (
                    f"Missing FEM MAT file for v={case.velocity_label}: "
                    f"{fem_mat}. Add this file before running this velocity."
                ),
            )

        if needs_postprocess and not fem_local_stress.exists():
            return (
                "missing_local_stress_reference",
                (
                    f"Missing FEM local-stress reference for v={case.velocity_label}: "
                    f"{fem_local_stress}."
                ),
            )

        return None, None

    def run(self):
        self.output_dir.mkdir(parents=True, exist_ok=True)
        summary_path = self.output_dir / "crack_velocity_sweep_summary.csv"
        cases = self.build_cases()
        const_editor = None

        if not self.dry_run and not self.postprocess_only:
            self.cleanup_workdirs()
            const_editor = ConstFileEditor(SCRIPT_DIR)

        try:
            with summary_path.open("w", newline="") as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=VELOCITY_FIELDNAMES)
                writer.writeheader()

                print("=" * 80)
                print("Dynamic crack velocity sweep")
                print("=" * 80)
                print(
                    f"Velocities: {', '.join(format_velocity(v) for v in self.velocities)}"
                )
                print(
                    "Baseline: "
                    f"rGL={BASE_RGL}, aL={baseline_aL()}, "
                    f"lL={baseline_lL()}, HL={baseline_HL()}"
                )
                print(
                    f"Steps: {self.step_start}..{self.max_step} "
                    f"({self.step_end - self.step_start} solver calls)"
                )
                print(f"hL: {self.hL:.12g}")
                print("Theta rule: Section 5.2-style uniform division")
                print(f"Theta reference radius: {self.theta_reference_radius:.12g}")
                print(f"Output directory: {self.output_dir}")
                print(f"Summary CSV: {summary_path}")
                print(f"Dry run: {self.dry_run}")
                print(f"Postprocess: {self.postprocess}")
                print(f"Postprocess only: {self.postprocess_only}")
                if self.allow_missing_fem_data:
                    print("Missing FEM MAT files are allowed; use with care.")

                planned_folders = set()
                for case in cases:
                    row = self.handle_velocity_case(case, const_editor, planned_folders)
                    writer.writerow(row)
                    csv_file.flush()
                    planned_folders.add(row["folder"])
        finally:
            if const_editor is not None:
                const_editor.restore()

        print(f"\nSummary written to: {summary_path}")

    def handle_velocity_case(self, case, const_editor, planned_folders):
        folder = self.output_dir / f"v{case.velocity_label}" / case.folder_name
        status, message = self.missing_references_before_run(case, folder)
        if status is not None:
            print(f"[{case.idx:02d}] MISS velocity {case.label}: {message}")
            row = case.csv_row(folder, status, 0.0, message)
            return self.add_reference_columns(row, case)

        row = self.handle_case(case, const_editor, planned_folders)
        return self.add_reference_columns(row, case)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run a dynamic crack-velocity sweep with the parameter-sweep "
            "baseline mesh."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--velocities",
        type=float,
        nargs="+",
        default=None,
        help=(
            "Velocity values to sweep. If omitted, uses "
            "200..1500 m/s with step 100."
        ),
    )
    parser.add_argument(
        "--velocity-start",
        type=int,
        default=DEFAULT_VELOCITY_START,
        help="First default velocity when --velocities is omitted.",
    )
    parser.add_argument(
        "--velocity-stop",
        type=int,
        default=DEFAULT_VELOCITY_STOP,
        help="Last default velocity when --velocities is omitted.",
    )
    parser.add_argument(
        "--velocity-step",
        type=int,
        default=DEFAULT_VELOCITY_STEP,
        help="Default velocity increment when --velocities is omitted.",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=None,
        help=(
            "Legacy step-count mode. If set, max_step becomes "
            "step_start + steps - 1."
        ),
    )
    parser.add_argument(
        "--max-step",
        type=int,
        default=100,
        help="Last dynamic step to run/check; main.py receives step_end=max_step+1.",
    )
    parser.add_argument(
        "--step-start",
        type=int,
        default=0,
        help="First step passed to main.py.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for copied result folders and the summary CSV.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write/check the sweep plan without running simulations or editing const files.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun even when a matching result folder already exists.",
    )
    parser.add_argument(
        "--postprocess",
        action="store_true",
        dest="postprocess",
        help="Run local-stress and DSIF post-processing after each completed case.",
    )
    parser.add_argument(
        "--no-postprocess",
        action="store_false",
        dest="postprocess",
        help="Do not run post-processing after each completed case.",
    )
    parser.add_argument(
        "--postprocess-skip-dsif",
        action="store_true",
        help="With postprocess modes, only calculate local-stress outputs.",
    )
    parser.add_argument(
        "--postprocess-only",
        action="store_true",
        help=(
            "Post-process existing velocity-sweep result folders only. "
            "Missing result folders are recorded in the summary and main.py is not run."
        ),
    )
    parser.add_argument(
        "--allow-missing-fem-data",
        action="store_true",
        help=(
            "Allow a simulation or DSIF postprocess to proceed when the matching "
            "FEM MAT file is missing. This preserves legacy zero-BC behavior and "
            "is normally not recommended."
        ),
    )
    parser.set_defaults(postprocess=True)
    return parser.parse_args()


def resolve_velocities(args):
    if args.velocities is not None:
        velocities = args.velocities
    else:
        if args.velocity_step <= 0:
            raise ValueError("--velocity-step must be positive")
        if args.velocity_stop < args.velocity_start:
            raise ValueError("--velocity-stop must be >= --velocity-start")
        velocities = list(
            range(args.velocity_start, args.velocity_stop + 1, args.velocity_step)
        )

    for velocity in velocities:
        velocity_file_label(velocity)
    return velocities


def main():
    args = parse_args()
    if args.steps is not None and args.steps <= 0:
        raise ValueError("--steps must be positive")
    if args.postprocess_only and not args.postprocess:
        raise ValueError("--postprocess-only cannot be combined with --no-postprocess")

    max_step = args.max_step
    if args.steps is not None:
        max_step = args.step_start + args.steps - 1
    if max_step < args.step_start:
        raise ValueError("--max-step must be >= --step-start")

    sweep = CrackVelocitySweep(
        velocities=resolve_velocities(args),
        output_dir=args.output_dir,
        step_start=args.step_start,
        max_step=max_step,
        dry_run=args.dry_run,
        force=args.force,
        postprocess=args.postprocess,
        postprocess_skip_dsif=args.postprocess_skip_dsif,
        postprocess_only=args.postprocess_only,
        allow_missing_fem_data=args.allow_missing_fem_data,
    )
    sweep.run()


if __name__ == "__main__":
    main()
