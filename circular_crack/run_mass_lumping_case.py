#!/usr/bin/env python3
"""
Run one existing dynamic circular-crack case again with an SFEM mass-lumping mode.

The script reads the source case's run_config.json, applies the same dynamic
parameters, runs main.py with the requested SFEM_MASS_LUMPING mode in the
environment, copies the generated step results to a separate output case
folder, and then runs the same local-stress/DSIF post-processing used by the
sweep drivers.
"""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SOURCE_CASE = (
    SCRIPT_DIR
    / "results"
    / "crack_velocity_sweep_dynamic"
    / "v500"
    / "v500_rGL8_aL20_lL12_HL15_dtheta1.139241"
)
DEFAULT_PHYSICAL_OUTPUT_ROOT = SCRIPT_DIR / "results" / "crack_velocity_sweep_dynamic_lumping"
DEFAULT_PREDICTOR_OUTPUT_ROOT = (
    SCRIPT_DIR / "results" / "crack_velocity_sweep_dynamic_lumped_predictor"
)

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from param_sweep_dynamic import (  # noqa: E402
    ConstFileEditor,
    SweepCase,
    add_j_integral_args,
    j_integral_overrides_from_args,
)


SUMMARY_FIELDS = [
    "source_case",
    "output_case",
    "status",
    "seconds",
    "step_start",
    "step_end",
    "mass_lumping",
    "mass_lumping_mode",
    "postprocess",
    "skip_dsif",
    "run_log",
    "message",
]


def load_json(path):
    with Path(path).open() as f:
        return json.load(f)


def default_output_case(source_case, mass_lumping_mode):
    source_case = Path(source_case).resolve()
    source_root = SCRIPT_DIR / "results" / "crack_velocity_sweep_dynamic"
    try:
        relative = source_case.relative_to(source_root)
    except ValueError:
        return source_case.parent / f"{source_case.name}_{mass_lumping_mode}"

    if mass_lumping_mode == "predictor":
        return DEFAULT_PREDICTOR_OUTPUT_ROOT / relative
    return DEFAULT_PHYSICAL_OUTPUT_ROOT / relative


def make_case_from_config(config):
    return SweepCase(
        idx=int(config.get("idx", 1)),
        group=str(config.get("group", "mass_lumping")),
        label=str(config.get("label", "mass_lumping")),
        v=float(config["velocity"]),
        rGL=int(config["rGL"]),
        aL=int(config["aL"]),
        lL=int(config["lL"]),
        HL=int(config["HL"]),
        d_theta=float(config["d_theta"]),
        nLtheta=int(config.get("nLtheta", round(90.0 / float(config["d_theta"])))),
        hG=float(config.get("hG", 0.0)),
        hL_theta_max=float(config.get("hL_theta_max", 0.0)),
        actual_ratio=float(config.get("actual_ratio", 0.0)),
        theta_reference_radius=float(config.get("theta_reference_radius", 0.0)),
    )


def step_name(step):
    return f"step{int(step):05d}"


def is_step_workdir(path):
    return path.is_dir() and re.fullmatch(r"step\d{5}", path.name) is not None


def clear_generated_workdirs():
    results_dir = SCRIPT_DIR / "results"
    if results_dir.exists():
        for step_dir in results_dir.iterdir():
            if is_step_workdir(step_dir):
                shutil.rmtree(step_dir)

    inputfiles_dir = SCRIPT_DIR / "inputfiles"
    if inputfiles_dir.exists():
        shutil.rmtree(inputfiles_dir)


def step_result_complete(step_dir):
    log_dir = step_dir / "log"
    visual_dir = step_dir / "visual"
    if not log_dir.is_dir() or not visual_dir.is_dir():
        return False

    required = [
        log_dir / "u.g.dat",
        log_dir / "u.l.dat",
        log_dir / "u_gl.l.dat",
    ]
    if any((not path.exists()) or path.stat().st_size == 0 for path in required):
        return False

    return any(path.suffix == ".vtu" for path in visual_dir.iterdir())


def missing_generated_results(step_start, step_end):
    missing = []
    for step in range(step_start, step_end):
        result_dir = SCRIPT_DIR / "results" / step_name(step)
        if not step_result_complete(result_dir):
            missing.append(step_name(step))
    return missing


def copy_generated_outputs(output_case, step_start, step_end):
    output_case.mkdir(parents=True, exist_ok=True)

    for step in range(step_start, step_end):
        name = step_name(step)
        source = SCRIPT_DIR / "results" / name
        destination = output_case / name
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(source, destination)

    source_input = SCRIPT_DIR / "inputfiles"
    if source_input.exists():
        destination_input = output_case / "inputfiles"
        if destination_input.exists():
            shutil.rmtree(destination_input)
        shutil.copytree(source_input, destination_input)


def mass_lumping_env_value(mode):
    if mode == "predictor":
        return "PREDICTOR"
    return "1"


def write_lumping_run_config(
    output_case,
    source_case,
    source_config,
    command,
    step_start,
    step_end,
    mass_lumping_mode,
):
    env_value = mass_lumping_env_value(mass_lumping_mode)
    config = dict(source_config)
    config.update(
        {
            "source_case": str(Path(source_case).resolve()),
            "mass_lumping": mass_lumping_mode == "physical",
            "mass_lumping_mode": mass_lumping_mode,
            "mass_lumped_predictor": mass_lumping_mode == "predictor",
            "environment": {"SFEM_MASS_LUMPING": env_value},
            "step_start": int(step_start),
            "step_end": int(step_end),
            "steps_run": int(step_end) - int(step_start),
            "command": command,
            "timestamp": datetime.now().isoformat(),
        }
    )
    path = output_case / "run_config.json"
    path.write_text(json.dumps(config, indent=2))
    return path


def apply_hl_from_config(const_editor, source_config):
    hL = source_config.get("hL")
    if hL is None:
        return

    local_path = const_editor.paths["local"]
    local_text = local_path.read_text()
    local_text = ConstFileEditor._replace_assignment(
        local_text,
        "hL",
        f"{float(hL):.12g}  # Length of local element [m]",
    )
    local_path.write_text(local_text)


def write_summary(path, row):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerow(row)


def summary_row(args, status, seconds, message, run_log):
    return {
        "source_case": str(args.source_case.resolve()),
        "output_case": str(args.output_case.resolve()),
        "status": status,
        "seconds": f"{seconds:.3f}",
        "step_start": args.step_start,
        "step_end": args.step_end,
        "mass_lumping": int(args.mass_lumping_mode == "physical"),
        "mass_lumping_mode": args.mass_lumping_mode,
        "postprocess": int(args.postprocess),
        "skip_dsif": int(args.skip_dsif),
        "run_log": str(run_log) if run_log else "",
        "message": message,
    }


def run_main_with_lumping(args, case, source_config):
    env_value = mass_lumping_env_value(args.mass_lumping_mode)
    command = [
        sys.executable,
        "main.py",
        "--step_start",
        str(args.step_start),
        "--step_end",
        str(args.step_end),
        "--delete",
    ]
    display_command = [f"SFEM_MASS_LUMPING={env_value}", *command]
    env = os.environ.copy()
    env["SFEM_MASS_LUMPING"] = env_value

    run_log = args.output_case / "run_lumping.log"
    const_editor = ConstFileEditor(SCRIPT_DIR)
    start = time.time()
    try:
        clear_generated_workdirs()
        const_editor.apply_case(case)
        apply_hl_from_config(const_editor, source_config)
        args.output_case.mkdir(parents=True, exist_ok=True)

        with run_log.open("w") as log_file:
            log_file.write(f"Command: {' '.join(display_command)}\n")
            log_file.write(f"Started: {datetime.now().isoformat()}\n")
            log_file.write(f"Source case: {args.source_case.resolve()}\n\n")
            completed = subprocess.run(
                command,
                cwd=SCRIPT_DIR,
                env=env,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                check=False,
            )

        seconds = time.time() - start
        if completed.returncode != 0:
            return (
                "failed",
                seconds,
                f"main.py exited with code {completed.returncode}; see {run_log}",
                run_log,
            )

        missing = missing_generated_results(args.step_start, args.step_end)
        if missing:
            message = f"Missing expected generated results: {', '.join(missing[:5])}"
            if len(missing) > 5:
                message += f" ... (+{len(missing) - 5} more)"
            return "failed", seconds, message, run_log

        copy_generated_outputs(args.output_case, args.step_start, args.step_end)
        write_lumping_run_config(
            args.output_case,
            args.source_case,
            source_config,
            display_command,
            args.step_start,
            args.step_end,
            args.mass_lumping_mode,
        )
        clear_generated_workdirs()
        return "done", seconds, "Mass-lumping simulation completed", run_log
    finally:
        const_editor.restore()


def run_postprocess(args):
    from postprocess_dynamic import postprocess_case

    outputs = postprocess_case(
        args.output_case,
        step=args.step_end - 1,
        output_dir=args.output_case / "postprocess",
        skip_dsif=args.skip_dsif,
        j_params=j_integral_overrides_from_args(args),
    )
    summary = Path(outputs["summary"])
    try:
        return summary.relative_to(args.output_case)
    except ValueError:
        return summary


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the v500 dynamic case with an SFEM mass-lumping mode.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--source-case",
        type=Path,
        default=DEFAULT_SOURCE_CASE,
        help="Existing non-lumped case folder whose run_config.json defines the case.",
    )
    parser.add_argument(
        "--output-case",
        type=Path,
        default=None,
        help=(
            "Destination for lumped results. Defaults to the same relative path "
            "under a mode-specific crack-velocity result directory."
        ),
    )
    parser.add_argument(
        "--mass-lumping-mode",
        choices=("physical", "predictor"),
        default="physical",
        help=(
            "physical uses SFEM_MASS_LUMPING=1 and changes the dynamic mass matrix; "
            "predictor uses SFEM_MASS_LUMPING=PREDICTOR, which keeps the final "
            "consistent-mass solve and only uses the lumped matrix as an initial guess."
        ),
    )
    parser.add_argument(
        "--step-start",
        type=int,
        default=None,
        help="First step to run; defaults to source run_config.json.",
    )
    parser.add_argument(
        "--step-end",
        type=int,
        default=None,
        help="Exclusive final step for main.py; defaults to source run_config.json.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Remove an existing output case before running.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the resolved paths and parameters without running main.py.",
    )
    parser.add_argument(
        "--postprocess",
        action="store_true",
        dest="postprocess",
        help="Run local-stress and DSIF post-processing after the simulation.",
    )
    parser.add_argument(
        "--no-postprocess",
        action="store_false",
        dest="postprocess",
        help="Skip post-processing.",
    )
    parser.add_argument(
        "--postprocess-only",
        action="store_true",
        help="Run post-processing on an existing output case without rerunning main.py.",
    )
    parser.add_argument(
        "--skip-dsif",
        action="store_true",
        help="Only calculate local-stress outputs during post-processing.",
    )
    add_j_integral_args(parser)
    parser.set_defaults(postprocess=True)
    return parser.parse_args()


def resolve_args(args):
    args.source_case = args.source_case.resolve()
    if not args.source_case.is_dir():
        raise FileNotFoundError(f"Source case folder not found: {args.source_case}")

    config_path = args.source_case / "run_config.json"
    if not config_path.exists():
        raise FileNotFoundError(f"Missing source run_config.json: {config_path}")

    source_config = load_json(config_path)
    if args.output_case is None:
        args.output_case = default_output_case(args.source_case, args.mass_lumping_mode)
    args.output_case = args.output_case.resolve()

    if args.step_start is None:
        args.step_start = int(source_config.get("step_start", 0))
    if args.step_end is None:
        args.step_end = int(source_config.get("step_end", 101))
    if args.step_start >= args.step_end:
        raise ValueError("--step-start must be smaller than --step-end")

    return args, source_config


def main():
    args, source_config = resolve_args(parse_args())
    case = make_case_from_config(source_config)
    summary_path = args.output_case / "mass_lumping_summary.csv"

    print("=" * 80)
    print("Dynamic mass-lumping single-case run")
    print("=" * 80)
    print(f"Source case: {args.source_case}")
    print(f"Output case: {args.output_case}")
    print(f"Steps: {args.step_start}..{args.step_end - 1}")
    print(
        "Case: "
        f"v={case.velocity_label}, rGL={case.rGL}, aL={case.aL}, "
        f"lL={case.lL}, HL={case.HL}, d_theta={case.d_theta:.12g}"
    )
    print(f"Mass-lumping mode: {args.mass_lumping_mode}")
    print(f"Environment: SFEM_MASS_LUMPING={mass_lumping_env_value(args.mass_lumping_mode)}")
    print(f"Postprocess: {args.postprocess}")
    print(f"Skip DSIF: {args.skip_dsif}")

    if args.dry_run:
        print("Dry run only; no files were changed.")
        return

    if args.postprocess_only:
        if not args.output_case.exists():
            raise FileNotFoundError(f"Output case does not exist: {args.output_case}")
        start = time.time()
        post_summary = run_postprocess(args)
        seconds = time.time() - start
        row = summary_row(
            args,
            "postprocessed",
            seconds,
            f"Postprocess completed: {post_summary}",
            None,
        )
        write_summary(summary_path, row)
        print(f"Postprocess completed: {post_summary}")
        print(f"Summary written to: {summary_path}")
        return

    if args.output_case.exists():
        if not args.force:
            raise FileExistsError(
                f"Output case already exists: {args.output_case}. "
                "Use --force to replace it or --postprocess-only to reuse it."
            )
        shutil.rmtree(args.output_case)

    status, seconds, message, run_log = run_main_with_lumping(args, case, source_config)

    if status == "done" and args.postprocess:
        try:
            post_summary = run_postprocess(args)
            status = "done_postprocessed"
            message = f"{message}; postprocess={post_summary}"
        except Exception as exc:
            status = "postprocess_failed"
            message = f"{message}; postprocess failed: {exc}"

    row = summary_row(args, status, seconds, message, run_log)
    write_summary(summary_path, row)

    print(message)
    print(f"Run log: {run_log}")
    print(f"Summary written to: {summary_path}")
    if status not in {"done", "done_postprocessed"}:
        sys.exit(1)


if __name__ == "__main__":
    main()
