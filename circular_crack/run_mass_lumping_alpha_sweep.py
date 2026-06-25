#!/usr/bin/env python3
"""
Run one dynamic circular-crack case for several blended mass-lumping alphas.

Each run uses

    SFEM_MASS_LUMPING_ALPHA=<alpha> python main.py ...

then copies the generated step results into an alpha-specific result folder
and runs the same dynamic local-stress/DSIF post-processing used by the sweep
drivers.
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
DEFAULT_OUTPUT_ROOT = (
    SCRIPT_DIR / "results" / "crack_velocity_sweep_dynamic_mass_lumping_alpha"
)
DEFAULT_ALPHAS = (0.01, 0.02, 0.05)

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from param_sweep_dynamic import (  # noqa: E402
    ConstFileEditor,
    SweepCase,
    add_j_integral_args,
    j_integral_overrides_from_args,
)


SUMMARY_FIELDS = [
    "idx",
    "alpha",
    "alpha_label",
    "source_case",
    "output_case",
    "status",
    "run_seconds",
    "postprocess_seconds",
    "step_start",
    "step_end",
    "postprocess",
    "skip_dsif",
    "run_log",
    "message",
]


def load_json(path):
    with Path(path).open() as f:
        return json.load(f)


def format_alpha_for_path(alpha):
    text = f"{float(alpha):.6g}".replace("-", "m").replace(".", "p")
    return f"alpha{text}"


def make_case_from_config(config):
    return SweepCase(
        idx=int(config.get("idx", 1)),
        group=str(config.get("group", "mass_lumping_alpha")),
        label=str(config.get("label", "mass_lumping_alpha")),
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


def source_relative_case(source_case):
    source_root = SCRIPT_DIR / "results" / "crack_velocity_sweep_dynamic"
    try:
        return Path(source_case).resolve().relative_to(source_root)
    except ValueError:
        return Path(source_case).resolve().name


def output_case_for_alpha(output_root, source_case, alpha):
    return Path(output_root).resolve() / format_alpha_for_path(alpha) / source_relative_case(source_case)


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


def result_exists(case_folder, final_step):
    return step_result_complete(Path(case_folder) / step_name(final_step))


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


def write_alpha_run_config(output_case, source_case, source_config, command, step_start, step_end, alpha):
    config = dict(source_config)
    config.update(
        {
            "source_case": str(Path(source_case).resolve()),
            "mass_lumping_alpha": float(alpha),
            "environment": {"SFEM_MASS_LUMPING_ALPHA": f"{float(alpha):.12g}"},
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


def run_postprocess(output_case, target_step, args):
    from postprocess_dynamic import postprocess_case

    start = time.time()
    outputs = postprocess_case(
        output_case,
        step=target_step,
        output_dir=output_case / "postprocess",
        skip_dsif=args.skip_dsif,
        j_params=j_integral_overrides_from_args(args),
    )
    seconds = time.time() - start
    summary = Path(outputs["summary"])
    try:
        summary = summary.relative_to(output_case)
    except ValueError:
        pass
    return seconds, summary


def run_alpha(idx, alpha, args, case, source_config):
    output_case = output_case_for_alpha(args.output_root, args.source_case, alpha)
    final_step = args.step_end - 1
    run_seconds = 0.0
    post_seconds = 0.0
    run_log = output_case / f"run_{format_alpha_for_path(alpha)}.log"

    if args.dry_run:
        return summary_row(
            idx,
            alpha,
            args,
            output_case,
            "would_run",
            run_seconds,
            post_seconds,
            run_log,
            "Dry run only; simulation was not executed",
        )

    if args.postprocess_only:
        if not result_exists(output_case, final_step):
            return summary_row(
                idx,
                alpha,
                args,
                output_case,
                "missing_result",
                run_seconds,
                post_seconds,
                run_log,
                "No existing final-step outputs; postprocess-only did not run simulation",
            )
        try:
            post_seconds, post_summary = run_postprocess(output_case, final_step, args)
            return summary_row(
                idx,
                alpha,
                args,
                output_case,
                "postprocessed",
                run_seconds,
                post_seconds,
                run_log,
                f"Postprocess completed: {post_summary}",
            )
        except Exception as exc:
            return summary_row(
                idx,
                alpha,
                args,
                output_case,
                "postprocess_failed",
                run_seconds,
                post_seconds,
                run_log,
                f"Postprocess failed: {exc}",
            )

    if output_case.exists():
        if args.force:
            shutil.rmtree(output_case)
        else:
            return summary_row(
                idx,
                alpha,
                args,
                output_case,
                "exists",
                run_seconds,
                post_seconds,
                run_log,
                "Output case already exists; use --force to replace it",
            )

    env_value = f"{float(alpha):.12g}"
    command = [
        sys.executable,
        "main.py",
        "--step_start",
        str(args.step_start),
        "--step_end",
        str(args.step_end),
        "--delete",
    ]
    display_command = [f"SFEM_MASS_LUMPING_ALPHA={env_value}", *command]
    env = os.environ.copy()
    env["SFEM_MASS_LUMPING_ALPHA"] = env_value
    env.pop("SFEM_MASS_LUMPING", None)

    const_editor = ConstFileEditor(SCRIPT_DIR)
    start = time.time()
    try:
        clear_generated_workdirs()
        const_editor.apply_case(case)
        apply_hl_from_config(const_editor, source_config)
        output_case.mkdir(parents=True, exist_ok=True)

        with run_log.open("w") as log_file:
            log_file.write(f"Command: {' '.join(display_command)}\n")
            log_file.write(f"Started: {datetime.now().isoformat()}\n")
            log_file.write(f"Source case: {args.source_case.resolve()}\n")
            log_file.write(f"Mass lumping alpha: {env_value}\n\n")
            completed = subprocess.run(
                command,
                cwd=SCRIPT_DIR,
                env=env,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                check=False,
            )

        run_seconds = time.time() - start
        if completed.returncode != 0:
            return summary_row(
                idx,
                alpha,
                args,
                output_case,
                "failed",
                run_seconds,
                post_seconds,
                run_log,
                f"main.py exited with code {completed.returncode}; see {run_log}",
            )

        missing = missing_generated_results(args.step_start, args.step_end)
        if missing:
            message = f"Missing expected generated results: {', '.join(missing[:5])}"
            if len(missing) > 5:
                message += f" ... (+{len(missing) - 5} more)"
            return summary_row(
                idx,
                alpha,
                args,
                output_case,
                "failed",
                run_seconds,
                post_seconds,
                run_log,
                message,
            )

        copy_generated_outputs(output_case, args.step_start, args.step_end)
        write_alpha_run_config(
            output_case,
            args.source_case,
            source_config,
            display_command,
            args.step_start,
            args.step_end,
            alpha,
        )
        clear_generated_workdirs()

        status = "done"
        message = "Simulation completed"
        if args.postprocess:
            try:
                post_seconds, post_summary = run_postprocess(output_case, final_step, args)
                status = "done_postprocessed"
                message = f"{message}; postprocess={post_summary}"
            except Exception as exc:
                status = "postprocess_failed"
                message = f"{message}; postprocess failed: {exc}"

        return summary_row(
            idx,
            alpha,
            args,
            output_case,
            status,
            run_seconds,
            post_seconds,
            run_log,
            message,
        )
    finally:
        const_editor.restore()


def summary_row(idx, alpha, args, output_case, status, run_seconds, post_seconds, run_log, message):
    return {
        "idx": idx,
        "alpha": f"{float(alpha):.12g}",
        "alpha_label": format_alpha_for_path(alpha),
        "source_case": str(args.source_case.resolve()),
        "output_case": str(Path(output_case).resolve()),
        "status": status,
        "run_seconds": f"{run_seconds:.3f}",
        "postprocess_seconds": f"{post_seconds:.3f}",
        "step_start": args.step_start,
        "step_end": args.step_end,
        "postprocess": int(args.postprocess),
        "skip_dsif": int(args.skip_dsif),
        "run_log": str(run_log) if run_log else "",
        "message": message,
    }


def write_summary_header(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()


def append_summary_row(path, row):
    with path.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS)
        writer.writerow(row)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run one dynamic case for several SFEM_MASS_LUMPING_ALPHA values.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--source-case",
        type=Path,
        default=DEFAULT_SOURCE_CASE,
        help="Existing reference case folder whose run_config.json defines the case.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Root directory for alpha-specific result folders and summary CSV.",
    )
    parser.add_argument(
        "--alphas",
        type=float,
        nargs="+",
        default=list(DEFAULT_ALPHAS),
        help="Mass-lumping alpha values to run sequentially.",
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
        help="Remove existing alpha result folders before running.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print and summarize the resolved run plan without running main.py.",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue to later alpha values after a failed run.",
    )
    parser.add_argument(
        "--postprocess",
        action="store_true",
        dest="postprocess",
        help="Run local-stress and DSIF post-processing after each simulation.",
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
        help="Post-process existing alpha result folders without rerunning main.py.",
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
    args.output_root = args.output_root.resolve()
    if not args.source_case.is_dir():
        raise FileNotFoundError(f"Source case folder not found: {args.source_case}")

    config_path = args.source_case / "run_config.json"
    if not config_path.exists():
        raise FileNotFoundError(f"Missing source run_config.json: {config_path}")

    source_config = load_json(config_path)
    if args.step_start is None:
        args.step_start = int(source_config.get("step_start", 0))
    if args.step_end is None:
        args.step_end = int(source_config.get("step_end", 101))
    if args.step_start >= args.step_end:
        raise ValueError("--step-start must be smaller than --step-end")

    normalized = []
    seen = set()
    for alpha in args.alphas:
        alpha = float(alpha)
        if alpha < 0.0 or alpha > 1.0:
            raise ValueError(f"Alpha must be in [0, 1], got {alpha}")
        label = format_alpha_for_path(alpha)
        if label not in seen:
            normalized.append(alpha)
            seen.add(label)
    if not normalized:
        raise ValueError("--alphas must include at least one alpha value")
    args.alphas = normalized

    if args.postprocess_only and not args.postprocess:
        raise ValueError("--postprocess-only cannot be combined with --no-postprocess")

    return args, source_config


def main():
    args, source_config = resolve_args(parse_args())
    case = make_case_from_config(source_config)
    summary_path = args.output_root / "mass_lumping_alpha_sweep_summary.csv"

    print("=" * 80)
    print("Dynamic mass-lumping alpha sweep")
    print("=" * 80)
    print(f"Source case: {args.source_case}")
    print(f"Output root: {args.output_root}")
    print(f"Alphas: {', '.join(f'{alpha:.12g}' for alpha in args.alphas)}")
    print(f"Steps: {args.step_start}..{args.step_end - 1}")
    print(
        "Case: "
        f"v={case.velocity_label}, rGL={case.rGL}, aL={case.aL}, "
        f"lL={case.lL}, HL={case.HL}, d_theta={case.d_theta:.12g}"
    )
    print(f"Postprocess: {args.postprocess}")
    print(f"Skip DSIF: {args.skip_dsif}")
    print(f"Dry run: {args.dry_run}")
    print(f"Summary CSV: {summary_path}")

    if not args.dry_run:
        write_summary_header(summary_path)

    failed = False
    for idx, alpha in enumerate(args.alphas, start=1):
        output_case = output_case_for_alpha(args.output_root, args.source_case, alpha)
        print("-" * 80)
        print(f"[{idx}/{len(args.alphas)}] alpha={alpha:.12g}")
        print(f"Output case: {output_case}")
        if not args.dry_run:
            print(f"Environment: SFEM_MASS_LUMPING_ALPHA={alpha:.12g}")

        row = run_alpha(idx, alpha, args, case, source_config)
        if not args.dry_run:
            append_summary_row(summary_path, row)
        print(f"Status: {row['status']}")
        print(f"Message: {row['message']}")

        if row["status"] in {"failed", "postprocess_failed", "missing_result", "exists"}:
            failed = True
            if not args.keep_going:
                print("Stopping after this alpha. Use --keep-going to continue after failures.")
                break

    if args.dry_run:
        print("\nDry run only; no files were changed.")
    else:
        print(f"\nSummary written to: {summary_path}")
    if failed and not args.dry_run:
        sys.exit(1)


if __name__ == "__main__":
    main()
