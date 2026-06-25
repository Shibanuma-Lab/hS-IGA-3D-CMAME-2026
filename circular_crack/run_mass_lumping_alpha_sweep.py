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
DEFAULT_VELOCITY = 500.0

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from param_sweep_dynamic import (  # noqa: E402
    BASE_RGL,
    ConstFileEditor,
    SweepCase,
    add_j_integral_args,
    baseline_HL,
    baseline_aL,
    baseline_lL,
    calculate_actual_hG,
    calculate_uniform_theta,
    format_velocity,
    j_integral_overrides_from_args,
    load_current_dynamic_constants,
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
    "timing_summary",
    "timing_by_step",
    "timing_steps_found",
    "timing_steps_missing",
    "timing_log_total_time_sec",
    "timing_requested_phase_time_sec",
    "timing_input_time_sec",
    "timing_matrix_generation_time_sec",
    "timing_solver_time_sec",
    "timing_output_time_sec",
    "timing_other_time_sec",
    "message",
]

TIME_FIELDS = [
    "input_time_sec",
    "nonzero_detection_time_sec",
    "K_GG_assembly_time_sec",
    "K_LL_assembly_time_sec",
    "K_GL_assembly_time_sec",
    "K_total_assembly_time_sec",
    "M_GL_assembly_time_sec",
    "matrix_generation_time_sec",
    "solver_time_sec",
    "stress_calculation_time_sec",
    "output_time_sec",
    "total_time_sec",
]

REQUESTED_TIME_FIELDS = [
    "input_time_sec",
    "matrix_generation_time_sec",
    "solver_time_sec",
    "output_time_sec",
]

CORE_TIME_FIELDS = REQUESTED_TIME_FIELDS + ["total_time_sec"]

TIME_PATTERNS = {
    "input_time_sec": re.compile(r"- input elapse time:\s*([0-9.Ee+-]+)"),
    "nonzero_detection_time_sec": re.compile(
        r"- nonzero detection elapse time:\s*([0-9.Ee+-]+)"
    ),
    "K_GG_assembly_time_sec": re.compile(
        r"- K_GG assembly elapse time:\s*([0-9.Ee+-]+)"
    ),
    "K_LL_assembly_time_sec": re.compile(
        r"- K_LL assembly elapse time:\s*([0-9.Ee+-]+)"
    ),
    "K_GL_assembly_time_sec": re.compile(
        r"- K_GL assembly elapse time:\s*([0-9.Ee+-]+)"
    ),
    "K_total_assembly_time_sec": re.compile(
        r"- K_total assembly elapse time:\s*([0-9.Ee+-]+)"
    ),
    "M_GL_assembly_time_sec": re.compile(
        r"- M_GL assembly elapse time:\s*([0-9.Ee+-]+)"
    ),
    "matrix_generation_time_sec": re.compile(
        r"- matrix generation elapse time:\s*([0-9.Ee+-]+)"
    ),
    "solver_time_sec": re.compile(r"- solver elapse time:\s*([0-9.Ee+-]+)"),
    "stress_calculation_time_sec": re.compile(
        r"- stress calculation elapse time:\s*([0-9.Ee+-]+)"
    ),
    "output_time_sec": re.compile(r"- output elapse time:\s*([0-9.Ee+-]+)"),
    "total_time_sec": re.compile(r"- total\s+elapse time:\s*([0-9.Ee+-]+)"),
}

STEP_PATTERN = re.compile(r"\* current time step:\s*([0-9]+)\s+([0-9.Ee+-]+)")
SOLVER_ITER_PATTERN = re.compile(r"- monolis converge iter\s*:\s*([0-9]+)")
SOLVER_RESIDUAL_PATTERN = re.compile(
    r"- monolis converge residual:\s*([0-9.Ee+-]+)"
)

STEP_TIMING_FIELDS = [
    "step",
    "step_name",
    "step_kind",
    "status",
    "analysis_log",
    "current_time_step",
    "physical_time",
    *TIME_FIELDS,
    "requested_phase_time_sec",
    "other_time_sec",
    "monolis_converge_iter",
    "monolis_converge_residual",
    "missing_fields",
]

TIMING_SUMMARY_FIELDS = [
    "alpha",
    "alpha_label",
    "source_case",
    "output_case",
    "step_start",
    "step_end",
    "steps_expected",
    "steps_found",
    "steps_missing",
    "wall_run_seconds",
    "postprocess_seconds",
    "log_total_time_sec",
    "requested_phase_time_sec",
    "input_time_sec",
    "matrix_generation_time_sec",
    "solver_time_sec",
    "output_time_sec",
    "nonzero_detection_time_sec",
    "stress_calculation_time_sec",
    "other_time_sec",
    "static_total_time_sec",
    "dynamic_total_time_sec",
    "timing_by_step",
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


def make_velocity_baseline_config(velocity=DEFAULT_VELOCITY):
    mesh_constants = load_current_dynamic_constants()
    hL = float(mesh_constants["hL"])
    hG = calculate_actual_hG(hL, BASE_RGL, mesh_constants)
    theta_info = calculate_uniform_theta(hL, mesh_constants["crack_radius"])
    v_label = format_velocity(velocity)

    case = SweepCase(
        idx=1,
        group="velocity",
        label=f"v={v_label}",
        v=float(velocity),
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

    config = case.run_config(
        step_start=0,
        step_end=101,
        command=["generated_from_velocity_sweep_baseline"],
        hL=hL,
    )
    config.update(
        {
            "source_case_generation": "velocity_sweep_baseline_fallback",
            "source_case_missing": True,
            "baseline_source": (
                "param_sweep_dynamic baseline: "
                f"rGL={BASE_RGL}, aL={baseline_aL()}, "
                f"lL={baseline_lL()}, HL={baseline_HL()}"
            ),
        }
    )
    return config


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


def format_seconds(value):
    if value == "":
        return ""
    return f"{float(value):.6f}"


def sum_numeric(rows, field):
    total = 0.0
    for row in rows:
        value = row.get(field, "")
        if value == "":
            continue
        total += float(value)
    return total


def relative_path(path, base=SCRIPT_DIR):
    path = Path(path)
    try:
        return str(path.resolve().relative_to(Path(base).resolve()))
    except ValueError:
        return str(path.resolve())


def parse_analysis_log(log_path):
    row = {field: "" for field in TIME_FIELDS}
    row.update(
        {
            "status": "missing_log",
            "current_time_step": "",
            "physical_time": "",
            "requested_phase_time_sec": "",
            "other_time_sec": "",
            "monolis_converge_iter": "",
            "monolis_converge_residual": "",
            "missing_fields": ",".join(CORE_TIME_FIELDS),
        }
    )

    if not log_path.exists():
        return row
    if log_path.stat().st_size == 0:
        row["status"] = "empty_log"
        return row

    text = log_path.read_text(errors="replace")
    missing = []
    for field, pattern in TIME_PATTERNS.items():
        values = [float(value) for value in pattern.findall(text)]
        if values:
            row[field] = format_seconds(sum(values))
        elif field in CORE_TIME_FIELDS:
            missing.append(field)

    step_matches = STEP_PATTERN.findall(text)
    if step_matches:
        step_id, physical_time = step_matches[-1]
        row["current_time_step"] = step_id
        row["physical_time"] = format_seconds(physical_time)

    iter_matches = SOLVER_ITER_PATTERN.findall(text)
    if iter_matches:
        row["monolis_converge_iter"] = iter_matches[-1]

    residual_matches = SOLVER_RESIDUAL_PATTERN.findall(text)
    if residual_matches:
        row["monolis_converge_residual"] = format_seconds(residual_matches[-1])

    requested = sum_numeric([row], "input_time_sec")
    requested += sum_numeric([row], "matrix_generation_time_sec")
    requested += sum_numeric([row], "solver_time_sec")
    requested += sum_numeric([row], "output_time_sec")
    row["requested_phase_time_sec"] = format_seconds(requested)

    if row["total_time_sec"] != "":
        other = float(row["total_time_sec"]) - requested
        row["other_time_sec"] = format_seconds(other)

    row["missing_fields"] = ",".join(missing)
    row["status"] = "ok" if not missing else "incomplete_log"
    return row


def collect_step_timings(case_folder, step_start, step_end):
    rows = []
    for step in range(step_start, step_end):
        name = step_name(step)
        log_path = Path(case_folder) / name / "log" / "analysis.log"
        row = parse_analysis_log(log_path)
        row.update(
            {
                "step": step,
                "step_name": name,
                "step_kind": "static" if step == 0 else "dynamic",
                "analysis_log": relative_path(log_path),
            }
        )
        rows.append(row)
    return rows


def write_step_timings(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=STEP_TIMING_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def timing_summary_row(
    rows,
    alpha,
    args,
    output_case,
    step_csv,
    run_seconds,
    post_seconds,
):
    steps_found = sum(
        1 for row in rows if row.get("status") not in {"missing_log", "empty_log"}
    )
    steps_missing = sum(
        1 for row in rows if row.get("status") in {"missing_log", "empty_log"}
    )
    log_total = sum_numeric(rows, "total_time_sec")
    requested = sum_numeric(rows, "requested_phase_time_sec")
    static_total = sum_numeric(
        [row for row in rows if row.get("step_kind") == "static"],
        "total_time_sec",
    )
    dynamic_total = sum_numeric(
        [row for row in rows if row.get("step_kind") == "dynamic"],
        "total_time_sec",
    )

    return {
        "alpha": f"{float(alpha):.12g}",
        "alpha_label": format_alpha_for_path(alpha),
        "source_case": str(args.source_case.resolve()),
        "output_case": str(Path(output_case).resolve()),
        "step_start": args.step_start,
        "step_end": args.step_end,
        "steps_expected": args.step_end - args.step_start,
        "steps_found": steps_found,
        "steps_missing": steps_missing,
        "wall_run_seconds": format_seconds(run_seconds),
        "postprocess_seconds": format_seconds(post_seconds),
        "log_total_time_sec": format_seconds(log_total),
        "requested_phase_time_sec": format_seconds(requested),
        "input_time_sec": format_seconds(sum_numeric(rows, "input_time_sec")),
        "matrix_generation_time_sec": format_seconds(
            sum_numeric(rows, "matrix_generation_time_sec")
        ),
        "solver_time_sec": format_seconds(sum_numeric(rows, "solver_time_sec")),
        "output_time_sec": format_seconds(sum_numeric(rows, "output_time_sec")),
        "nonzero_detection_time_sec": format_seconds(
            sum_numeric(rows, "nonzero_detection_time_sec")
        ),
        "stress_calculation_time_sec": format_seconds(
            sum_numeric(rows, "stress_calculation_time_sec")
        ),
        "other_time_sec": format_seconds(sum_numeric(rows, "other_time_sec")),
        "static_total_time_sec": format_seconds(static_total),
        "dynamic_total_time_sec": format_seconds(dynamic_total),
        "timing_by_step": relative_path(step_csv),
    }


def write_timing_summary(row, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=TIMING_SUMMARY_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerow(row)


def collect_and_write_timings(output_case, alpha, args, run_seconds, post_seconds):
    output_case = Path(output_case)
    timing_dir = output_case / "timing"
    step_csv = timing_dir / "timing_by_step.csv"
    summary_csv = timing_dir / "timing_summary.csv"

    rows = collect_step_timings(output_case, args.step_start, args.step_end)
    write_step_timings(rows, step_csv)
    summary = timing_summary_row(
        rows,
        alpha,
        args,
        output_case,
        step_csv,
        run_seconds,
        post_seconds,
    )
    write_timing_summary(summary, summary_csv)
    return summary_csv, step_csv, summary


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
            row = summary_row(
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
            return add_timing_to_summary_row(
                row,
                output_case,
                alpha,
                args,
                run_seconds,
                post_seconds,
            )
        except Exception as exc:
            row = summary_row(
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
            return add_timing_to_summary_row(
                row,
                output_case,
                alpha,
                args,
                run_seconds,
                post_seconds,
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
            log_file.write(f"Source config origin: {args.source_config_origin}\n")
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

        row = summary_row(
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
        return add_timing_to_summary_row(
            row,
            output_case,
            alpha,
            args,
            run_seconds,
            post_seconds,
        )
    finally:
        const_editor.restore()


def add_timing_to_summary_row(row, output_case, alpha, args, run_seconds, post_seconds):
    if not args.collect_timing:
        return row

    try:
        summary_csv, step_csv, timing = collect_and_write_timings(
            output_case,
            alpha,
            args,
            run_seconds,
            post_seconds,
        )
    except Exception as exc:
        row["message"] = f"{row['message']}; timing collection failed: {exc}"
        return row

    row.update(
        {
            "timing_summary": relative_path(summary_csv),
            "timing_by_step": relative_path(step_csv),
            "timing_steps_found": timing["steps_found"],
            "timing_steps_missing": timing["steps_missing"],
            "timing_log_total_time_sec": timing["log_total_time_sec"],
            "timing_requested_phase_time_sec": timing["requested_phase_time_sec"],
            "timing_input_time_sec": timing["input_time_sec"],
            "timing_matrix_generation_time_sec": timing["matrix_generation_time_sec"],
            "timing_solver_time_sec": timing["solver_time_sec"],
            "timing_output_time_sec": timing["output_time_sec"],
            "timing_other_time_sec": timing["other_time_sec"],
        }
    )
    return row


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
        "timing_summary": "",
        "timing_by_step": "",
        "timing_steps_found": "",
        "timing_steps_missing": "",
        "timing_log_total_time_sec": "",
        "timing_requested_phase_time_sec": "",
        "timing_input_time_sec": "",
        "timing_matrix_generation_time_sec": "",
        "timing_solver_time_sec": "",
        "timing_output_time_sec": "",
        "timing_other_time_sec": "",
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
        help=(
            "Existing reference case folder whose run_config.json defines the case. "
            "If the default v500 folder is missing, the script generates the same "
            "velocity-sweep baseline config from the current const files."
        ),
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
    parser.add_argument(
        "--collect-timing",
        action="store_true",
        dest="collect_timing",
        help="Write timing/timing_by_step.csv and timing/timing_summary.csv for each alpha result.",
    )
    parser.add_argument(
        "--no-collect-timing",
        action="store_false",
        dest="collect_timing",
        help="Do not parse per-step analysis.log files after each run.",
    )
    add_j_integral_args(parser)
    parser.set_defaults(postprocess=True, collect_timing=True)
    args = parser.parse_args()
    args.source_case_explicit = any(
        arg == "--source-case" or arg.startswith("--source-case=")
        for arg in sys.argv[1:]
    )
    return args


def resolve_args(args):
    args.source_case = args.source_case.resolve()
    args.output_root = args.output_root.resolve()
    config_path = args.source_case / "run_config.json"
    args.source_config_origin = "source_case"

    if config_path.exists():
        source_config = load_json(config_path)
    elif args.source_case_explicit:
        if not args.source_case.is_dir():
            raise FileNotFoundError(f"Source case folder not found: {args.source_case}")
        raise FileNotFoundError(f"Missing source run_config.json: {config_path}")
    else:
        source_config = make_velocity_baseline_config(DEFAULT_VELOCITY)
        args.source_config_origin = "generated_velocity_baseline"

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
    print(f"Source config: {args.source_config_origin}")
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
    print(f"Collect timing: {args.collect_timing}")
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
