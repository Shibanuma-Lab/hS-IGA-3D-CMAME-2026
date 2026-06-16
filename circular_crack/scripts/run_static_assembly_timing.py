#!/usr/bin/env python3
"""
Run static assembly-only timing for verification_5_2 cases.

The script reuses saved input files from results/verification_5_2, runs the
Fortran solver with SFEM_ASSEMBLY_ONLY=1, and stores logs plus a summary CSV.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
CIRCULAR_DIR = SCRIPT_DIR.parent
PROJECT_ROOT = CIRCULAR_DIR.parent
SFEM_DIR = PROJECT_ROOT / "sfem_linear"
EXECUTABLE = SFEM_DIR / "bin" / "sfem_linear"
RESULTS_DIR = CIRCULAR_DIR / "results" / "verification_5_2"
OUTPUT_DIR = RESULTS_DIR / "assembly_timing"
TMP_ROOT = PROJECT_ROOT / ".tmp_static_assembly_timing"

PAPER_CASE_DOFS = {
    4: {
        46260, 61776, 80424, 102492, 129402, 159345, 193584, 234090,
        277992, 327066, 381600, 444465, 511038, 583947, 666864, 753600,
        847548, 948996,
    },
    8: {
        253626, 345465, 454059, 586992, 739332, 830388, 920943, 1020741,
        1124421, 1244214, 1362294, 1491297, 1631547, 1776768,
    },
}


TIME_PATTERNS = {
    "K_GG_time_sec": re.compile(r"- K_GG assembly elapse time:\s*([0-9.Ee+-]+)"),
    "K_LL_time_sec": re.compile(r"- K_LL assembly elapse time:\s*([0-9.Ee+-]+)"),
    "K_GL_time_sec": re.compile(r"- K_GL assembly elapse time:\s*([0-9.Ee+-]+)"),
    "K_total_time_sec": re.compile(r"- K_total assembly elapse time:\s*([0-9.Ee+-]+)"),
    "matrix_generation_time_sec": re.compile(r"- matrix generation elapse time:\s*([0-9.Ee+-]+)"),
    "input_time_sec": re.compile(r"- input elapse time:\s*([0-9.Ee+-]+)"),
    "nonzero_detection_time_sec": re.compile(r"- nonzero detection elapse time:\s*([0-9.Ee+-]+)"),
    "total_time_sec": re.compile(r"- total\s+elapse time:\s*([0-9.Ee+-]+)"),
}


def format_case_dir(rgl: int, h_l: float, h_g: float) -> Path:
    return RESULTS_DIR / f"rGL{rgl}_0.25" / f"hL_{h_l:.6f}_hG_{h_g:.6f}"


def step_string(h_l: float) -> str:
    return f"{int(round(1.0 / h_l)):05d}"


def read_l2_rows(rgl: int) -> list[dict[str, str]]:
    csv_path = RESULTS_DIR / f"L2_norm_gpu_rGL{rgl}.csv"
    with csv_path.open(newline="") as file:
        return list(csv.DictReader(file))


def read_gauss_counts(rgl: int) -> dict[int, int]:
    csv_path = RESULTS_DIR / f"rGL{rgl}_0.25" / "gauss_count.csv"
    counts: dict[int, int] = {}
    if not csv_path.exists():
        return counts

    with csv_path.open(newline="") as file:
        for row in csv.reader(file):
            if len(row) < 2:
                continue
            counts[int(float(row[0]))] = int(float(row[1]))
    return counts


def parse_analysis_log(log_path: Path) -> dict[str, str]:
    values = {key: "" for key in TIME_PATTERNS}
    if not log_path.exists():
        return values

    text = log_path.read_text(errors="replace")
    for key, pattern in TIME_PATTERNS.items():
        match = pattern.search(text)
        if match:
            values[key] = match.group(1)
    return values


def copy_input(case_dir: Path, h_l: float, run_dir: Path) -> Path:
    step = step_string(h_l)
    input_dir = case_dir / f"inputfiles_step{step}"
    if not input_dir.exists():
        input_dir = case_dir / f"step{step}"
    if not input_dir.exists():
        raise FileNotFoundError(f"input directory not found for {case_dir}")

    if run_dir.exists():
        shutil.rmtree(run_dir)
    shutil.copytree(input_dir, run_dir)
    return input_dir


def run_case(row: dict[str, str], threads: int, force: bool) -> dict[str, str]:
    rgl = int(float(row["rGL"]))
    h_l = float(row["hL"])
    h_g = float(row["hG"])
    dof = int(float(row["dof"]))
    step = step_string(h_l)

    case_dir = format_case_dir(rgl, h_l, h_g)
    label = f"rGL{rgl}_hL_{h_l:.6f}_threads_{threads}"
    run_dir = TMP_ROOT / label
    log_dir = OUTPUT_DIR / f"threads_{threads}" / f"rGL{rgl}" / f"hL_{h_l:.6f}_hG_{h_g:.6f}"
    analysis_log = log_dir / "analysis.log"
    stdout_log = log_dir / "solver.stdout.log"

    result = {
        "rGL": str(rgl),
        "hL": f"{h_l:.16g}",
        "hG": f"{h_g:.16g}",
        "step": step,
        "dof": str(dof),
        "relative_L2_norm": row.get("relative_L2_norm", ""),
        "gauss_points": "",
        "threads": str(threads),
        "case_dir": str(case_dir.relative_to(PROJECT_ROOT)),
        "analysis_log": str(analysis_log.relative_to(PROJECT_ROOT)),
        "status": "skipped",
        "returncode": "",
    }

    if analysis_log.exists() and not force:
        result["status"] = "ok_existing"
        result["returncode"] = "0"
        result.update(parse_analysis_log(analysis_log))
        return result

    if not EXECUTABLE.exists():
        raise FileNotFoundError(f"solver executable not found: {EXECUTABLE}")

    source_input = copy_input(case_dir, h_l, run_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["SFEM_ASSEMBLY_ONLY"] = "1"
    env["SFEM_OMP_THREADS"] = str(threads)
    env["OMP_NUM_THREADS"] = str(threads)

    with stdout_log.open("w") as stdout_file:
        completed = subprocess.run(
            [str(EXECUTABLE), "input.dat"],
            cwd=run_dir,
            env=env,
            stdout=stdout_file,
            stderr=subprocess.STDOUT,
            check=False,
        )

    produced_log = run_dir / "log" / "analysis.log"
    if produced_log.exists():
        shutil.copy2(produced_log, analysis_log)

    result["status"] = "ok" if completed.returncode == 0 else "failed"
    result["returncode"] = str(completed.returncode)
    result["source_input"] = str(source_input.relative_to(PROJECT_ROOT))
    result.update(parse_analysis_log(analysis_log))
    return result


def write_summary(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "rGL",
        "hL",
        "hG",
        "step",
        "dof",
        "relative_L2_norm",
        "gauss_points",
        "threads",
        "K_GG_time_sec",
        "K_LL_time_sec",
        "K_GL_time_sec",
        "K_total_time_sec",
        "matrix_generation_time_sec",
        "input_time_sec",
        "nonzero_detection_time_sec",
        "total_time_sec",
        "status",
        "returncode",
        "case_dir",
        "source_input",
        "analysis_log",
    ]

    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rgl", type=int, nargs="+", default=[4, 8])
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 24])
    parser.add_argument("--case-set", choices=["all", "paper"], default="all")
    parser.add_argument("--output-subdir", default="assembly_timing")
    parser.add_argument("--force", action="store_true", help="rerun cases even when logs already exist")
    parser.add_argument("--summarize-only", action="store_true", help="only summarize existing logs; do not run solver")
    parser.add_argument("--max-cases", type=int, default=0, help="limit cases per rGL for smoke tests")
    parser.add_argument("--start-index", type=int, default=0, help="zero-based start index per rGL")
    parser.add_argument("--end-index", type=int, default=0, help="exclusive end index per rGL; 0 means no limit")
    return parser.parse_args()


def summarize_existing_case(row: dict[str, str], threads: int) -> dict[str, str]:
    rgl = int(float(row["rGL"]))
    h_l = float(row["hL"])
    h_g = float(row["hG"])
    dof = int(float(row["dof"]))
    step = step_string(h_l)
    case_dir = format_case_dir(rgl, h_l, h_g)
    analysis_log = (
        OUTPUT_DIR
        / f"threads_{threads}"
        / f"rGL{rgl}"
        / f"hL_{h_l:.6f}_hG_{h_g:.6f}"
        / "analysis.log"
    )

    result = {
        "rGL": str(rgl),
        "hL": f"{h_l:.16g}",
        "hG": f"{h_g:.16g}",
        "step": step,
        "dof": str(dof),
        "relative_L2_norm": row.get("relative_L2_norm", ""),
        "gauss_points": "",
        "threads": str(threads),
        "case_dir": str(case_dir.relative_to(PROJECT_ROOT)),
        "source_input": "",
        "analysis_log": str(analysis_log.relative_to(PROJECT_ROOT)),
        "status": "missing",
        "returncode": "",
    }
    result.update(parse_analysis_log(analysis_log))
    if analysis_log.exists():
        result["status"] = "ok_existing" if result.get("K_total_time_sec") else "incomplete_log"
        result["returncode"] = "0" if result["status"] == "ok_existing" else ""
    return result


def main() -> int:
    global OUTPUT_DIR

    args = parse_args()
    OUTPUT_DIR = RESULTS_DIR / args.output_subdir
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    TMP_ROOT.mkdir(parents=True, exist_ok=True)

    all_results: list[dict[str, str]] = []
    for rgl in args.rgl:
        l2_rows = read_l2_rows(rgl)
        gauss_counts = read_gauss_counts(rgl)

        start = args.start_index
        end = args.end_index if args.end_index > 0 else len(l2_rows)
        rows = l2_rows[start:end]
        if args.case_set == "paper":
            wanted_dofs = PAPER_CASE_DOFS.get(rgl, set())
            rows = [row for row in rows if int(float(row["dof"])) in wanted_dofs]
        if args.max_cases > 0:
            rows = rows[: args.max_cases]

        for row in rows:
            dof = int(float(row["dof"]))
            for threads in args.threads:
                if args.summarize_only:
                    result = summarize_existing_case(row, threads)
                else:
                    print(
                        f"Running rGL={rgl}, hL={float(row['hL']):.6f}, "
                        f"threads={threads}",
                        flush=True,
                    )
                    result = run_case(row, threads, args.force)
                result["gauss_points"] = str(gauss_counts.get(dof, ""))
                all_results.append(result)

    summary_path = OUTPUT_DIR / "static_assembly_timing_summary.csv"
    write_summary(all_results, summary_path)
    completed_path = OUTPUT_DIR / "static_assembly_timing_completed.csv"
    completed_rows = [
        row for row in all_results
        if row.get("status", "").startswith("ok") and row.get("K_total_time_sec")
    ]
    write_summary(completed_rows, completed_path)
    print(f"Summary written to {summary_path}")
    print(f"Completed-case summary written to {completed_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
