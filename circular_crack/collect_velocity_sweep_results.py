#!/usr/bin/env python3
"""
Collect post-processed crack-velocity sweep CSV files.

By default this script reads:

    results/crack_velocity_sweep_dynamic

and writes four CSV files under:

    results/crack_velocity_sweep_dynamic/collected_postprocess

Outputs:

1. local_stress_norm_ave.csv
2. local_stress_norm_max_min.csv
3. dsif_norm_ave.csv
4. dsif_norm_max_min.csv
"""

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RESULTS_DIR = SCRIPT_DIR / "results" / "crack_velocity_sweep_dynamic"
DEFAULT_OUTPUT_DIRNAME = "collected_postprocess"
DEFAULT_VELOCITY_START = 200
DEFAULT_VELOCITY_STOP = 1500
DEFAULT_VELOCITY_STEP = 100

SUMMARY_FILE = "crack_velocity_sweep_summary.csv"
PROFILE_FILE = "local_stress_normalized_profile.csv"
DSIF_FILE = "dsif_normalized.csv"
TARGET_LOCAL_STRESS_DISTANCE_HL = 5.0


@dataclass
class VelocityResult:
    velocity: float
    label: str
    folder: Path | None
    postprocess_dir: Path | None
    profile_rows: list[dict[str, str]]
    dsif_rows: list[dict[str, str]]


def format_velocity(value):
    value_float = float(value)
    if value_float.is_integer():
        return str(int(value_float))
    return f"{value_float:g}".replace(".", "p")


def velocity_header(value):
    return f"V{format_velocity(value)}"


def parse_velocity_label(text):
    if text is None:
        return None
    match = re.search(r"[-+]?\d+(?:\.\d+)?", str(text))
    if match is None:
        return None
    return format_velocity(match.group(0))


def default_velocities():
    return list(
        range(
            DEFAULT_VELOCITY_START,
            DEFAULT_VELOCITY_STOP + 1,
            DEFAULT_VELOCITY_STEP,
        )
    )


def read_csv_dicts(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))


def read_summary(results_dir):
    summary_path = results_dir / SUMMARY_FILE
    if not summary_path.exists():
        return {}, None

    rows = read_csv_dicts(summary_path)
    by_velocity = {}
    for row in rows:
        label = (
            parse_velocity_label(row.get("v"))
            or parse_velocity_label(row.get("velocity"))
            or parse_velocity_label(row.get("label"))
        )
        if label is not None:
            by_velocity[label] = row
    return by_velocity, summary_path


def unique_paths(paths):
    seen = set()
    unique = []
    for path in paths:
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        unique.append(path)
    return unique


def folder_candidates_from_summary(row, results_dir, v_label):
    folder_text = (row or {}).get("folder", "").strip()
    if not folder_text:
        return []

    raw = Path(folder_text)
    candidates = [raw]
    v_dirname = f"v{v_label}"

    if raw.is_absolute():
        parts = raw.parts
        if v_dirname in parts:
            v_index = parts.index(v_dirname)
            candidates.append(results_dir.joinpath(*parts[v_index:]))
        candidates.append(results_dir / v_dirname / raw.name)
    else:
        candidates.append(results_dir / raw)
        candidates.append(results_dir / v_dirname / raw)

    return unique_paths(candidates)


def postprocess_dir_for(folder, step=None):
    postprocess_dir = folder / "postprocess"
    if step is not None:
        postprocess_dir = postprocess_dir / f"step{int(step):05d}"
    return postprocess_dir


def postprocess_score(folder, step=None):
    postprocess_dir = postprocess_dir_for(folder, step=step)
    score = 0
    if (postprocess_dir / PROFILE_FILE).exists():
        score += 1
    if (postprocess_dir / DSIF_FILE).exists():
        score += 1
    return score


def choose_existing_folder(candidates, step=None):
    existing = [path for path in candidates if path.is_dir()]
    if not existing:
        return None
    return max(existing, key=lambda path: (postprocess_score(path, step=step), str(path)))


def discover_case_folder(results_dir, v_label, step=None):
    v_dir = results_dir / f"v{v_label}"
    if not v_dir.is_dir():
        return None

    candidates = []
    if (v_dir / "postprocess").is_dir():
        candidates.append(v_dir)

    for child in sorted(v_dir.iterdir()):
        if child.is_dir():
            candidates.append(child)

    return choose_existing_folder(candidates, step=step)


def find_case_folder(results_dir, v_label, summary_row, step=None):
    folder = choose_existing_folder(
        folder_candidates_from_summary(summary_row, results_dir, v_label),
        step=step,
    )
    if folder is not None:
        return folder
    return discover_case_folder(results_dir, v_label, step=step)


def find_postprocess_dir(folder, step=None):
    postprocess_dir = postprocess_dir_for(folder, step=step)
    if step is not None:
        return postprocess_dir

    # If the run was post-processed for exactly one non-default step, accept it.
    root = folder / "postprocess"
    if root.is_dir() and (
        (root / PROFILE_FILE).exists() or (root / DSIF_FILE).exists()
    ):
        return root

    if root.is_dir():
        step_dirs = [
            child
            for child in root.iterdir()
            if child.is_dir() and re.fullmatch(r"step\d{5}", child.name)
        ]
        scored = [
            (int((child / PROFILE_FILE).exists()) + int((child / DSIF_FILE).exists()), child)
            for child in step_dirs
        ]
        scored = [item for item in scored if item[0] > 0]
        if len(scored) == 1:
            return scored[0][1]

    return postprocess_dir


def load_velocity_result(results_dir, velocity, summary_row=None, step=None):
    v_label = format_velocity(velocity)
    folder = find_case_folder(results_dir, v_label, summary_row, step=step)
    if folder is None:
        return VelocityResult(
            velocity=float(velocity),
            label=velocity_header(velocity),
            folder=None,
            postprocess_dir=None,
            profile_rows=[],
            dsif_rows=[],
        )

    postprocess_dir = find_postprocess_dir(folder, step=step)
    profile_path = postprocess_dir / PROFILE_FILE
    dsif_path = postprocess_dir / DSIF_FILE
    profile_rows = read_csv_dicts(profile_path) if profile_path.exists() else []
    dsif_rows = read_csv_dicts(dsif_path) if dsif_path.exists() else []

    return VelocityResult(
        velocity=float(velocity),
        label=velocity_header(velocity),
        folder=folder,
        postprocess_dir=postprocess_dir,
        profile_rows=profile_rows,
        dsif_rows=dsif_rows,
    )


def csv_value(row, column):
    value = row.get(column, "")
    return "" if value is None else value


def finite_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    return number


def format_number(value):
    if value is None:
        return ""
    return f"{value:.17g}"


def write_single_column_csv(output_path, values):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as f:
        writer = csv.writer(f)
        for value in values:
            writer.writerow([value])


def write_vertical_max_min_csv(output_path, pairs):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as f:
        writer = csv.writer(f)
        for index, (max_value, min_value) in enumerate(pairs):
            if index > 0:
                writer.writerow([])
            writer.writerow([max_value])
            writer.writerow([min_value])


def local_stress_target_profile_row(result):
    for row in result.profile_rows:
        distance = finite_float(row.get("distance_hL"))
        if distance is None:
            continue
        if math.isclose(
            distance,
            TARGET_LOCAL_STRESS_DISTANCE_HL,
            rel_tol=0.0,
            abs_tol=1e-9,
        ):
            return row

    target_index = int(TARGET_LOCAL_STRESS_DISTANCE_HL) - 1
    if target_index >= 0 and target_index < len(result.profile_rows):
        return result.profile_rows[target_index]
    return None


def local_stress_5hL_norm_ave(result):
    row = local_stress_target_profile_row(result)
    if row is None:
        return ""
    return csv_value(row, "norm_ave")


def local_stress_5hL_norm_max_min(result):
    row = local_stress_target_profile_row(result)
    if row is None:
        return "", ""
    return csv_value(row, "norm_max"), csv_value(row, "norm_min")


def dsif_normalized_values(result):
    return [
        value
        for value in (finite_float(row.get("normalized")) for row in result.dsif_rows)
        if value is not None
    ]


def dsif_stats(result):
    values = dsif_normalized_values(result)
    if not values:
        return None, None, None
    return sum(values) / len(values), max(values), min(values)


def write_local_stress_outputs(output_dir, results):
    norm_ave_path = output_dir / "local_stress_norm_ave.csv"
    write_single_column_csv(
        norm_ave_path,
        [local_stress_5hL_norm_ave(result) for result in results],
    )

    max_min_path = output_dir / "local_stress_norm_max_min.csv"
    write_vertical_max_min_csv(
        max_min_path,
        [local_stress_5hL_norm_max_min(result) for result in results],
    )

    return [norm_ave_path, max_min_path]


def write_dsif_outputs(output_dir, results):
    ave_path = output_dir / "dsif_norm_ave.csv"
    max_min_path = output_dir / "dsif_norm_max_min.csv"
    output_dir.mkdir(parents=True, exist_ok=True)

    stats_by_result = {result.label: dsif_stats(result) for result in results}

    write_single_column_csv(
        ave_path,
        [
            format_number(stats_by_result[result.label][0])
            for result in results
        ],
    )

    write_vertical_max_min_csv(
        max_min_path,
        [
            (
                format_number(stats_by_result[result.label][1]),
                format_number(stats_by_result[result.label][2]),
            )
            for result in results
        ],
    )

    return [ave_path, max_min_path]


def collect_warnings(results, summary_path=None):
    warnings = []
    if summary_path is None:
        warnings.append(
            f"Summary CSV not found; discovered folders directly under the results directory."
        )

    for result in results:
        if result.folder is None:
            warnings.append(f"{result.label}: result folder not found")
            continue
        if result.postprocess_dir is None or not result.postprocess_dir.is_dir():
            warnings.append(f"{result.label}: postprocess directory not found")
            continue
        if not result.profile_rows:
            warnings.append(f"{result.label}: missing or empty {result.postprocess_dir / PROFILE_FILE}")
        elif local_stress_target_profile_row(result) is None:
            warnings.append(
                f"{result.label}: no local-stress row for distance_hL={TARGET_LOCAL_STRESS_DISTANCE_HL:g}"
            )
        if not result.dsif_rows:
            warnings.append(f"{result.label}: missing or empty {result.postprocess_dir / DSIF_FILE}")
    return warnings


def parse_args():
    parser = argparse.ArgumentParser(
        description="Collect crack-velocity sweep postprocess CSV files into four summary CSV files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help="Directory containing crack_velocity_sweep_summary.csv and v*/case folders.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Directory for collected CSV files. "
            f"Default: RESULTS_DIR/{DEFAULT_OUTPUT_DIRNAME}"
        ),
    )
    parser.add_argument(
        "--velocities",
        type=float,
        nargs="+",
        default=None,
        help="Velocities to collect. Default is 200..1500 with step 100.",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=None,
        help=(
            "Optional postprocess step. If set, reads postprocess/stepNNNNN. "
            "If omitted, reads the default postprocess directory."
        ),
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return a non-zero exit code when any requested output is missing.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    results_dir = args.results_dir.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else results_dir / DEFAULT_OUTPUT_DIRNAME
    )
    velocities = args.velocities if args.velocities is not None else default_velocities()

    summary_rows, summary_path = read_summary(results_dir)
    results = [
        load_velocity_result(
            results_dir,
            velocity,
            summary_row=summary_rows.get(format_velocity(velocity)),
            step=args.step,
        )
        for velocity in velocities
    ]

    written = []
    written.extend(write_local_stress_outputs(output_dir, results))
    written.extend(write_dsif_outputs(output_dir, results))

    warnings = collect_warnings(results, summary_path=summary_path)

    print(f"Collected velocity-sweep CSV files written to: {output_dir}")
    for path in written:
        print(f"  - {path}")

    if warnings:
        print("\nWarnings:")
        for warning in warnings:
            print(f"  - {warning}")
        if args.strict:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
