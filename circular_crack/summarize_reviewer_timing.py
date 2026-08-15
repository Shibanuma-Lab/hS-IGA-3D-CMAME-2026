#!/usr/bin/env python3
"""Create reviewer-facing timing tables from hierarchical hS-IGA timing CSVs."""

import argparse
import csv
from decimal import Decimal
from pathlib import Path


CUMULATIVE_TARGETS = (20, 40, 60, 80, 100)


def read_rows(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))


def write_rows(path, fieldnames, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def seconds(row):
    return float(row["time_sec"])


def minutes(value):
    return f"{float(value) / 60.0:.9f}"


def fraction(value, parent):
    return f"{float(value) / float(parent):.9f}" if parent else ""


def percent(value, parent):
    return f"{100.0 * float(value) / float(parent):.6f}" if parent else ""


def reconcile_displayed_children(rows, parent_phase_id, residual_phase_id):
    parent = next(row for row in rows if row["phase_id"] == parent_phase_id)
    children = [
        row for row in rows if row["parent_phase_id"] == parent_phase_id
    ]
    residual = next(row for row in children if row["phase_id"] == residual_phase_id)
    other_children = [
        row for row in children if row["phase_id"] != residual_phase_id
    ]

    residual["time_min"] = (
        f"{Decimal(parent['time_min']) - sum(Decimal(row['time_min']) for row in other_children):.9f}"
    )
    if "fraction_of_parent" in residual:
        residual["fraction_of_parent"] = (
            f"{Decimal('1') - sum(Decimal(row['fraction_of_parent']) for row in other_children):.9f}"
        )
        residual["percent_of_parent"] = (
            f"{Decimal('100') - sum(Decimal(row['percent_of_parent']) for row in other_children):.6f}"
        )
    if parent_phase_id == "fortran_total":
        residual["fraction_of_total"] = (
            f"{Decimal('1') - sum(Decimal(row['fraction_of_total']) for row in other_children):.9f}"
        )
        residual["percent_of_total"] = (
            f"{Decimal('100') - sum(Decimal(row['percent_of_total']) for row in other_children):.6f}"
        )


def unique_phase(rows, tree_id, phase_id, target_step=None):
    matches = [
        row
        for row in rows
        if row["tree_id"] == tree_id
        and row["phase_id"] == phase_id
        and (
            target_step is None
            or int(row["target_step"]) == int(target_step)
        )
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one {tree_id}/{phase_id} row"
            f" for target_step={target_step}; found {len(matches)}"
        )
    return matches[0]


def make_cumulative_table(cumulative_rows):
    output = []
    previous = 0.0
    for target in CUMULATIVE_TARGETS:
        root = unique_phase(
            cumulative_rows,
            "wall_clock",
            "fortran_total",
            target,
        )
        total = seconds(root)
        output.append(
            {
                "target_step": target,
                "steps_included": target + 1,
                "hSIGA_cumulative_fortran_time_min": minutes(total),
                "increment_since_previous_target_min": minutes(total - previous),
                "timing_basis": "measured_wall_clock",
                "source_phase_id": "fortran_total",
                "notes": (
                    f"Includes step 0 through step {target}; excludes Python "
                    "orchestration and post-processing."
                ),
            }
        )
        previous = total
    return output


def make_coarse_table(summary_rows):
    total = seconds(unique_phase(summary_rows, "wall_clock", "fortran_total"))
    values = {
        "input": seconds(unique_phase(summary_rows, "wall_clock", "input")),
        "matrix_generation": seconds(
            unique_phase(summary_rows, "wall_clock", "matrix_generation")
        ),
        "solver": seconds(unique_phase(summary_rows, "wall_clock", "solver")),
        "output": seconds(unique_phase(summary_rows, "wall_clock", "output")),
    }
    values["other"] = total - sum(values.values())

    phases = [
        (
            "fortran_total",
            "",
            "Fortran total",
            total,
            "measured_wall_clock",
            "Sum of analysis.log total elapsed times for step 0 through step 100.",
        ),
        (
            "input",
            "fortran_total",
            "Input",
            values["input"],
            "measured_wall_clock",
            "",
        ),
        (
            "matrix_generation",
            "fortran_total",
            "Matrix generation",
            values["matrix_generation"],
            "measured_wall_clock",
            "",
        ),
        (
            "solver",
            "fortran_total",
            "Linear solver",
            values["solver"],
            "measured_wall_clock",
            "",
        ),
        (
            "output",
            "fortran_total",
            "Output",
            values["output"],
            "measured_wall_clock",
            "",
        ),
        (
            "other",
            "fortran_total",
            "Other",
            values["other"],
            "derived_remainder",
            (
                "Fortran total minus input, matrix generation, solver and output; "
                "includes nonzero detection, stress calculation and residual overhead."
            ),
        ),
    ]

    output = [
        {
            "phase_order": order,
            "level": 0 if phase_id == "fortran_total" else 1,
            "phase_id": phase_id,
            "parent_phase_id": parent_id,
            "phase_label": label,
            "time_min": minutes(value),
            "fraction_of_total": fraction(value, total),
            "percent_of_total": percent(value, total),
            "timing_basis": basis,
            "notes": notes,
        }
        for order, (phase_id, parent_id, label, value, basis, notes)
        in enumerate(phases)
    ]
    reconcile_displayed_children(output, "fortran_total", "other")
    return output


def allocated_coupling_phase(summary_rows, prefix, phase_suffix):
    wall = seconds(
        unique_phase(summary_rows, "wall_clock", f"{prefix}_assembly")
    )
    work_total = seconds(
        unique_phase(
            summary_rows,
            f"{prefix}_thread_work",
            f"{prefix}_thread_work_total",
        )
    )
    work = seconds(
        unique_phase(
            summary_rows,
            f"{prefix}_thread_work",
            f"{prefix}_{phase_suffix}",
        )
    )
    return wall * work / work_total


def make_reviewer_detail_table(summary_rows):
    total = seconds(unique_phase(summary_rows, "wall_clock", "fortran_total"))
    input_time = seconds(unique_phase(summary_rows, "wall_clock", "input"))
    matrix_time = seconds(
        unique_phase(summary_rows, "wall_clock", "matrix_generation")
    )
    solver_time = seconds(unique_phase(summary_rows, "wall_clock", "solver"))
    output_time = seconds(unique_phase(summary_rows, "wall_clock", "output"))
    other_time = total - input_time - matrix_time - solver_time - output_time

    candidate = sum(
        allocated_coupling_phase(summary_rows, prefix, "candidate_search")
        for prefix in ("K_GL", "M_GL")
    )
    mapping = sum(
        allocated_coupling_phase(summary_rows, prefix, "inverse_mapping")
        for prefix in ("K_GL", "M_GL")
    )
    basis = sum(
        allocated_coupling_phase(summary_rows, prefix, "IGA_basis_gradient")
        for prefix in ("K_GL", "M_GL")
    )
    assembly = matrix_time - candidate - mapping - basis
    if assembly < 0.0:
        raise ValueError("Allocated matrix-assembly remainder is negative")

    phases = [
        (
            0,
            0,
            "fortran_total",
            "",
            "Fortran total",
            total,
            total,
            "measured_wall_clock",
            "Step 0 through step 100.",
        ),
        (
            10,
            1,
            "input",
            "fortran_total",
            "Input",
            input_time,
            total,
            "measured_wall_clock",
            "",
        ),
        (
            20,
            1,
            "matrix_generation",
            "fortran_total",
            "Matrix generation",
            matrix_time,
            total,
            "measured_wall_clock",
            "Its four child rows reconcile exactly to this parent before rounding.",
        ),
        (
            21,
            2,
            "candidate_search",
            "matrix_generation",
            "Candidate search/bookkeeping",
            candidate,
            matrix_time,
            "allocated_wall_clock_from_thread_work",
            "K_GL and M_GL merged; coupling wall times allocated using measured thread-work shares.",
        ),
        (
            22,
            2,
            "inverse_mapping",
            "matrix_generation",
            "Inverse mapping (Newton-Raphson)",
            mapping,
            matrix_time,
            "allocated_wall_clock_from_thread_work",
            "K_GL and M_GL merged; includes initial physical-element mapping and IGA Newton-Raphson mapping.",
        ),
        (
            23,
            2,
            "IGA_basis_gradient",
            "matrix_generation",
            "B-spline basis/gradient evaluation",
            basis,
            matrix_time,
            "allocated_wall_clock_from_thread_work",
            "K_GL and M_GL merged; includes basis/derivative, Jacobian and physical-gradient evaluation.",
        ),
        (
            24,
            2,
            "matrix_assembly",
            "matrix_generation",
            "Matrix assembly and remaining matrix work",
            assembly,
            matrix_time,
            "derived_wall_clock_remainder",
            (
                "Matrix-generation remainder after search, mapping and basis allocation; "
                "includes K_GG/K_LL, coupling arithmetic, sparse accumulation, other mass "
                "assembly, mass lumping, dynamic-matrix combination, RHS generation and "
                "unclassified matrix overhead."
            ),
        ),
        (
            30,
            1,
            "solver",
            "fortran_total",
            "Linear solver",
            solver_time,
            total,
            "measured_wall_clock",
            "",
        ),
        (
            40,
            1,
            "output",
            "fortran_total",
            "Output",
            output_time,
            total,
            "measured_wall_clock",
            "",
        ),
        (
            50,
            1,
            "other",
            "fortran_total",
            "Other",
            other_time,
            total,
            "derived_remainder",
            "Includes nonzero detection, stress calculation and residual Fortran overhead.",
        ),
    ]

    output = []
    for order, level, phase_id, parent_id, label, value, parent, basis_name, notes in phases:
        output.append(
            {
                "phase_order": order,
                "level": level,
                "phase_id": phase_id,
                "parent_phase_id": parent_id,
                "phase_label": label,
                "time_min": minutes(value),
                "fraction_of_parent": fraction(value, parent),
                "percent_of_parent": percent(value, parent),
                "fraction_of_total": fraction(value, total),
                "percent_of_total": percent(value, total),
                "timing_basis": basis_name,
                "notes": notes,
            }
        )

    reconcile_displayed_children(output, "matrix_generation", "matrix_assembly")
    reconcile_displayed_children(output, "fortran_total", "other")
    matrix_children = [
        row
        for row in output
        if row["parent_phase_id"] == "matrix_generation"
    ]
    displayed_difference = float(
        next(row for row in output if row["phase_id"] == "matrix_generation")[
            "time_min"
        ]
    ) - sum(float(row["time_min"]) for row in matrix_children)
    if abs(displayed_difference) > 5.0e-9:
        raise ValueError(
            "Displayed matrix hierarchy does not reconcile: "
            f"{displayed_difference:.12e} min"
        )
    return output


def main():
    parser = argparse.ArgumentParser(
        description="Create final reviewer-facing hS-IGA timing CSV tables."
    )
    parser.add_argument(
        "--timing-dir",
        type=Path,
        required=True,
        help="Directory containing timing_summary.csv and timing_cumulative.csv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory; defaults to TIMING_DIR/final_statistics.",
    )
    args = parser.parse_args()

    timing_dir = args.timing_dir.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else timing_dir / "final_statistics"
    )
    summary_rows = read_rows(timing_dir / "timing_summary.csv")
    cumulative_rows = read_rows(timing_dir / "timing_cumulative.csv")

    cumulative = make_cumulative_table(cumulative_rows)
    coarse = make_coarse_table(summary_rows)
    detail = make_reviewer_detail_table(summary_rows)

    write_rows(
        output_dir / "01_hSIGA_cumulative_total_time.csv",
        [
            "target_step",
            "steps_included",
            "hSIGA_cumulative_fortran_time_min",
            "increment_since_previous_target_min",
            "timing_basis",
            "source_phase_id",
            "notes",
        ],
        cumulative,
    )
    write_rows(
        output_dir / "02_hSIGA_coarse_phase_time_101_steps.csv",
        [
            "phase_order",
            "level",
            "phase_id",
            "parent_phase_id",
            "phase_label",
            "time_min",
            "fraction_of_total",
            "percent_of_total",
            "timing_basis",
            "notes",
        ],
        coarse,
    )
    write_rows(
        output_dir / "03_hSIGA_reviewer_detailed_time_101_steps.csv",
        [
            "phase_order",
            "level",
            "phase_id",
            "parent_phase_id",
            "phase_label",
            "time_min",
            "fraction_of_parent",
            "percent_of_parent",
            "fraction_of_total",
            "percent_of_total",
            "timing_basis",
            "notes",
        ],
        detail,
    )
    print(output_dir)


if __name__ == "__main__":
    main()
