#!/usr/bin/env python3
"""
Post-process one dynamic parameter-sweep case.

Outputs:
- local_stress_cal.csv
- local_stress_normalized_profile.csv
- local_stress_ring_5hL.csv
- hsiga_J.csv / hsiga_J_KId.csv
- fem_reference_J.csv / fem_reference_J_KId.csv
- dsif_normalized.csv
"""

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from const import const_jintegral as cji
from const import const_local_mesh as clm
from const import material_property as mp
from const import simulation_params as sp
from jintegral import FEMReferenceJIntegralCalculator, JIntegralCalculator

FEM_LOCAL_STRESS_DIR = SCRIPT_DIR / "data" / "FEM_local_stress"
FEM_LOCAL_STRESS_PROFILE_STEP = 100


@dataclass
class CaseConfig:
    velocity: float
    rGL: int
    aL: int
    lL: int
    HL: int
    d_theta: float
    hL: float
    step_start: int
    step_end: int

    @property
    def final_step(self):
        return self.step_end - 1


@dataclass
class FEMLocalStressReference:
    profile_distance_mm: Optional[np.ndarray]
    profile_stress: Optional[np.ndarray]
    profile_source: Optional[str]
    ring5_distance_mm: float
    ring5_stress: float
    ring5_source: str

    @property
    def has_profile(self):
        return self.profile_distance_mm is not None and self.profile_stress is not None


@dataclass(frozen=True)
class JIntegralParams:
    Rj0: float
    Rj1: float
    Wj0: float
    Wj1: float
    ngp: int

    def calculator_kwargs(self):
        return {
            "Rj0": self.Rj0,
            "Rj1": self.Rj1,
            "Wj0": self.Wj0,
            "Wj1": self.Wj1,
            "ngp": self.ngp,
        }


def default_j_integral_params():
    return {
        "Rj0": float(cji.Rj0),
        "Rj1": float(cji.Rj1),
        "Wj0": float(cji.Wj0),
        "Wj1": float(cji.Wj1),
        "ngp": int(cji.ngp),
    }


def resolve_j_integral_params(overrides=None):
    if isinstance(overrides, JIntegralParams):
        params = overrides
    else:
        data = default_j_integral_params()
        if overrides:
            for key in data:
                value = overrides.get(key)
                if value is not None:
                    data[key] = value
        params = JIntegralParams(
            Rj0=float(data["Rj0"]),
            Rj1=float(data["Rj1"]),
            Wj0=float(data["Wj0"]),
            Wj1=float(data["Wj1"]),
            ngp=int(data["ngp"]),
        )

    if params.Rj0 <= 0.0 or params.Rj1 < params.Rj0:
        raise ValueError(
            f"Require 0 < Rj0 <= Rj1, got Rj0={params.Rj0}, Rj1={params.Rj1}"
        )
    if params.Wj0 <= 0.0 or params.Wj1 < params.Wj0:
        raise ValueError(
            f"Require 0 < Wj0 <= Wj1, got Wj0={params.Wj0}, Wj1={params.Wj1}"
        )
    if params.ngp < 1:
        raise ValueError(f"Require ngp >= 1, got {params.ngp}")

    return params


def load_case_config(case_folder):
    config_path = Path(case_folder) / "run_config.json"
    if not config_path.exists():
        raise FileNotFoundError(f"Missing run_config.json: {config_path}")

    data = json.loads(config_path.read_text())
    return CaseConfig(
        velocity=float(data["velocity"]),
        rGL=int(data["rGL"]),
        aL=int(data["aL"]),
        lL=int(data["lL"]),
        HL=int(data["HL"]),
        d_theta=float(data["d_theta"]),
        hL=float(data.get("hL", clm.hL)),
        step_start=int(data.get("step_start", 0)),
        step_end=int(data.get("step_end", 101)),
    )


def read_table_after_count(path, dtype=float):
    return np.loadtxt(path, skiprows=1, dtype=dtype)


def local_stress_from_reaction(case_folder, config, step):
    step_dir = Path(case_folder) / f"step{step:05d}"
    log_dir = step_dir / "log"
    nodeL = read_table_after_count(step_dir / "node.l.dat")[:, 1:4]
    elemL = read_table_after_count(step_dir / "elem.l.dat", dtype=int)[:, 1:]
    rfL = read_table_after_count(log_dir / "reaction.l.dat")[:, 1:4]

    nLr = config.aL + config.lL
    nLtheta = int(round(90.0 / config.d_theta))
    nlrt = (nLr + 1) * (nLtheta + 1)

    nodeLxy = nodeL[:nlrt, :2]
    elemLxy = elemL[: nLr * nLtheta, :4]

    div_area = np.zeros(nlrt)
    for elem in elemLxy:
        pts = nodeLxy[elem - 1]
        vec1 = pts[1] - pts[0]
        vec2 = pts[3] - pts[0]
        vec3 = pts[1] - pts[2]
        vec4 = pts[3] - pts[2]
        area = 0.5 * abs(vec1[0] * vec2[1] - vec1[1] * vec2[0])
        area += 0.5 * abs(vec3[0] * vec4[1] - vec3[1] * vec4[0])
        div_area[elem - 1] += area / 4.0

    rfLz = rfL[:nlrt, 2]
    stress_nodes = np.full(nlrt, np.nan)
    valid = div_area > 0.0
    stress_nodes[valid] = -(rfLz[valid] / div_area[valid])

    stress_grid = stress_nodes.reshape(nLtheta + 1, nLr + 1)
    if step == 0:
        ap_1based = 1
    elif step <= config.aL:
        ap_1based = step + 2
    else:
        ap_1based = config.aL + 2

    stress_theta_distance = stress_grid[:, ap_1based - 1 :]
    stress_distance_theta = stress_theta_distance.T
    theta_deg = np.linspace(0.0, 90.0, nLtheta + 1)
    distance_hL = np.arange(1, stress_distance_theta.shape[0] + 1)

    return distance_hL, theta_deg, stress_distance_theta


def fem_velocity_label(velocity):
    v_label = int(round(float(velocity)))
    return v_label


def fem_local_stress_profile_path(velocity):
    v_label = fem_velocity_label(velocity)
    return FEM_LOCAL_STRESS_DIR / f"FEM_v_{v_label}_last.xlsx"


def fem_local_stress_history_path(velocity):
    v_label = fem_velocity_label(velocity)
    return FEM_LOCAL_STRESS_DIR / f"FEM_v_{v_label}_list_5.xlsx"


def load_fem_local_stress_profile(velocity):
    path = fem_local_stress_profile_path(velocity)
    if not path.exists():
        raise FileNotFoundError(f"Missing FEM local stress reference: {path}")

    frame = pd.read_excel(path, header=None)
    frame = frame.dropna(how="any")
    return (
        frame.iloc[:, 0].to_numpy(dtype=float),
        frame.iloc[:, 1].to_numpy(dtype=float),
        path,
    )


def load_fem_ring5_history(velocity):
    path = fem_local_stress_history_path(velocity)
    if not path.exists():
        raise FileNotFoundError(f"Missing FEM local-stress 5hL history: {path}")

    frame = pd.read_excel(path, header=None)
    distance_mm = math.nan
    if frame.shape[0] > 0 and frame.shape[1] > 1:
        maybe_distance = frame.iloc[0, 1]
        if pd.notna(maybe_distance):
            distance_mm = float(maybe_distance)

    data = frame.dropna(how="any")
    if data.empty or data.shape[1] < 2:
        raise ValueError(f"FEM local-stress 5hL history has no data rows: {path}")

    return (
        distance_mm,
        data.iloc[:, 0].to_numpy(dtype=float),
        data.iloc[:, 1].to_numpy(dtype=float),
        path,
    )


def load_fem_ring5_stress(velocity, step, hL):
    distance_mm, steps, stress, path = load_fem_ring5_history(velocity)
    matches = np.where(np.isclose(steps, float(step), rtol=0.0, atol=1e-9))[0]
    if len(matches) == 0:
        raise ValueError(
            f"FEM local-stress 5hL history has no row for step {step}: {path}"
        )

    if not np.isfinite(distance_mm):
        distance_mm = 5 * hL * 1000.0
    return distance_mm, float(stress[matches[0]]), path


def load_fem_local_stress_reference(velocity, step, hL):
    profile_distance_mm = None
    profile_stress = None
    profile_source = None

    # The *_last.xlsx files are full spatial profiles for the final dynamic
    # step.  Intermediate-step FEM local-stress data are only available at 5hL.
    if int(step) == FEM_LOCAL_STRESS_PROFILE_STEP:
        profile_distance_mm, profile_stress, profile_path = load_fem_local_stress_profile(
            velocity
        )
        profile_source = str(profile_path)

    try:
        ring5_distance_mm, ring5_stress, ring5_path = load_fem_ring5_stress(
            velocity,
            step,
            hL,
        )
    except FileNotFoundError:
        if profile_distance_mm is None or profile_stress is None:
            raise
        ring5_distance_mm = 5 * hL * 1000.0
        ring5_stress = interpolate_fem_stress(
            ring5_distance_mm,
            profile_distance_mm,
            profile_stress,
        )
        ring5_path = Path(profile_source)

    return FEMLocalStressReference(
        profile_distance_mm=profile_distance_mm,
        profile_stress=profile_stress,
        profile_source=profile_source,
        ring5_distance_mm=ring5_distance_mm,
        ring5_stress=ring5_stress,
        ring5_source=str(ring5_path),
    )


def interpolate_fem_stress(distance_mm, fem_distance_mm, fem_stress):
    return float(np.interp(distance_mm, fem_distance_mm, fem_stress))


def safe_ratio(numerator, denominator):
    if not np.isfinite(denominator) or abs(denominator) < 1e-30:
        return math.nan
    return numerator / denominator


def write_local_stress_outputs(case_folder, config, step, output_dir):
    distance_hL, theta_deg, stress = local_stress_from_reaction(case_folder, config, step)
    fem_reference = load_fem_local_stress_reference(config.velocity, step, config.hL)

    output_dir.mkdir(parents=True, exist_ok=True)

    cal_path = output_dir / "local_stress_cal.csv"
    with cal_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["distance_hL", "distance_mm", "theta_deg", "stress_cal"])
        for i, dist_idx in enumerate(distance_hL):
            distance_mm = dist_idx * config.hL * 1000.0
            for theta, stress_value in zip(theta_deg, stress[i]):
                writer.writerow([dist_idx, distance_mm, theta, stress_value])

    profile_path = output_dir / "local_stress_normalized_profile.csv"
    with profile_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "distance_hL",
            "distance_mm",
            "stress_FEM",
            "stress_cal_ave",
            "stress_cal_max",
            "stress_cal_min",
            "norm_ave",
            "norm_max",
            "norm_min",
        ])
        max_distance = min(config.lL - 1, len(distance_hL))
        for i in range(max_distance):
            dist_idx = int(distance_hL[i])
            distance_mm = dist_idx * config.hL * 1000.0
            if fem_reference.has_profile:
                fem_value = interpolate_fem_stress(
                    distance_mm,
                    fem_reference.profile_distance_mm,
                    fem_reference.profile_stress,
                )
            else:
                fem_value = math.nan
            ring = stress[i]
            stress_ave = float(np.nanmean(ring))
            stress_max = float(np.nanmax(ring))
            stress_min = float(np.nanmin(ring))
            writer.writerow([
                dist_idx,
                distance_mm,
                fem_value,
                stress_ave,
                stress_max,
                stress_min,
                safe_ratio(stress_ave, fem_value),
                safe_ratio(stress_max, fem_value),
                safe_ratio(stress_min, fem_value),
            ])

    ring5_path = output_dir / "local_stress_ring_5hL.csv"
    with ring5_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["theta_deg", "stress_cal", "stress_FEM_5hL", "normalized"])
        if len(distance_hL) >= 5:
            for theta, stress_value in zip(theta_deg, stress[4]):
                writer.writerow([
                    theta,
                    stress_value,
                    fem_reference.ring5_stress,
                    safe_ratio(stress_value, fem_reference.ring5_stress),
                ])

    return {
        "local_stress_cal": str(cal_path),
        "local_stress_profile": str(profile_path),
        "local_stress_ring_5hL": str(ring5_path),
        "fem_local_stress_profile_source": fem_reference.profile_source,
        "fem_local_stress_profile_available": fem_reference.has_profile,
        "fem_local_stress_ring5_source": fem_reference.ring5_source,
        "fem_local_stress_ring5_distance_mm": fem_reference.ring5_distance_mm,
    }


def run_j_integral_outputs(case_folder, config, step, output_dir, j_params=None):
    j_params = resolve_j_integral_params(j_params)
    common = {
        "step_start": step,
        "step_end": step,
        "v": config.velocity,
        "aL": config.aL,
        "lL": config.lL,
        "HL": config.HL,
        "d_theta": config.d_theta,
        "hL": config.hL,
        "c": sp.c,
        "nu": mp.Nu,
        "EE": mp.EE,
        "rho": mp.Rho,
    }
    common.update(j_params.calculator_kwargs())

    hsiga_j = output_dir / "hsiga_J.csv"
    hsiga_calc = JIntegralCalculator(result_root=case_folder, **common)
    hsiga_calc.run(output_file=hsiga_j)

    fem_j = output_dir / "fem_reference_J.csv"
    fem_calc = FEMReferenceJIntegralCalculator(**common)
    fem_calc.run(output_file=fem_j)

    return {
        "hsiga_j": str(hsiga_j),
        "hsiga_kid": str(output_dir / "hsiga_J_KId.csv"),
        "fem_j": str(fem_j),
        "fem_kid": str(output_dir / "fem_reference_J_KId.csv"),
        "j_integral": j_params.calculator_kwargs(),
    }


def read_last_kid(path):
    with Path(path).open(newline="") as f:
        rows = list(csv.reader(f))
    if len(rows) < 2:
        raise ValueError(f"KId file has no data row: {path}")
    angles = np.array([float(x) for x in rows[0][1:]], dtype=float)
    values = np.array([float(x) for x in rows[-1][1:]], dtype=float)
    return angles, values


def write_dsif_normalized(paths, output_dir):
    hsiga_angles, hsiga_k = read_last_kid(paths["hsiga_kid"])
    fem_angles, fem_k = read_last_kid(paths["fem_kid"])
    if len(hsiga_angles) != len(fem_angles) or not np.allclose(hsiga_angles, fem_angles):
        fem_k = np.interp(hsiga_angles, fem_angles, fem_k)

    out_path = output_dir / "dsif_normalized.csv"
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["theta_deg", "KId_cal", "KId_FEM", "normalized"])
        for theta, k_cal, k_fem in zip(hsiga_angles, hsiga_k, fem_k):
            writer.writerow([theta, k_cal, k_fem, safe_ratio(k_cal, k_fem)])

    return {"dsif_normalized": str(out_path)}


def postprocess_case(
    case_folder,
    step=None,
    output_dir=None,
    skip_dsif=False,
    j_params=None,
):
    case_folder = Path(case_folder).resolve()
    config = load_case_config(case_folder)
    target_step = config.final_step if step is None else int(step)
    output_dir = case_folder / "postprocess" if output_dir is None else Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "case_folder": str(case_folder),
        "step": target_step,
    }
    outputs.update(write_local_stress_outputs(case_folder, config, target_step, output_dir))

    if not skip_dsif:
        j_outputs = run_j_integral_outputs(
            case_folder,
            config,
            target_step,
            output_dir,
            j_params=j_params,
        )
        outputs.update(j_outputs)
        outputs.update(write_dsif_normalized(j_outputs, output_dir))

    summary_path = output_dir / "postprocess_summary.json"
    summary_path.write_text(json.dumps(outputs, indent=2))
    outputs["summary"] = str(summary_path)
    return outputs


def postprocess_steps(
    case_folder,
    steps,
    output_dir=None,
    skip_dsif=False,
    j_params=None,
):
    case_folder = Path(case_folder).resolve()
    root_output_dir = case_folder / "postprocess" if output_dir is None else Path(output_dir)
    steps = [int(step) for step in steps]
    outputs_by_step = {}
    for step in steps:
        step_output_dir = root_output_dir / f"step{step:05d}"
        outputs_by_step[str(step)] = postprocess_case(
            case_folder,
            step=step,
            output_dir=step_output_dir,
            skip_dsif=skip_dsif,
            j_params=j_params,
        )

    summary_path = root_output_dir / "postprocess_steps_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "case_folder": str(case_folder),
        "steps": steps,
        "outputs": outputs_by_step,
    }
    summary_path.write_text(json.dumps(summary, indent=2))
    summary["summary"] = str(summary_path)
    return summary


def parse_args():
    parser = argparse.ArgumentParser(
        description="Post-process one dynamic parameter-sweep case.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("case_folder", type=Path)
    parser.add_argument(
        "--step",
        type=int,
        default=None,
        help="Step to post-process; defaults to run_config final step",
    )
    parser.add_argument(
        "--steps",
        type=int,
        nargs="+",
        default=None,
        help=(
            "Post-process multiple steps. Outputs are written under "
            "OUTPUT_DIR/stepNNNNN; cannot be combined with --step."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory; defaults to case_folder/postprocess",
    )
    parser.add_argument(
        "--skip-dsif",
        action="store_true",
        help="Only calculate local-stress outputs",
    )
    parser.add_argument(
        "--rj0",
        type=float,
        default=None,
        help="Override const_jintegral.Rj0",
    )
    parser.add_argument(
        "--rj1",
        type=float,
        default=None,
        help="Override const_jintegral.Rj1",
    )
    parser.add_argument(
        "--wj0",
        type=float,
        default=None,
        help="Override const_jintegral.Wj0",
    )
    parser.add_argument(
        "--wj1",
        type=float,
        default=None,
        help="Override const_jintegral.Wj1",
    )
    parser.add_argument(
        "--ngp",
        type=int,
        default=None,
        help="Override const_jintegral.ngp",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.step is not None and args.steps is not None:
        raise ValueError("--step cannot be combined with --steps")

    j_params = {
        "Rj0": args.rj0,
        "Rj1": args.rj1,
        "Wj0": args.wj0,
        "Wj1": args.wj1,
        "ngp": args.ngp,
    }

    if args.steps is not None:
        outputs = postprocess_steps(
            args.case_folder,
            steps=args.steps,
            output_dir=args.output_dir,
            skip_dsif=args.skip_dsif,
            j_params=j_params,
        )
    else:
        outputs = postprocess_case(
            args.case_folder,
            step=args.step,
            output_dir=args.output_dir,
            skip_dsif=args.skip_dsif,
            j_params=j_params,
        )
    print(json.dumps(outputs, indent=2))


if __name__ == "__main__":
    main()
