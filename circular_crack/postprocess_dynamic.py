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

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from const import const_local_mesh as clm
from const import material_property as mp
from const import simulation_params as sp
from jintegral import FEMReferenceJIntegralCalculator, JIntegralCalculator


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


def load_fem_local_stress(velocity):
    v_label = int(round(float(velocity)))
    path = SCRIPT_DIR / "data" / "FEM_local_stress" / f"FEM_v_{v_label}_last.xlsx"
    if not path.exists():
        raise FileNotFoundError(f"Missing FEM local stress reference: {path}")

    frame = pd.read_excel(path, header=None)
    frame = frame.dropna(how="any")
    return frame.iloc[:, 0].to_numpy(dtype=float), frame.iloc[:, 1].to_numpy(dtype=float)


def interpolate_fem_stress(distance_mm, fem_distance_mm, fem_stress):
    return float(np.interp(distance_mm, fem_distance_mm, fem_stress))


def safe_ratio(numerator, denominator):
    if not np.isfinite(denominator) or abs(denominator) < 1e-30:
        return math.nan
    return numerator / denominator


def write_local_stress_outputs(case_folder, config, step, output_dir):
    distance_hL, theta_deg, stress = local_stress_from_reaction(case_folder, config, step)
    fem_distance_mm, fem_stress = load_fem_local_stress(config.velocity)

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
            fem_value = interpolate_fem_stress(distance_mm, fem_distance_mm, fem_stress)
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
            fem_value = interpolate_fem_stress(5 * config.hL * 1000.0, fem_distance_mm, fem_stress)
            for theta, stress_value in zip(theta_deg, stress[4]):
                writer.writerow([theta, stress_value, fem_value, safe_ratio(stress_value, fem_value)])

    return {
        "local_stress_cal": str(cal_path),
        "local_stress_profile": str(profile_path),
        "local_stress_ring_5hL": str(ring5_path),
    }


def run_j_integral_outputs(case_folder, config, step, output_dir):
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


def postprocess_case(case_folder, step=None, output_dir=None, skip_dsif=False):
    case_folder = Path(case_folder).resolve()
    config = load_case_config(case_folder)
    final_step = config.final_step if step is None else int(step)
    output_dir = case_folder / "postprocess" if output_dir is None else Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "case_folder": str(case_folder),
        "step": final_step,
    }
    outputs.update(write_local_stress_outputs(case_folder, config, final_step, output_dir))

    if not skip_dsif:
        j_outputs = run_j_integral_outputs(case_folder, config, final_step, output_dir)
        outputs.update(j_outputs)
        outputs.update(write_dsif_normalized(j_outputs, output_dir))

    summary_path = output_dir / "postprocess_summary.json"
    summary_path.write_text(json.dumps(outputs, indent=2))
    outputs["summary"] = str(summary_path)
    return outputs


def parse_args():
    parser = argparse.ArgumentParser(
        description="Post-process one dynamic parameter-sweep case.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("case_folder", type=Path)
    parser.add_argument("--step", type=int, default=None, help="Step to post-process; defaults to run_config final step")
    parser.add_argument("--output-dir", type=Path, default=None, help="Output directory; defaults to case_folder/postprocess")
    parser.add_argument("--skip-dsif", action="store_true", help="Only calculate local-stress outputs")
    return parser.parse_args()


def main():
    args = parse_args()
    outputs = postprocess_case(
        args.case_folder,
        step=args.step,
        output_dir=args.output_dir,
        skip_dsif=args.skip_dsif,
    )
    print(json.dumps(outputs, indent=2))


if __name__ == "__main__":
    main()
