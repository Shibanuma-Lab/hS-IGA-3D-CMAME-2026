#!/usr/bin/env python3
"""
GPU-assisted L2 norm calculation for verification 5.2 results.

This is intentionally separate from calculate_L2_norm_v2.py. It reuses the
existing data loading and numerical displacement interpolators, while moving
batch quadrature geometry, Sneddon interpolation, and error reduction to CuPy.
"""

import argparse
import csv
import os
import site
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np


def configure_cuda_environment():
    """Point CuPy at CUDA runtime wheels when no system CUDA toolkit exists."""
    if os.environ.get("CUDA_PATH"):
        return

    candidate_roots = []
    try:
        candidate_roots.append(Path(site.getusersitepackages()) / "nvidia" / "cuda_runtime")
    except Exception:
        pass
    for site_path in site.getsitepackages():
        candidate_roots.append(Path(site_path) / "nvidia" / "cuda_runtime")

    for cuda_root in candidate_roots:
        if (cuda_root / "include" / "cuda.h").exists() and (cuda_root / "lib").exists():
            os.environ["CUDA_PATH"] = str(cuda_root)
            break


configure_cuda_environment()

try:
    import cupy as cp
except ImportError:  # pragma: no cover - exercised only when CuPy is missing
    cp = None

from calculate_L2_norm_v2 import GaussQuadrature, L2NormCalculator


def require_cupy():
    """Return the CuPy module or raise a helpful setup error."""
    if cp is None:
        raise RuntimeError(
            "CuPy is not installed in this Python environment. "
            "For an NVIDIA RTX 5070, install a CUDA 12 build, for example: "
            "pipenv run pip install cupy-cuda12x"
        )
    return cp


def check_gpu_runtime():
    """Run a tiny CuPy operation before the expensive L2 setup starts."""
    xp = require_cupy()
    try:
        x = xp.arange(4, dtype=xp.float64)
        float(xp.sum(x * x).get())
    except Exception as exc:
        raise RuntimeError(
            "CuPy is installed, but a test GPU operation failed. "
            "If the message mentions 'libnvrtc.so', install NVRTC with "
            "'python3 -m pip install nvidia-cuda-nvrtc-cu12'. "
            "If it mentions 'Failed to auto-detect CUDA root directory', install "
            "'nvidia-cuda-runtime-cu12' or set CUDA_PATH to the CUDA runtime root. "
            "If it mentions 'cudaErrorInsufficientDriver', update the NVIDIA "
            "Windows/WSL driver or install a CuPy/CUDA runtime version that is "
            "compatible with your current driver."
        ) from exc


class L2NormGPUCalculator(L2NormCalculator):
    """GPU-assisted L2 calculator using the CPU version's mesh readers."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._init_fast_global_grid()
        self._fast_global_points = 0
        self._slow_local_points = 0
        self._fallback_global_points = 0

    def _init_gpu_sneddon(self):
        xp = require_cupy()

        self._sn_r_grid = xp.asarray(self.sneddon.ur1_interp.r_grid, dtype=xp.float64)
        self._sn_z_grid = xp.asarray(self.sneddon.ur1_interp.z_grid, dtype=xp.float64)
        self._sn_ur1 = xp.asarray(self.sneddon.ur1_interp.grid_values, dtype=xp.float64)
        self._sn_ur2 = xp.asarray(self.sneddon.ur2_interp.grid_values, dtype=xp.float64)
        self._sn_uz1 = xp.asarray(self.sneddon.uz1_interp.grid_values, dtype=xp.float64)
        self._sn_uz2 = xp.asarray(self.sneddon.uz2_interp.grid_values, dtype=xp.float64)

    def _init_fast_global_grid(self):
        """Prepare vectorized trilinear interpolation for the structured global mesh."""
        x_coords = np.unique(np.round(self.node_g[:, 0], 12))
        y_coords = np.unique(np.round(self.node_g[:, 1], 12))
        z_coords = np.unique(np.round(self.node_g[:, 2], 12))

        nx, ny, nz = len(x_coords), len(y_coords), len(z_coords)
        if nx * ny * nz != len(self.node_g):
            print("  Fast global interpolation disabled: global nodes are not a full structured grid")
            self._fast_global_enabled = False
            return

        self._fast_global_enabled = True
        self._g_x = x_coords
        self._g_y = y_coords
        self._g_z = z_coords
        self._g_u = self.u_g.reshape(nz, ny, nx, 3)
        print(f"  Fast global interpolation enabled: nx={nx}, ny={ny}, nz={nz}")

    def _bilinear_gpu(self, r, z, values):
        xp = require_cupy()

        r_grid = self._sn_r_grid
        z_grid = self._sn_z_grid
        nr = r_grid.size
        nz = z_grid.size

        inside = (
            (r >= r_grid[0]) & (r <= r_grid[-1]) &
            (z >= z_grid[0]) & (z <= z_grid[-1])
        )

        ir = xp.searchsorted(r_grid, r) - 1
        iz = xp.searchsorted(z_grid, z) - 1
        ir = xp.clip(ir, 0, nr - 2)
        iz = xp.clip(iz, 0, nz - 2)

        r0 = r_grid[ir]
        r1 = r_grid[ir + 1]
        z0 = z_grid[iz]
        z1 = z_grid[iz + 1]

        xi = xp.where(r1 > r0, (r - r0) / (r1 - r0), 0.0)
        eta = xp.where(z1 > z0, (z - z0) / (z1 - z0), 0.0)

        v00 = values[iz, ir]
        v10 = values[iz, ir + 1]
        v11 = values[iz + 1, ir + 1]
        v01 = values[iz + 1, ir]

        out = (
            (1.0 - xi) * (1.0 - eta) * v00 +
            xi * (1.0 - eta) * v10 +
            xi * eta * v11 +
            (1.0 - xi) * eta * v01
        )
        return xp.where(inside, out, 0.0)

    def _sneddon_cartesian_gpu(self, xyz):
        xp = require_cupy()

        x = xyz[:, 0]
        y = xyz[:, 1]
        z = xyz[:, 2]
        r = xp.sqrt(x * x + y * y)
        theta = xp.where(r == 0.0, xp.pi / 2.0, xp.arctan2(y, x))

        p0 = self.sneddon.p0
        c = self.sneddon.c
        ee = self.sneddon.EE
        nu = self.sneddon.nu

        ur0 = -(nu * p0 / ee) * r
        uz0 = (p0 / ee) * z

        ur1 = self._bilinear_gpu(r, z, self._sn_ur1)
        ur2 = self._bilinear_gpu(r, z, self._sn_ur2)
        uz1 = self._bilinear_gpu(r, z, self._sn_uz1)
        uz2 = self._bilinear_gpu(r, z, self._sn_uz2)

        ur_corr = (
            (2.0 * p0 * c * (1.0 + nu) / (xp.pi * ee)) *
            ((1.0 - 2.0 * nu) * ur1 - ur2)
        )
        ur = xp.where(r == 0.0, 0.0, ur_corr)

        uz_corr = -(
            4.0 * p0 * c * (1.0 - nu * nu) / (xp.pi * ee)
        ) * (uz1 + uz2 / (2.0 * (1.0 - nu)))
        uz = xp.where((r >= c) & (z == 0.0), 0.0, uz_corr)

        ux = (ur0 + ur) * xp.cos(theta)
        uy = (ur0 + ur) * xp.sin(theta)
        uz_total = uz0 + uz

        return xp.stack((ux, uy, uz_total), axis=1)

    def _numerical_displacement_cpu_batch(self, xyz_cpu):
        out = np.empty((xyz_cpu.shape[0], 3), dtype=np.float64)
        local_mask = self._is_in_local_region_cpu_batch(xyz_cpu)

        global_mask = ~local_mask
        if np.any(global_mask):
            global_points = xyz_cpu[global_mask]
            if self._fast_global_enabled:
                out[global_mask] = self._trilinear_global_cpu_batch(global_points)
                self._fast_global_points += int(global_points.shape[0])
            else:
                for i, (x, y, z) in zip(np.where(global_mask)[0], global_points):
                    out[i] = np.array([self.ugx(x, y, z), self.ugy(x, y, z), self.ugz(x, y, z)])
                self._fallback_global_points += int(global_points.shape[0])

        if np.any(local_mask):
            local_points = xyz_cpu[local_mask]
            for i, (x, y, z) in zip(np.where(local_mask)[0], local_points):
                out[i] = np.array([self.uglx(x, y, z), self.ugly(x, y, z), self.uglz(x, y, z)])
            self._slow_local_points += int(local_points.shape[0])

        return out

    def _is_in_local_region_cpu_batch(self, xyz):
        x = xyz[:, 0]
        y = xyz[:, 1]
        z = xyz[:, 2]
        z_mask = (0.0 <= z) & (z <= 0.25)

        phi = np.arctan2(y, x)
        r = np.sqrt(x * x + y * y)
        nL_theta = self.nL_theta
        if nL_theta == 0:
            return z_mask & (0.75 <= r) & (r <= 1.25)

        half_pi = 0.5 * np.pi
        delta_theta = half_pi / nL_theta
        phi_mod = np.mod(phi, delta_theta)
        cos_term1 = np.cos(half_pi / (2 * nL_theta))
        cos_term2 = np.cos(phi_mod - half_pi / (2 * nL_theta))

        regular = np.abs(cos_term2) >= 1e-10
        factor = np.empty_like(r)
        factor[regular] = cos_term1 / cos_term2[regular]
        factor[~regular] = 1.0

        r_min = 0.75 * factor
        r_max = 1.25 * factor
        region = (r_min <= r) & (r <= r_max)
        fallback_region = (0.75 <= r) & (r <= 1.25)
        return z_mask & np.where(regular, region, fallback_region)

    def _trilinear_global_cpu_batch(self, xyz):
        x = xyz[:, 0]
        y = xyz[:, 1]
        z = xyz[:, 2]
        xg, yg, zg = self._g_x, self._g_y, self._g_z
        grid = self._g_u

        inside = (
            (x >= xg[0]) & (x <= xg[-1]) &
            (y >= yg[0]) & (y <= yg[-1]) &
            (z >= zg[0]) & (z <= zg[-1])
        )

        out = np.zeros((xyz.shape[0], 3), dtype=np.float64)
        if not np.any(inside):
            return out

        idx = np.where(inside)[0]
        xi = np.searchsorted(xg, x[idx], side="right") - 1
        yi = np.searchsorted(yg, y[idx], side="right") - 1
        zi = np.searchsorted(zg, z[idx], side="right") - 1

        xi = np.clip(xi, 0, len(xg) - 2)
        yi = np.clip(yi, 0, len(yg) - 2)
        zi = np.clip(zi, 0, len(zg) - 2)

        x0, x1 = xg[xi], xg[xi + 1]
        y0, y1 = yg[yi], yg[yi + 1]
        z0, z1 = zg[zi], zg[zi + 1]

        tx = np.where(x1 > x0, (x[idx] - x0) / (x1 - x0), 0.0)[:, None]
        ty = np.where(y1 > y0, (y[idx] - y0) / (y1 - y0), 0.0)[:, None]
        tz = np.where(z1 > z0, (z[idx] - z0) / (z1 - z0), 0.0)[:, None]

        c000 = grid[zi, yi, xi]
        c100 = grid[zi, yi, xi + 1]
        c110 = grid[zi, yi + 1, xi + 1]
        c010 = grid[zi, yi + 1, xi]
        c001 = grid[zi + 1, yi, xi]
        c101 = grid[zi + 1, yi, xi + 1]
        c111 = grid[zi + 1, yi + 1, xi + 1]
        c011 = grid[zi + 1, yi + 1, xi]

        c00 = c000 * (1.0 - tx) + c100 * tx
        c10 = c010 * (1.0 - tx) + c110 * tx
        c01 = c001 * (1.0 - tx) + c101 * tx
        c11 = c011 * (1.0 - tx) + c111 * tx
        c0 = c00 * (1.0 - ty) + c10 * ty
        c1 = c01 * (1.0 - ty) + c11 * ty
        out[idx] = c0 * (1.0 - tz) + c1 * tz
        return out

    def calculate_gpu(self, quadrature_order=8, batch_elements=128):
        xp = require_cupy()
        if batch_elements < 1:
            raise ValueError("batch_elements must be >= 1")

        device = xp.cuda.Device()
        props = xp.cuda.runtime.getDeviceProperties(device.id)
        device_name = props["name"].decode() if isinstance(props["name"], bytes) else props["name"]
        print(f"\nCalculating L2 norm with GPU-assisted quadrature...")
        print(f"  CUDA device: {device.id} ({device_name})")
        print(f"  Quadrature order: {quadrature_order}")
        print(f"  Batch elements: {batch_elements}")

        start_time = time.time()
        self._init_gpu_sneddon()

        quad = GaussQuadrature(quadrature_order)
        q_points = xp.asarray(quad.points_3d, dtype=xp.float64)
        q_weights = xp.asarray(quad.weights_3d, dtype=xp.float64)

        xi = q_points[:, 0]
        eta = q_points[:, 1]
        zeta = q_points[:, 2]

        shape = xp.stack((
            0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta),
            0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta),
            0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta),
            0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta),
            0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + zeta),
            0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + zeta),
            0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + zeta),
            0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + zeta),
        ), axis=1)

        dshape = xp.empty((q_points.shape[0], 3, 8), dtype=xp.float64)
        dshape[:, 0, 0] = -0.125 * (1.0 - eta) * (1.0 - zeta)
        dshape[:, 0, 1] = 0.125 * (1.0 - eta) * (1.0 - zeta)
        dshape[:, 0, 2] = 0.125 * (1.0 + eta) * (1.0 - zeta)
        dshape[:, 0, 3] = -0.125 * (1.0 + eta) * (1.0 - zeta)
        dshape[:, 0, 4] = -0.125 * (1.0 - eta) * (1.0 + zeta)
        dshape[:, 0, 5] = 0.125 * (1.0 - eta) * (1.0 + zeta)
        dshape[:, 0, 6] = 0.125 * (1.0 + eta) * (1.0 + zeta)
        dshape[:, 0, 7] = -0.125 * (1.0 + eta) * (1.0 + zeta)

        dshape[:, 1, 0] = -0.125 * (1.0 - xi) * (1.0 - zeta)
        dshape[:, 1, 1] = -0.125 * (1.0 + xi) * (1.0 - zeta)
        dshape[:, 1, 2] = 0.125 * (1.0 + xi) * (1.0 - zeta)
        dshape[:, 1, 3] = 0.125 * (1.0 - xi) * (1.0 - zeta)
        dshape[:, 1, 4] = -0.125 * (1.0 - xi) * (1.0 + zeta)
        dshape[:, 1, 5] = -0.125 * (1.0 + xi) * (1.0 + zeta)
        dshape[:, 1, 6] = 0.125 * (1.0 + xi) * (1.0 + zeta)
        dshape[:, 1, 7] = 0.125 * (1.0 - xi) * (1.0 + zeta)

        dshape[:, 2, 0] = -0.125 * (1.0 - xi) * (1.0 - eta)
        dshape[:, 2, 1] = -0.125 * (1.0 + xi) * (1.0 - eta)
        dshape[:, 2, 2] = -0.125 * (1.0 + xi) * (1.0 + eta)
        dshape[:, 2, 3] = -0.125 * (1.0 - xi) * (1.0 + eta)
        dshape[:, 2, 4] = 0.125 * (1.0 - xi) * (1.0 - eta)
        dshape[:, 2, 5] = 0.125 * (1.0 + xi) * (1.0 - eta)
        dshape[:, 2, 6] = 0.125 * (1.0 + xi) * (1.0 + eta)
        dshape[:, 2, 7] = 0.125 * (1.0 - xi) * (1.0 + eta)

        integral_error = xp.asarray(0.0, dtype=xp.float64)
        integral_exact = xp.asarray(0.0, dtype=xp.float64)

        elements = self.bg_mesh.elements
        nodes = self.bg_mesh.nodes
        n_elements = len(elements)

        for start in range(0, n_elements, batch_elements):
            end = min(start + batch_elements, n_elements)
            elem_nodes = elements[start:end]
            coords = xp.asarray(nodes[elem_nodes], dtype=xp.float64)

            # Broadcasted reductions avoid CuPy einsum/matmul dispatch, which
            # would require cuBLASLt for these small fixed-size contractions.
            xyz = xp.sum(shape[None, :, :, None] * coords[:, None, :, :], axis=2)
            jac = xp.sum(dshape[None, :, :, :, None] * coords[:, None, None, :, :], axis=3)
            det_j = xp.abs(
                jac[:, :, 0, 0] * (jac[:, :, 1, 1] * jac[:, :, 2, 2] - jac[:, :, 1, 2] * jac[:, :, 2, 1]) -
                jac[:, :, 0, 1] * (jac[:, :, 1, 0] * jac[:, :, 2, 2] - jac[:, :, 1, 2] * jac[:, :, 2, 0]) +
                jac[:, :, 0, 2] * (jac[:, :, 1, 0] * jac[:, :, 2, 1] - jac[:, :, 1, 1] * jac[:, :, 2, 0])
            )

            xyz_flat = xyz.reshape((-1, 3))
            weight_flat = (det_j * q_weights[None, :]).reshape(-1)

            # Numerical interpolation still uses the trusted CPU path.
            xyz_cpu = xp.asnumpy(xyz_flat)
            u_num_cpu = self._numerical_displacement_cpu_batch(xyz_cpu)
            u_num = xp.asarray(u_num_cpu, dtype=xp.float64)

            u_ana = self._sneddon_cartesian_gpu(xyz_flat)
            error = u_num - u_ana

            integral_error += xp.sum(weight_flat * xp.sum(error * error, axis=1))
            integral_exact += xp.sum(weight_flat * xp.sum(u_ana * u_ana, axis=1))

            if end == n_elements or (start // batch_elements) % 10 == 0:
                progress = 100.0 * end / n_elements
                print(f"  Progress: {progress:.1f}% ({end}/{n_elements} elements)")

        integral_error = float(integral_error.get())
        integral_exact = float(integral_exact.get())
        relative_l2 = float(np.sqrt(integral_error / integral_exact))
        elapsed_time = time.time() - start_time
        total_interp_points = self._fast_global_points + self._slow_local_points + self._fallback_global_points

        print(f"\n  Integration complete in {elapsed_time:.1f}s")
        print(f"  Integral error = {integral_error:.6e}")
        print(f"  Integral exact = {integral_exact:.6e}")
        print(f"  Relative L2 norm = {relative_l2:.6e}")
        if total_interp_points:
            print(
                "  Numerical interpolation points: "
                f"fast global={self._fast_global_points} "
                f"({100.0 * self._fast_global_points / total_interp_points:.1f}%), "
                f"slow local={self._slow_local_points} "
                f"({100.0 * self._slow_local_points / total_interp_points:.1f}%), "
                f"fallback global={self._fallback_global_points}"
            )

        return {
            "hL": self.hL,
            "hG": self.config["hG"],
            "rGL": self.config["rGL"],
            "dof": (len(self.node_g) + len(self.node_l)) * 3,
            "sneddon_file": self.sneddon_file,
            "integral_error": integral_error,
            "integral_exact": integral_exact,
            "relative_L2_norm": relative_l2,
            "computation_time": elapsed_time,
            "gpu_batch_elements": batch_elements,
        }


def process_all_results_gpu(base_dir="results/verification_5_2", rGL=2,
                            sneddon_file="sneddon_python.mat",
                            output_file=None, batch_elements=128,
                            max_folders=None):
    check_gpu_runtime()

    base_path = Path(base_dir)
    rGL_folder = base_path / f"rGL{rGL}_0.25"
    if not rGL_folder.exists():
        print(f"Error: Folder not found: {rGL_folder}")
        return []

    def extract_hL(folder):
        try:
            return float(folder.name.split("_")[1])
        except Exception:
            return 0.0

    result_folders = [
        d for d in rGL_folder.iterdir()
        if d.is_dir() and d.name.startswith("hL_")
    ]
    result_folders = sorted(result_folders, key=extract_hL, reverse=True)
    if max_folders is not None:
        result_folders = result_folders[:max_folders]

    if output_file is None:
        output_file = base_path / f"L2_norm_gpu_rGL{rGL}.csv"
    output_file = Path(output_file)

    fieldnames = [
        "rGL", "hL", "hG", "dof", "relative_L2_norm", "sneddon_file",
        "computation_time", "gpu_batch_elements", "timestamp"
    ]
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

    print(f"\nProcessing rGL = {rGL} with GPU-assisted L2 calculation")
    print(f"Found {len(result_folders)} result folders")
    print(f"CSV file initialized: {output_file}")

    results = []
    for i, folder in enumerate(result_folders):
        print(f"\n--- Processing {i + 1}/{len(result_folders)}: {folder.name} ---")
        folder_start = time.time()
        try:
            calc = L2NormGPUCalculator(folder, sneddon_file=sneddon_file)
            result = calc.calculate_gpu(quadrature_order=8, batch_elements=batch_elements)
            results.append(result)
            print(f"\nFolder {folder.name} completed in {time.time() - folder_start:.1f}s")

            with open(output_file, "a", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writerow({
                    "rGL": result["rGL"],
                    "hL": result["hL"],
                    "hG": result["hG"],
                    "dof": result["dof"],
                    "relative_L2_norm": result["relative_L2_norm"],
                    "sneddon_file": result["sneddon_file"],
                    "computation_time": result["computation_time"],
                    "gpu_batch_elements": result["gpu_batch_elements"],
                    "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                })
        except Exception as exc:
            print(f"Error processing {folder.name}: {exc}")
            raise

    print(f"\nAll results saved to: {output_file}")
    print(f"Total completed: {len(results)}/{len(result_folders)} folders")
    return results


def main():
    parser = argparse.ArgumentParser(
        description="GPU-assisted L2 norm calculation for S-IGA verification results"
    )
    parser.add_argument("--rGL", type=int, default=2,
                        help="rGL value to process (default: 2)")
    parser.add_argument("--base-dir", type=str, default="results/verification_5_2",
                        help="Base results directory (default: results/verification_5_2)")
    parser.add_argument("--sneddon-file", type=str, default="sneddon_python.mat",
                        help="Sneddon data file (default: sneddon_python.mat)")
    parser.add_argument("--output", type=str,
                        help="Output CSV file (default: auto-generated)")
    parser.add_argument("--batch-elements", type=int, default=128,
                        help="Background elements per GPU batch (default: 128)")
    parser.add_argument("--max-folders", type=int,
                        help="Process only the first N result folders (debug/testing)")
    args = parser.parse_args()

    if args.batch_elements < 1:
        parser.error("--batch-elements must be >= 1")
    if args.max_folders is not None and args.max_folders < 1:
        parser.error("--max-folders must be >= 1")

    process_all_results_gpu(
        base_dir=args.base_dir,
        rGL=args.rGL,
        sneddon_file=args.sneddon_file,
        output_file=args.output,
        batch_elements=args.batch_elements,
        max_folders=args.max_folders,
    )


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
