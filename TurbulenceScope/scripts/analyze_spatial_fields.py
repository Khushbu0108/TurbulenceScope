#!/usr/bin/env python3
"""Quantitative analysis of a 3-D velocity cutout: speed, vorticity, Q, lambda_2, threshold sensitivity.

Real data (after you obtain a JHTDB cutout; see docs/reproducibility.md):
    python scripts/analyze_spatial_fields.py --input cutout.h5 --dataset <name> --axis-order zyxc \
        --dx 0.006135923 --out-dir results/spatial

Code-validation mode (SYNTHETIC analytic flow; NOT turbulence, never a research result):
    python scripts/analyze_spatial_fields.py --synthetic taylor_green --n 64 \
        --out-dir results/synthetic_validation

Outputs: <out-dir>/spatial_summary.json, threshold_sensitivity.csv, two PNG figures, and
<out-dir>/fields.vtk (velocity, vorticity, Q_numpy, lambda2_numpy, speed) for ParaView.

--dx is REQUIRED for real data: the grid spacing must come from the file/dataset metadata,
it is not guessed here.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from turbulencescope import spatial as sp

BLUE, ORANGE, GREY, RED = "#1f6fb4", "#d9731a", "#6b6b6b", "#b3261e"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--input", type=Path, help="velocity cutout (.npy or HDF5)")
    src.add_argument("--synthetic", choices=["taylor_green", "abc"], help="analytic validation flow")
    ap.add_argument("--dataset", help="HDF5 dataset name of the velocity array")
    ap.add_argument("--axis-order", default="zyxc", choices=["zyxc", "xyzc"])
    ap.add_argument("--dx", type=float, help="grid spacing (required for --input)")
    ap.add_argument("--n", type=int, default=64, help="grid points per side for --synthetic")
    ap.add_argument("--factors", type=float, nargs="+", default=[0.5, 1.0, 2.0, 4.0, 8.0],
                    help="Q thresholds as multiples of Q_rms")
    ap.add_argument("--out-dir", type=Path, default=REPO / "results" / "spatial")
    a = ap.parse_args()

    a.out_dir.mkdir(parents=True, exist_ok=True)
    if a.synthetic:
        u, dx = (sp.taylor_green if a.synthetic == "taylor_green" else sp.abc_flow)(a.n)
        periodic, tag = True, f"SYNTHETIC TEST FIELD ({a.synthetic}) - code validation only"
        src_desc = {"synthetic": a.synthetic, "n": a.n}
    else:
        if a.dx is None:
            ap.error("--dx is required with --input (take it from the dataset metadata)")
        u, dx = sp.load_velocity(a.input, a.dataset, axis_order=a.axis_order), a.dx
        periodic, tag = False, f"{a.input.name}"
        src_desc = {"input": str(a.input), "dataset": a.dataset, "axis_order": a.axis_order}

    u = u.astype(np.float32 if not a.synthetic else np.float64)
    print(f"velocity shape {u.shape}, dtype {u.dtype}, dx={dx:.9g}, periodic_gradient={periodic}")

    A = sp.velocity_gradient(u, dx, periodic=periodic)
    omega, Q, L2 = sp.vorticity(A), sp.q_criterion(A), sp.lambda2(A)
    spd, wmag, div = sp.speed(u), np.linalg.norm(omega, axis=-1), sp.divergence(A)

    # Edge cells of a non-periodic cutout use one-sided stencils: report interior-only stats too.
    interior = (slice(2, -2),) * 3 if not periodic else (slice(None),) * 3
    rows = sp.q_threshold_sensitivity(Q[interior], factors=a.factors)
    pd.DataFrame(rows).to_csv(a.out_dir / "threshold_sensitivity.csv", index=False)

    summary = {
        "source": src_desc, "shape": list(u.shape[:3]), "dx": dx, "periodic_gradient": periodic,
        "speed": sp.field_statistics(spd), "vorticity_magnitude": sp.field_statistics(wmag),
        "Q": sp.field_statistics(Q), "lambda2": sp.field_statistics(L2),
        "divergence_rms": float(np.sqrt(np.mean(div**2))),
        "divergence_rms_over_vorticity_rms": float(np.sqrt(np.mean(div**2)) / np.sqrt(np.mean(wmag**2))),
        "fraction_Q_positive": float((Q > 0).mean()), "fraction_lambda2_negative": float((L2 < 0).mean()),
        "fraction_Q_pos_and_lambda2_neg": float(((Q > 0) & (L2 < 0)).mean()),
        "q_stats_region": "full volume" if periodic else "interior (2-voxel margin excluded from threshold table only)",
        "threshold_sensitivity": rows,
    }
    (a.out_dir / "spatial_summary.json").write_text(json.dumps(summary, indent=2))

    sp.write_legacy_vtk(a.out_dir / "fields.vtk",
                        {"velocity": u, "vorticity": omega, "Q_numpy": Q, "lambda2_numpy": L2, "speed": spd},
                        spacing=(dx, dx, dx), title=tag)

    # ---- figures
    plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.25, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 200})
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].hist(spd.ravel(), bins=80, color=BLUE)
    ax[0].set(xlabel="|u| (velocity units of input)", ylabel="voxel count", title="Speed distribution")
    ax[1].hist(wmag.ravel(), bins=80, color=ORANGE)
    ax[1].set(xlabel="|ω| (1/time units of input)", ylabel="voxel count", title="Vorticity magnitude")
    fig.suptitle(tag, fontsize=9, color=RED if a.synthetic else "black")
    fig.tight_layout()
    fig.savefig(a.out_dir / "fig_spatial_distributions.png")
    plt.close(fig)

    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    pos = df[df["volume_fraction"] > 0]
    ax[0].semilogy(pos["factor"], pos["volume_fraction"], "o-", color=BLUE)
    for _, r in df[df["volume_fraction"] == 0].iterrows():
        ax[0].annotate("no voxels", (r["factor"], ax[0].get_ylim()[0]), ha="center", va="bottom",
                       fontsize=8, color=GREY)
    ax[0].set(xlabel="threshold / Q_rms", ylabel="volume fraction of {Q > threshold}",
              title="Volume enclosed vs threshold")
    ax[1].plot(df["factor"], df["n_components"], "s-", color=ORANGE)
    ax[1].set(xlabel="threshold / Q_rms", ylabel="connected components (6-connectivity)",
              title="Fragmentation vs threshold")
    fig.suptitle(tag, fontsize=9, color=RED if a.synthetic else "black")
    fig.tight_layout()
    fig.savefig(a.out_dir / "fig_q_threshold_sensitivity.png")
    plt.close(fig)

    print(json.dumps({k: summary[k] for k in ("divergence_rms_over_vorticity_rms", "fraction_Q_positive",
                                              "fraction_Q_pos_and_lambda2_neg")}, indent=2))
    print(pd.DataFrame(rows)[["factor", "threshold", "volume_fraction", "n_components"]].to_string(index=False))
    print("wrote ->", a.out_dir)


if __name__ == "__main__":
    main()
