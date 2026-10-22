#!/usr/bin/env python3
"""Analyse the two supplied JHTDB summary files and write figures + tables.

Usage (from the repository root):
    python scripts/analyze_summaries.py
    python scripts/analyze_summaries.py --data-dir data --out-dir results

Source files are only read, never modified. A SHA-256 of each is recorded in the output JSON.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from turbulencescope import spectrum_analysis as sa
from turbulencescope import summary_io as sio

BLUE, ORANGE, GREY, RED = "#1f6fb4", "#d9731a", "#6b6b6b", "#b3261e"


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def style():
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "figure.dpi": 100, "savefig.dpi": 200})


def fig_history(h, qc, out):
    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for a, col, c, lab in ((ax[0], "energy", BLUE, "Kinetic energy (code units, as supplied)"),
                           (ax[1], "Re_lambda", ORANGE, r"Taylor-scale Reynolds number $Re_\lambda$")):
        a.plot(h["time"], h[col], color=c, lw=1.4)
        a.set_ylabel(lab, fontsize=9)
        dup_t = qc.duplicate_rows["time"].unique()
        a.plot(dup_t, h.set_index("time").loc[dup_t, col].groupby(level=0).first(), "o",
               mfc="none", mec=RED, ms=8, label="duplicated record (kept in file)")
    ax[0].legend(loc="lower right", frameon=False)
    ax[1].set_xlabel("Simulation time (code units)")
    ax[0].set_title("Supplied time history, ener_Re_time.txt (t = 0 to 10.04)")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_spectrum(s, ranges, kup, out):
    k, E = s["k"].to_numpy(), s["E"].to_numpy()
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.loglog(k, E, color=BLUE, lw=1.5, label="E(k), supplied")
    # reference line anchored at the first fitted range's centre, drawn lightly
    if len(ranges):
        r0 = ranges.iloc[0]
        kk = np.array([max(r0.kmin / 1.5, 1), r0.kmax * 1.5])
        kc = np.sqrt(r0.kmin * r0.kmax)
        Ec = r0.amplitude_A * kc ** r0.slope
        ax.loglog(kk, Ec * (kk / kc) ** sa.KOLMOGOROV_SLOPE, "--", color=GREY, lw=1.2,
                  label=r"slope $-5/3$ (reference only)")
        for _, r in ranges.iterrows():
            ax.axvspan(r.kmin, r.kmax, color=ORANGE, alpha=0.12)
    ax.axvline(kup, color=RED, ls=":", lw=1.2, label=f"minimum of E(k) at k = {kup}")
    ax.axvspan(1, 2, color=GREY, alpha=0.15)
    ax.text(1.05, 2e-3, "forced\nshells\n(k ≤ 2)", fontsize=8, color=GREY)
    ax.set_xlabel("Wavenumber k (integer shells)")
    ax.set_ylabel("E(k) (code units, as supplied)")
    ax.set_title("Radial energy spectrum (shaded: candidate near -5/3 ranges)")
    ax.text(0.98, 0.97, "dashed line is a visual reference anchored to the k = 3–23 fit;\nit is not evidence of a -5/3 law",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color=GREY)
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_local_slope(s, out, window=21):
    k, E = s["k"].to_numpy(), s["E"].to_numpy()
    sl = sa.local_loglog_slope(k, E, window)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.semilogx(k, sl, color=BLUE, lw=1.4)
    ax.axhline(sa.KOLMOGOROV_SLOPE, color=GREY, ls="--", label="-5/3")
    ax.axhspan(sa.KOLMOGOROV_SLOPE - 0.25, sa.KOLMOGOROV_SLOPE + 0.25, color=GREY, alpha=0.15,
               label="±0.25 screening band")
    ax.set_ylim(-8, 6)
    ax.set_xlabel("k")
    ax.set_ylabel(f"local slope d ln E / d ln k (window = {window})")
    ax.set_title("Local log-log slope of E(k)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_compensated(s, out):
    k, E = s["k"].to_numpy(), s["E"].to_numpy()
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.loglog(k, sa.compensated(k, E), color=BLUE, lw=1.4)
    ax.set_xlabel("k")
    ax.set_ylabel(r"compensated spectrum $k^{5/3}E(k)$")
    ax.set_title(r"Compensated spectrum: a plateau would indicate $k^{-5/3}$")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, default=REPO / "data")
    ap.add_argument("--out-dir", type=Path, default=REPO / "results")
    a = ap.parse_args()
    figs, tabs = a.out_dir / "figures", a.out_dir / "tables"
    figs.mkdir(parents=True, exist_ok=True)
    tabs.mkdir(parents=True, exist_ok=True)
    style()

    p_hist, p_spec = a.data_dir / "ener_Re_time.txt", a.data_dir / "spectrum.txt"
    h, s = sio.read_energy_history(p_hist), sio.read_spectrum(p_spec)
    hq, sq = sio.qc_energy_history(h), sio.qc_spectrum(s)
    hs = sa.history_summary(h, dedupe=True)

    # --- tables -----------------------------------------------------------
    hq.duplicate_rows.to_csv(tabs / "history_duplicate_records.csv", index=False)
    hq.irregular_steps.to_csv(tabs / "history_irregular_time_steps.csv", index=False)
    rate = sa.energy_rate_of_change(h)
    rate.to_csv(tabs / "energy_rate_of_change.csv", index=False)
    ranges = sa.candidate_scaling_ranges(s["k"].to_numpy(), s["E"].to_numpy())
    ranges.to_csv(tabs / "spectrum_candidate_scaling_ranges.csv", index=False)

    # sensitivity of the fitted exponent to the chosen window
    sens = []
    for kmin, kmax in [(3, 10), (3, 23), (5, 15), (5, 30), (10, 30), (10, 50), (20, 60), (30, 57), (30, 100)]:
        sens.append(sa.fit_power_law(s["k"].to_numpy(), s["E"].to_numpy(), kmin, kmax))
    import pandas as pd
    pd.DataFrame(sens).to_csv(tabs / "spectrum_fit_sensitivity.csv", index=False)

    kup, has_up = sa.spectral_upturn_onset(s["k"], s["E"])
    e_from_spec = sa.integrate_spectrum(s["k"].to_numpy(), s["E"].to_numpy())
    e_mean = hs["energy"]["mean_time_weighted"]
    k_peak_half = int(s.loc[s["E"] >= 0.5 * s["E"].iloc[0], "k"].max())

    early = h.drop_duplicates()
    early = early[early["time"] <= 2.048]
    te = early["time"].to_numpy()
    re_early = float(np.trapezoid(early["Re_lambda"], te) / (te[-1] - te[0]))

    result = {
        "re_lambda_check_vs_jhtdb_readme": {
            "jhtdb_readme_value": 433, "window": "t <= 2.048 (last record at %.2f)" % te[-1],
            "file_time_weighted_mean": re_early,
            "full_period_time_weighted_mean": hs["Re_lambda"]["mean_time_weighted"]},
        "source_files": {p_hist.name: {"sha256": sha256(p_hist), "rows": hq.n_rows},
                         p_spec.name: {"sha256": sha256(p_spec), "rows": sq.n_rows}},
        "history_qc": {
            "rows": hq.n_rows, "non_finite_values": hq.n_nonfinite,
            "time_monotonic_nondecreasing": hq.is_monotonic_nondecreasing,
            "duplicated_time_values": sorted(hq.duplicate_rows["time"].unique().tolist()),
            "duplicate_rows_total": int(len(hq.duplicate_rows)),
            "duplicates_with_conflicting_values": hq.n_duplicate_times_conflicting,
            "positive_dt_counts": hq.dt_unique,
            "irregular_steps": hq.irregular_steps.to_dict("records"),
            "last_time": float(h["time"].iloc[-1]),
        },
        "history_summary_deduplicated": hs,
        "energy_dEdt_sign_change_times": sa.sign_change_times(rate["time"], rate["dE_dt"]),
        "spectrum_qc": sq.__dict__,
        "spectrum_checks": {
            "E1": float(s["E"].iloc[0]), "E_last": float(s["E"].iloc[-1]),
            "E_decay_k1_to_k100_factor": float(s["E"].iloc[0] / s["E"].iloc[99]),
            "k_of_minimum_E": kup, "upturn_present_after_minimum": has_up,
            "sum_E_over_k": e_from_spec,
            "time_weighted_mean_energy": e_mean,
            "relative_difference_pct": 100 * (e_from_spec - e_mean) / e_mean,
            "monotonic_decreasing_up_to_k": int(
                np.argmax(np.diff(s["E"].to_numpy()) > 0) + 1) if (np.diff(s["E"].to_numpy()) > 0).any() else None,
            "k_below_which_E_exceeds_half_of_E1": k_peak_half,
        },
        "candidate_scaling_ranges": ranges.to_dict("records"),
        "fit_sensitivity": sens,
    }
    (tabs / "summary_results.json").write_text(json.dumps(result, indent=2))

    # --- figures ----------------------------------------------------------
    fig_history(h, hq, figs / "fig_energy_Re_history.png")
    fig_spectrum(s, ranges, kup, figs / "fig_energy_spectrum_loglog.png")
    fig_local_slope(s, figs / "fig_spectrum_local_slope.png")
    fig_compensated(s, figs / "fig_spectrum_compensated.png")

    print(json.dumps({k: result[k] for k in ("history_qc", "spectrum_checks", "candidate_scaling_ranges")}, indent=2))
    print("\nenergy:", json.dumps(hs["energy"], indent=2))
    print("Re_lambda:", json.dumps(hs["Re_lambda"], indent=2))
    print("figures ->", figs)
    print("tables  ->", tabs)


if __name__ == "__main__":
    main()
