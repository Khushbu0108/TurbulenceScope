"""Quantitative analysis of the radial spectrum E(k) and the energy/Re_lambda history."""
from __future__ import annotations

import numpy as np
import pandas as pd

KOLMOGOROV_SLOPE = -5.0 / 3.0


def integrate_spectrum(k: np.ndarray, E: np.ndarray) -> float:
    """Total energy implied by the spectrum, sum_k E(k) * dk with dk = 1 (shell width).

    Definition: u'^2-type kinetic energy per unit mass, e = (1/2)<|u|^2> = integral E(k) dk.
    The JHTDB shells are integer-spaced, so the integral is a plain sum. This is a
    consistency check only: it assumes E(k) is the shell-integrated spectrum in the
    same normalisation as the ``energy`` column; that normalisation is not documented
    in the supplied files and is therefore *tested*, not assumed (see results).
    """
    k = np.asarray(k, float)
    E = np.asarray(E, float)
    dk = np.gradient(k)
    return float(np.sum(E * dk))


def local_loglog_slope(k: np.ndarray, E: np.ndarray, window: int = 21) -> np.ndarray:
    """Local power-law exponent d ln E / d ln k from a sliding least-squares fit.

    ``window`` must be odd. Points near the edges use a truncated (shrinking) window
    that is still centred, so the first/last ``window//2`` values are less reliable.
    """
    if window % 2 == 0 or window < 3:
        raise ValueError("window must be an odd integer >= 3")
    lk, lE = np.log(np.asarray(k, float)), np.log(np.asarray(E, float))
    n, h = len(lk), window // 2
    out = np.full(n, np.nan)
    for i in range(n):
        r = min(h, i, n - 1 - i)
        if r < 2:
            continue
        sl = slice(i - r, i + r + 1)
        out[i] = np.polyfit(lk[sl], lE[sl], 1)[0]
    return out


def fit_power_law(k: np.ndarray, E: np.ndarray, kmin: int, kmax: int) -> dict:
    """Least-squares fit of ln E = ln A + s ln k over kmin <= k <= kmax."""
    m = (k >= kmin) & (k <= kmax)
    if m.sum() < 3:
        raise ValueError("need at least 3 points in the fitting range")
    x, y = np.log(k[m].astype(float)), np.log(E[m].astype(float))
    s, lnA = np.polyfit(x, y, 1)
    resid = y - (s * x + lnA)
    ss_res, ss_tot = float(np.sum(resid**2)), float(np.sum((y - y.mean()) ** 2))
    n = int(m.sum())
    se = float(np.sqrt(ss_res / (n - 2) / np.sum((x - x.mean()) ** 2)))
    return {"kmin": int(kmin), "kmax": int(kmax), "n_points": n, "slope": float(s),
            "slope_stderr": se, "amplitude_A": float(np.exp(lnA)),
            "r_squared": 1 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
            "decades": float(np.log10(kmax / kmin))}


def candidate_scaling_ranges(k, E, window=21, tol=0.25, min_width=10) -> pd.DataFrame:
    """Find contiguous k-ranges where the local slope is within ``tol`` of -5/3.

    This is a *screening* tool. A returned range is a candidate, not evidence of a
    Kolmogorov inertial range; see docs/scientific_limitations.md.
    """
    s = local_loglog_slope(k, E, window)
    ok = np.isfinite(s) & (np.abs(s - KOLMOGOROV_SLOPE) <= tol)
    ranges, start = [], None
    for i, flag in enumerate(ok):
        if flag and start is None:
            start = i
        if (not flag or i == len(ok) - 1) and start is not None:
            end = i if flag else i - 1
            if end - start + 1 >= min_width:
                ranges.append((int(k[start]), int(k[end])))
            start = None
    rows = []
    for a, b in ranges:
        try:
            r = fit_power_law(np.asarray(k), np.asarray(E), a, b)
        except ValueError:
            continue
        rows.append(r)
    return pd.DataFrame(rows)


def compensated(k, E, slope=KOLMOGOROV_SLOPE):
    """Compensated spectrum k^{-slope} E(k); a plateau indicates that power law."""
    return np.asarray(k, float) ** (-slope) * np.asarray(E, float)


def spectral_upturn_onset(k, E):
    """Wavenumber of the minimum of E(k) (start of the high-k upturn, if any)."""
    i = int(np.argmin(E))
    return int(np.asarray(k)[i]), bool(i < len(E) - 1)


def history_summary(df: pd.DataFrame, dedupe: bool = True) -> dict:
    """Summary statistics of energy and Re_lambda.

    With ``dedupe=True`` exact duplicate (time, energy, Re) rows are excluded from the
    statistics only (they are still reported by ``qc_energy_history``) so that repeated
    records do not bias means. Time-averages use trapezoidal integration on actual times,
    which is correct for the irregular spacing present in the file.
    """
    d = df.drop_duplicates().reset_index(drop=True) if dedupe else df
    t = d["time"].to_numpy()
    out = {"n_records_used": int(len(d)), "t_start": float(t[0]), "t_end": float(t[-1])}
    for col in ("energy", "Re_lambda"):
        v = d[col].to_numpy()
        T = t[-1] - t[0]
        i_max, i_min = int(np.argmax(v)), int(np.argmin(v))
        out[col] = {
            "initial": float(v[0]), "final": float(v[-1]),
            "min": float(v[i_min]), "t_at_min": float(t[i_min]),
            "max": float(v[i_max]), "t_at_max": float(t[i_max]),
            "mean_time_weighted": float(np.trapezoid(v, t) / T),
            "mean_unweighted": float(v.mean()),
            "std_unweighted": float(v.std(ddof=1)),
            "peak_to_peak": float(v.max() - v.min()),
            "relative_peak_to_peak_pct": float(100 * (v.max() - v.min()) / v.mean()),
            "net_change_pct": float(100 * (v[-1] - v[0]) / v[0]),
        }
    out["pearson_corr_energy_Re"] = float(np.corrcoef(d["energy"], d["Re_lambda"])[0, 1])
    return out


def energy_rate_of_change(df: pd.DataFrame) -> pd.DataFrame:
    """Finite-difference dE/dt on the de-duplicated record using actual time spacing."""
    d = df.drop_duplicates().reset_index(drop=True)
    return pd.DataFrame({"time": d["time"], "dE_dt": np.gradient(d["energy"], d["time"])})


def sign_change_times(t, v):
    """Times (linear interpolation) where a derivative series changes sign."""
    t, v = np.asarray(t, float), np.asarray(v, float)
    idx = np.where(np.sign(v[:-1]) * np.sign(v[1:]) < 0)[0]
    return [float(t[i] - v[i] * (t[i + 1] - t[i]) / (v[i + 1] - v[i])) for i in idx]
