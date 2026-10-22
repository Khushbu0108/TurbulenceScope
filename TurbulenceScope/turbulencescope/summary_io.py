"""Parsers and validation for the two JHTDB summary text files.

The files are *global summaries*, not spatial fields:

* ``ener_Re_time.txt``: columns ``time energy Re_lambda`` (time history).
* ``spectrum.txt``: columns ``k E(k)`` (radial energy spectrum; the JHTDB page
  describes it as time-averaged).

Parsers never modify the source files and never silently drop rows.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd


class SummaryFormatError(ValueError):
    """Raised when a summary file does not have the expected structure."""


def _data_lines(path: Path) -> tuple[list[str], list[str]]:
    """Split a file into (header_lines, data_lines).

    Header lines are everything up to and including the dashed separator line.
    Blank lines in the data section are ignored.
    """
    text = Path(path).read_text(encoding="ascii", errors="strict")
    lines = text.splitlines()
    sep = next((i for i, ln in enumerate(lines) if ln.strip() and set(ln.strip()) == {"-"}), None)
    if sep is None:
        raise SummaryFormatError(f"{path}: no '-----' separator line found")
    header = lines[: sep + 1]
    data = [ln for ln in lines[sep + 1 :] if ln.strip()]
    return header, data


def _parse(path: Path, columns: list[str]) -> pd.DataFrame:
    header, data = _data_lines(path)
    names = " ".join(header[:-1]).split()
    rows, bad = [], []
    for n, ln in enumerate(data, start=len(header) + 1):
        parts = ln.split()
        if len(parts) != len(columns):
            bad.append((n, ln))
            continue
        try:
            rows.append([float(p) for p in parts])
        except ValueError:
            bad.append((n, ln))
    if bad:
        raise SummaryFormatError(f"{path}: {len(bad)} malformed data line(s), first: {bad[0]}")
    df = pd.DataFrame(rows, columns=columns)
    df.attrs["header_names"] = names
    df.attrs["source"] = str(path)
    return df


def read_energy_history(path: str | Path) -> pd.DataFrame:
    """Read ``ener_Re_time.txt`` into columns ``time, energy, Re_lambda``."""
    return _parse(Path(path), ["time", "energy", "Re_lambda"])


def read_spectrum(path: str | Path) -> pd.DataFrame:
    """Read ``spectrum.txt`` into columns ``k, E`` (``k`` is returned as int)."""
    df = _parse(Path(path), ["k", "E"])
    df["k"] = df["k"].astype(int)
    return df


@dataclass
class HistoryQC:
    n_rows: int
    n_nonfinite: int
    duplicate_rows: pd.DataFrame
    n_duplicate_times_conflicting: int
    is_monotonic_nondecreasing: bool
    dt_unique: dict = field(default_factory=dict)
    irregular_steps: pd.DataFrame = field(default_factory=pd.DataFrame)


def qc_energy_history(df: pd.DataFrame, nominal_dt: float = 0.02, tol: float = 1e-9) -> HistoryQC:
    """Report data-quality issues without altering the data.

    Duplicates are *reported* (with their row indices), never removed here.
    ``n_duplicate_times_conflicting`` counts repeated times whose other columns differ.
    """
    nonfinite = int((~np.isfinite(df[["time", "energy", "Re_lambda"]].to_numpy())).sum())
    dup_mask = df["time"].duplicated(keep=False)
    dups = df[dup_mask].copy()
    dups["source_row_index"] = dups.index
    conflicting = 0
    for _, g in df[dup_mask].groupby("time"):
        if len(g[["energy", "Re_lambda"]].drop_duplicates()) > 1:
            conflicting += 1
    dt = np.diff(df["time"].to_numpy())
    mono = bool(np.all(dt >= -tol))
    pos = dt[dt > tol]
    rounded = pd.Series(np.round(pos, 6)).value_counts().sort_index()
    irregular = []
    for i, d in enumerate(dt):
        if d > tol and abs(d - nominal_dt) > 1e-6:
            irregular.append({"row_before": i, "t_before": df["time"].iloc[i],
                              "t_after": df["time"].iloc[i + 1], "dt": d})
    return HistoryQC(
        n_rows=len(df),
        n_nonfinite=nonfinite,
        duplicate_rows=dups,
        n_duplicate_times_conflicting=conflicting,
        is_monotonic_nondecreasing=mono,
        dt_unique={float(k): int(v) for k, v in rounded.items()},
        irregular_steps=pd.DataFrame(irregular),
    )


@dataclass
class SpectrumQC:
    n_rows: int
    n_nonfinite: int
    n_nonpositive: int
    k_contiguous: bool
    k_min: int
    k_max: int
    duplicated_k: int


def qc_spectrum(df: pd.DataFrame) -> SpectrumQC:
    k = df["k"].to_numpy()
    return SpectrumQC(
        n_rows=len(df),
        n_nonfinite=int((~np.isfinite(df[["k", "E"]].to_numpy(dtype=float))).sum()),
        n_nonpositive=int((df["E"] <= 0).sum()),
        k_contiguous=bool(np.all(np.diff(k) == 1)),
        k_min=int(k.min()),
        k_max=int(k.max()),
        duplicated_k=int(df["k"].duplicated().sum()),
    )
