#!/usr/bin/env python3
"""Inspect a data file and print what is actually in it (no assumptions about field names).

    python scripts/inspect_dataset.py data/spectrum.txt
    python scripts/inspect_dataset.py data/ener_Re_time.txt
    python scripts/inspect_dataset.py path/to/cutout.h5     # HDF5: lists groups/datasets/attrs/shapes/ranges
    python scripts/inspect_dataset.py path/to/cutout.xmf    # XDMF: prints the XML structure (Grid/Topology/Attribute)
    python scripts/inspect_dataset.py path/to/u.npy

For JHTDB cutouts this is the first step: record the dataset names, shapes, axis order and any
coordinate/time attributes it reports, and write them to docs/methodology.md before analysing.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import numpy as np

from turbulencescope import summary_io as sio


def inspect_txt(p: Path):
    first = p.read_text().splitlines()[0].split()
    if "time" in first:
        df = sio.read_energy_history(p)
        qc = sio.qc_energy_history(df)
        print(f"{p.name}: time history, {qc.n_rows} rows, columns {list(df.columns)}")
        print(f"  time {df.time.min()} .. {df.time.max()}; duplicate time values: "
              f"{sorted(qc.duplicate_rows.time.unique().tolist())}; non-finite: {qc.n_nonfinite}")
        print(f"  positive time-step counts: {qc.dt_unique}")
        print(df.describe().to_string())
    else:
        df = sio.read_spectrum(p)
        q = sio.qc_spectrum(df)
        print(f"{p.name}: radial spectrum, {q.n_rows} rows, k = {q.k_min}..{q.k_max}, contiguous={q.k_contiguous}")
        print(f"  E range {df.E.min():.3e} .. {df.E.max():.3e}; non-positive: {q.n_nonpositive}")
        print("  NOTE: a 1-D spectral summary. It contains no spatial or phase information.")


def inspect_h5(p: Path):
    try:
        import h5py
    except ImportError:
        sys.exit("h5py not installed: pip install h5py")
    with h5py.File(p, "r") as f:
        print(f"{p.name}: HDF5; root attrs: {dict(f.attrs)}")

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                info = f"dataset shape={obj.shape} dtype={obj.dtype}"
                if obj.size <= 50_000_000:
                    a = np.asarray(obj)
                    if np.issubdtype(a.dtype, np.number):
                        info += f" min={a.min():.6g} max={a.max():.6g} mean={a.mean():.6g}"
                print(f"  /{name}: {info}; attrs={dict(obj.attrs)}")
            else:
                print(f"  /{name}/: group; attrs={dict(obj.attrs)}")

        f.visititems(visit)


def inspect_xdmf(p: Path):
    import xml.etree.ElementTree as ET
    root = ET.parse(p).getroot()
    for el in root.iter():
        if el.tag in {"Grid", "Topology", "Geometry", "Attribute", "DataItem", "Time"}:
            print(f"  <{el.tag} {dict(el.attrib)}> {(el.text or '').strip()[:100]}")


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    p = Path(sys.argv[1])
    if not p.exists():
        sys.exit(f"not found: {p}")
    ext = p.suffix.lower()
    if ext == ".txt":
        inspect_txt(p)
    elif ext in {".h5", ".hdf5"}:
        inspect_h5(p)
    elif ext in {".xmf", ".xdmf"}:
        inspect_xdmf(p)
    elif ext == ".npy":
        a = np.load(p, mmap_mode="r")
        print(f"{p.name}: npy shape={a.shape} dtype={a.dtype}")
    else:
        sys.exit(f"unsupported extension {ext}")


if __name__ == "__main__":
    main()
