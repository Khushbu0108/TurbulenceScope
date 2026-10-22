#!/usr/bin/env python3
"""Locate pvpython/pvbatch and run paraview/turbulence_pipeline.py with the given arguments.

    python scripts/run_paraview.py --input results/synthetic_validation/fields.vtk \
        --out-dir results/paraview_synthetic --tag "SYNTHETIC Taylor-Green (pipeline smoke test)"

Search order: $PVPYTHON, pvbatch/pvpython on PATH, then common macOS/Linux/Windows install
locations (which can differ between versions: set PVPYTHON if yours is elsewhere).
pvbatch is preferred for headless/offscreen runs; pvpython is used if pvbatch is absent.
"""
from __future__ import annotations

import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PIPELINE = REPO / "paraview" / "turbulence_pipeline.py"


def find_paraview() -> str | None:
    env = os.environ.get("PVPYTHON")
    if env and Path(env).exists():
        return env
    for name in ("pvbatch", "pvpython"):
        p = shutil.which(name)
        if p:
            return p
    patterns = [
        "/Applications/ParaView*.app/Contents/bin/pvbatch",
        "/Applications/ParaView*.app/Contents/bin/pvpython",
        "/opt/ParaView*/bin/pvbatch", "/opt/ParaView*/bin/pvpython",
        "C:/Program Files/ParaView*/bin/pvbatch.exe", "C:/Program Files/ParaView*/bin/pvpython.exe",
    ]
    for pat in patterns:
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[-1]
    return None


def main():
    exe = find_paraview()
    if not exe:
        sys.exit("ParaView's pvbatch/pvpython was not found. Install ParaView or set PVPYTHON=/path/to/pvpython.\n"
                 "In the GUI you can instead use Tools > Python Shell > Run Script on "
                 f"{PIPELINE} (pass options by editing sys.argv in the shell first).")
    cmd = [exe, str(PIPELINE), *sys.argv[1:]]
    print("running:", " ".join(cmd))
    sys.exit(subprocess.call(cmd, cwd=REPO))


if __name__ == "__main__":
    main()
