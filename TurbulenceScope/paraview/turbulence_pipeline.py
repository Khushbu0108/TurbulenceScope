"""TurbulenceScope ParaView pipeline (pvpython / pvbatch / ParaView Python Shell).

STATUS: written against the ParaView 6.x ``paraview.simple`` API but NOT EXECUTED by the
author of this repository at the time of writing (ParaView was unavailable in the build
sandbox). Treat every ParaView-specific call as untested until you have run it; the
results table in docs/reproducibility.md records what has actually been verified.

Run (from the repository root):
    pvpython paraview/turbulence_pipeline.py --input results/synthetic_validation/fields.vtk \
             --velocity-array velocity --out-dir results/paraview_synthetic --tag "SYNTHETIC Taylor-Green"
or, inside the ParaView GUI: Tools > Python Shell > Run Script, with sys.argv set via --args in the shell,
or use scripts/run_paraview.py which finds pvpython for you.

Stages
  A  inspect: dimensions, bounds, arrays, centering, time steps  -> inspection.json
  B  velocity: speed colour map, mid-plane slice, stream tracer (+tubes), arrow glyphs on the slice
  C  vortices: Gradient filter (vorticity + Q), Q iso-surfaces at several multiples of Q_rms
  D  time: if the reader reports >1 time step, one image per step with a FIXED colour range/camera
  E  export: PNGs and a .pvsm state file

Why each filter is used (and not more):
  * Slice          - interior cross-section; surface rendering alone hides the volume.
  * Stream Tracer  - trajectories of the instantaneous field; shows topology, not material paths.
  * Glyph (arrows) - direction on a plane where colour already encodes magnitude.
  * Gradient       - computes velocity-gradient derived quantities (vorticity, Q) with VTK's operator.
  * Contour        - iso-surface of Q; the threshold is the scientifically sensitive choice, so it is
                     swept rather than fixed.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

from paraview import servermanager as sm
from paraview.simple import *  # noqa: F401,F403  (standard pvpython idiom)

paraview.simple._DisableFirstRenderCameraReset()


# ----------------------------------------------------------------------------- helpers
def parse_args(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help=".vtk/.vti/.vtk-series, or .xmf/.xdmf (JHTDB cutout)")
    ap.add_argument("--velocity-array", default="velocity", help="POINT-data vector array name (check inspection.json)")
    ap.add_argument("--out-dir", default="results/paraview")
    ap.add_argument("--tag", default="", help="annotation text; use it to label synthetic/test data")
    ap.add_argument("--q-factors", type=float, nargs="+", default=[1.0, 2.0, 4.0, 8.0],
                    help="iso-surface thresholds as multiples of Q_rms")
    ap.add_argument("--size", type=int, nargs=2, default=[1600, 1200])
    ap.add_argument("--no-state", action="store_true")
    return ap.parse_args(argv)


def open_reader(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xmf", ".xdmf"):
        return XDMFReader(registrationName="cutout", FileNames=[path])
    return OpenDataFile(path)


def array_ranges(src, assoc="POINTS"):
    coll = src.PointData if assoc == "POINTS" else src.CellData
    return {k: {"components": coll[k].GetNumberOfComponents(), "range": list(coll[k].GetRange(-1))}
            for k in coll.keys()}


def fixed_camera(view, bounds, direction=(1, 0.8, 0.6)):
    cx, cy, cz = [(bounds[2 * i] + bounds[2 * i + 1]) / 2 for i in range(3)]
    diag = ((bounds[1] - bounds[0]) ** 2 + (bounds[3] - bounds[2]) ** 2 + (bounds[5] - bounds[4]) ** 2) ** 0.5
    d = (direction[0] ** 2 + direction[1] ** 2 + direction[2] ** 2) ** 0.5
    view.CameraFocalPoint = [cx, cy, cz]
    view.CameraPosition = [cx + 1.6 * diag * direction[0] / d, cy + 1.6 * diag * direction[1] / d,
                           cz + 1.6 * diag * direction[2] / d]
    view.CameraViewUp = [0, 0, 1]
    view.CameraParallelScale = 0.6 * diag


def new_view(size, tag, background=(1, 1, 1)):
    v = CreateView("RenderView")
    v.ViewSize = size
    v.OrientationAxesVisibility = 1
    v.Background = list(background)
    try:  # ParaView >= 5.10 palette handling; harmless if unavailable
        v.UseColorPaletteForBackground = 0
    except Exception:
        pass
    if tag:
        t = Text(registrationName="tag")
        t.Text = tag
        td = Show(t, v)
        td.WindowLocation = "Upper Left Corner"
        td.FontSize = 14
        td.Color = [0.7, 0.1, 0.1]
    return v


def shot(view, path, size):
    Render(view)
    SaveScreenshot(path, view, ImageResolution=size)
    print("saved", path)


# ------------------------------------------------------------------------------- main
def main(argv):
    a = parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    out = lambda n: os.path.join(a.out_dir, n)  # noqa: E731

    rdr = open_reader(os.path.abspath(a.input))
    rdr.UpdatePipeline()
    info = rdr.GetDataInformation()

    # ---- A. inspection ---------------------------------------------------------------
    point_arrays, cell_arrays = array_ranges(rdr, "POINTS"), array_ranges(rdr, "CELLS")
    if a.velocity_array not in point_arrays:
        raise SystemExit(f"point array '{a.velocity_array}' not found. Point arrays: {list(point_arrays)}; "
                         f"cell arrays: {list(cell_arrays)}. If velocity is CELL data, apply CellDatatoPointData "
                         f"first or pass the right name.")
    if point_arrays[a.velocity_array]["components"] != 3:
        raise SystemExit(f"'{a.velocity_array}' has {point_arrays[a.velocity_array]['components']} components, expected 3")
    bounds = list(info.GetBounds())
    insp = {"input": a.input, "data_class": rdr.GetXMLName(), "n_points": info.GetNumberOfPoints(),
            "n_cells": info.GetNumberOfCells(), "bounds": bounds, "extent": list(info.GetExtent()),
            "point_arrays": point_arrays, "cell_arrays": cell_arrays,
            "time_steps": list(rdr.TimestepValues) if hasattr(rdr, "TimestepValues") and rdr.TimestepValues else [],
            "velocity_array_centering": "point"}
    with open(out("inspection.json"), "w") as f:
        json.dump(insp, f, indent=2)
    print(json.dumps(insp, indent=2))

    vname = a.velocity_array
    vmin, vmax = point_arrays[vname]["range"][0], point_arrays[vname]["range"][1]
    # 'range' of the vector array with component -1 is the magnitude range.
    center = [(bounds[2 * i] + bounds[2 * i + 1]) / 2 for i in range(3)]
    diag = ((bounds[1] - bounds[0]) ** 2 + (bounds[3] - bounds[2]) ** 2 + (bounds[5] - bounds[4]) ** 2) ** 0.5

    # ---- B. velocity -----------------------------------------------------------------
    view = new_view(a.size, a.tag)
    fixed_camera(view, bounds)

    d_vol = Show(rdr, view)
    d_vol.SetRepresentationType("Outline")
    ColorBy(d_vol, None)
    lut = GetColorTransferFunction(vname)
    lut.ApplyPreset("Viridis (matplotlib)", True)
    lut.AutomaticRescaleRangeMode = "Never"      # fixed range => comparable across time steps/figures
    lut.RescaleTransferFunction(vmin, vmax)

    sl = Slice(registrationName="midplane_z", Input=rdr)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = center
    sl.SliceType.Normal = [0, 0, 1]
    d_sl = Show(sl, view)
    ColorBy(d_sl, ("POINTS", vname, "Magnitude"))
    d_sl.SetScalarBarVisibility(view, True)
    lut.RescaleTransferFunction(vmin, vmax)
    bar = GetScalarBar(lut, view)
    bar.Title = f"|{vname}|"
    bar.ComponentTitle = ""
    shot(view, out("fig01_velocity_magnitude_slice.png"), a.size)   # velocity-magnitude figure

    # glyphs on the slice, coloured by magnitude
    gl = Glyph(registrationName="arrows", Input=sl, GlyphType="Arrow")
    gl.OrientationArray = ["POINTS", vname]
    gl.ScaleArray = ["POINTS", "No scale array"]
    gl.ScaleFactor = 0.04 * diag
    gl.GlyphMode = "Every Nth Point"
    n_slice = sl.GetDataInformation().GetNumberOfPoints()
    gl.Stride = max(1, int(n_slice / 600))
    d_gl = Show(gl, view)
    ColorBy(d_gl, ("POINTS", vname, "Magnitude"))
    shot(view, out("fig02_velocity_slice_with_glyphs.png"), a.size)
    Hide(gl, view)

    # stream tracer from a point cloud near the centre
    st = StreamTracer(registrationName="streamlines", Input=rdr, SeedType="Point Cloud")
    st.Vectors = ["POINTS", vname]
    st.MaximumStreamlineLength = diag
    st.SeedType.Center = center
    st.SeedType.Radius = 0.25 * min(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
    st.SeedType.NumberOfPoints = 150
    tube = Tube(registrationName="streamtubes", Input=st)
    tube.Radius = 0.003 * diag
    d_t = Show(tube, view)
    ColorBy(d_t, ("POINTS", vname, "Magnitude"))
    Hide(sl, view)
    shot(view, out("fig03_streamlines.png"), a.size)
    Hide(tube, view)
    Show(sl, view)

    # ---- C. vorticity and Q ----------------------------------------------------------
    grad = Gradient(registrationName="velocity_gradient", Input=rdr)
    grad.ScalarArray = ["POINTS", vname]
    grad.ComputeVorticity = 1
    grad.ComputeQCriterion = 1
    grad.VorticityArrayName = "Vorticity_pv"
    grad.QCriterionArrayName = "Q_pv"
    grad.UpdatePipeline()

    # Q_rms and an independent cross-check against the NumPy-computed Q, when present.
    q_rms, xcheck = None, {}
    try:
        from vtkmodules.numpy_interface import dataset_adapter as dsa
        import numpy as np
        data = dsa.WrapDataObject(sm.Fetch(grad))
        Qpv = np.asarray(data.PointData["Q_pv"])
        q_rms = float(np.sqrt(np.mean(Qpv ** 2)))
        xcheck = {"Q_pv_rms": q_rms, "Q_pv_min": float(Qpv.min()), "Q_pv_max": float(Qpv.max())}
        if "Q_numpy" in data.PointData.keys():
            Qnp = np.asarray(data.PointData["Q_numpy"])
            interior = np.ones_like(Qnp, dtype=bool)
            # ParaView and NumPy differ at boundaries (stencil choice) -> compare on the bulk
            ext = rdr.GetDataInformation().GetExtent()
            nx, ny, nz = ext[1] - ext[0] + 1, ext[3] - ext[2] + 1, ext[5] - ext[4] + 1
            m = np.zeros((nz, ny, nx), bool)
            m[2:-2, 2:-2, 2:-2] = True
            interior = m.ravel()
            diff = np.abs(Qpv[interior] - Qnp[interior])
            xcheck.update({"max_abs_diff_interior": float(diff.max()),
                           "rel_to_Q_rms": float(diff.max() / q_rms) if q_rms else None,
                           "note": "NumPy path uses 2nd-order central differences; agreement to O(dx^2) is expected"})
    except Exception as e:  # pragma: no cover
        xcheck = {"skipped": f"{type(e).__name__}: {e}"}
    with open(out("q_crosscheck.json"), "w") as f:
        json.dump(xcheck, f, indent=2)
    print("Q cross-check:", xcheck)

    if q_rms:
        Hide(sl, view)
        Show(rdr, view).SetRepresentationType("Outline")
        for fac in a.q_factors:
            con = Contour(registrationName=f"Q_{fac:g}rms", Input=grad)
            con.ContourBy = ["POINTS", "Q_pv"]
            con.Isosurfaces = [fac * q_rms]
            con.ComputeNormals = 1
            con.PointMergeMethod = "Uniform Binning"
            d_c = Show(con, view)
            ColorBy(d_c, ("POINTS", vname, "Magnitude"))
            d_c.SetScalarBarVisibility(view, True)
            lut.RescaleTransferFunction(vmin, vmax)
            txt = Text(registrationName=f"cap{fac:g}")
            txt.Text = f"Q = {fac:g} x Q_rms = {fac * q_rms:.4g}"
            d_txt = Show(txt, view)
            d_txt.WindowLocation = "Upper Right Corner"
            d_txt.Color = [0, 0, 0]
            name = f"fig04_Q_isosurface_{fac:g}rms.png"
            shot(view, out(name), a.size)
            Hide(con, view)
            Hide(txt, view)
        # vorticity magnitude on the mid-plane
        Show(sl, view)
        sl2 = Slice(registrationName="midplane_vort", Input=grad)
        sl2.SliceType = "Plane"
        sl2.SliceType.Origin = center
        sl2.SliceType.Normal = [0, 0, 1]
        Hide(sl, view)
        d2 = Show(sl2, view)
        ColorBy(d2, ("POINTS", "Vorticity_pv", "Magnitude"))
        d2.RescaleTransferFunctionToDataRange(True, False)
        d2.SetScalarBarVisibility(view, True)
        shot(view, out("fig05_vorticity_magnitude_slice.png"), a.size)

    # ---- D. time ---------------------------------------------------------------------
    steps = insp["time_steps"]
    if len(steps) > 1:
        scene = GetAnimationScene()
        scene.UpdateAnimationUsingDataTimeSteps()
        Hide(sl2, view) if q_rms else None
        Show(sl, view)
        for i, t in enumerate(steps):          # LUT range and camera stay fixed (set above)
            scene.AnimationTime = t
            shot(view, out(f"fig06_time_{i:03d}.png"), a.size)
    else:
        print("single time step: temporal visualisation skipped (documented extension, nothing fabricated)")

    # ---- E. state --------------------------------------------------------------------
    if not a.no_state:
        SaveState(out("turbulence_pipeline.pvsm"))
        print("saved", out("turbulence_pipeline.pvsm"))


if __name__ == "__main__":
    main(sys.argv[1:])
