# Reproducibility

## 1. Verification status (what has actually been run)

Environment where the tested parts were run: Linux, Python 3.13, NumPy 2.5, pandas 3.0, Matplotlib 3.11, SciPy 1.18.

| Component | Status |
|---|---|
| Parsing, QC, summary statistics, spectrum analysis, figures | **Executed**; outputs in `results/` |
| Spatial maths (gradient, vorticity, Q, λ₂, thresholds, VTK writer) | **Executed and tested on analytic flows** |
| Unit tests | **41 tests, all passing** (`python -m unittest discover -s tests -v`); not run under pytest (pytest was not installable in the build sandbox, but the tests are plain `unittest` and pytest-collectable) |
| `scripts/analyze_spatial_fields.py` | **Executed on a synthetic Taylor–Green field only** (`results/synthetic_validation/`) |
| `scripts/inspect_dataset.py` | **Executed on the two `.txt` files**; HDF5/XDMF branches **untested** (no h5py, no sample file) |
| Real JHTDB cutout | **Not obtained.** The build sandbox could not reach JHTDB (policy 403) and a token is required |
| `paraview/turbulence_pipeline.py`, `scripts/run_paraview.py` | **Syntax-checked only. Not executed in ParaView.** No `.pvsm` is committed because none has been generated and opened |
| ParaView-generated figures | **None exist yet** |

The ParaView 6.1.0-RC1 app is installed on the owner's Mac (visible in the screenshots supplied), so the pipeline can be run
there; see §4.

## 2. Reproduce the tested analysis

```bash
git clone <your-repo-url> TurbulenceScope && cd TurbulenceScope
python -m venv .venv && source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python -m unittest discover -s tests -v                  # or: pytest -q
python scripts/inspect_dataset.py data/spectrum.txt
python scripts/inspect_dataset.py data/ener_Re_time.txt
python scripts/analyze_summaries.py                      # writes results/figures/*.png and results/tables/*
python scripts/analyze_spatial_fields.py --synthetic taylor_green --n 64 --out-dir results/synthetic_validation
```

Expected: tests print `OK`; `results/tables/summary_results.json` contains the SHA-256 of the data files
(`edde46fd…` for `ener_Re_time.txt`, `b45c6438…` for `spectrum.txt` for the copies used here).

## 3. Obtain a real JHTDB cutout (you must do this yourself)

I could not do this, and I will not guess endpoints, variable names or credentials. What the official pages state
(read 2026-10-21):

1. Open the dataset page, https://turbulence.pha.jhu.edu/Forced_isotropic_turbulence.aspx, and follow its **authorization token**
   link to request a personal turbulence token. Cutouts with the default token are limited to **4,096 cells**
   (about 16³) — enough for a pipeline smoke test, too small to show vortex tubes. A personal token raises the limit.
   Keep the token private; do not commit it (`.gitignore` excludes `*token*`, `.env`).
2. Use the **New Cutout** form (https://turbulence.idies.jhu.edu/database/cutout) or scripted `getcutout` calls in a
   SciServer notebook (the page recommends scripting for cutouts above about 2 GB). Results are written to your private
   SciServer scratch storage (listed under *Cutout History*), then you download them.
3. Request a **small** region (e.g. 128³ or 256³ at stride 1) at **one** time step first. Do not request the 1024³ dataset.
4. Whether the service offers HDF5 + XDMF for your request, and the dataset/variable names and axis order, are **not stated
   on the page I could read**. Inspect what you receive:
   ```bash
   python scripts/inspect_dataset.py path/to/your_cutout.h5
   ```
   and record the names, shapes, attributes, time stamp and index range in `docs/methodology.md`.
5. Analyse (grid spacing from the metadata; a stride > 1 multiplies dx):
   ```bash
   python scripts/analyze_spatial_fields.py --input path/to/cutout.h5 --dataset <velocity-dataset-name> \
       --axis-order zyxc --dx <spacing> --out-dir results/spatial
   ```
   Check `spatial_summary.json → divergence_rms_over_vorticity_rms`. If it is not small, the axis order or `--dx` is wrong.
6. Cite JHTDB (see `data/README.md`).

If the file is XDMF/HDF5 you can also open the `.xmf` directly in ParaView; the pipeline accepts `.xmf/.xdmf` via
`XDMFReader` (untested), and `--velocity-array` must be the *point* array name shown by `inspect_dataset.py`.
If the velocity is cell-centred, apply *Cell Data to Point Data* first.

## 4. Run the ParaView pipeline (on a machine with ParaView 6.x)

Smoke test on the synthetic field produced in §2:

```bash
python scripts/run_paraview.py --input results/synthetic_validation/fields.vtk --velocity-array velocity \
    --out-dir results/paraview_synthetic --tag "SYNTHETIC Taylor-Green (pipeline smoke test)"
```

`run_paraview.py` looks for `$PVPYTHON`, then `pvbatch`/`pvpython` on `PATH`, then common install folders. On macOS the
binaries are typically inside the app bundle, e.g. `/Applications/ParaView-6.1.0-RC1.app/Contents/bin/pvpython`
(verify the path on your machine; set `PVPYTHON` if different). In the GUI: *Tools → Python Shell → Run Script*.

Expected outputs in `--out-dir`: `inspection.json`, `fig01…fig05*.png`, `q_crosscheck.json`, `turbulence_pipeline.pvsm`.
`q_crosscheck.json` compares ParaView's Q with this repository's NumPy Q on the interior — a large relative difference means a
mismatch of array order or spacing.

**Expect to fix API details.** The script targets the 6.x `paraview.simple` API but was never run. Likely trouble spots:
`Gradient` property names (`ComputeQCriterion`, `QCriterionArrayName`), `StreamTracer` `SeedType` strings, `Glyph.GlyphMode`.
Run it, copy any traceback, and fix or report it; update the table in §1 afterwards.

## 5. Opening a state file

After a successful run, `File → Load State…` and choose `turbulence_pipeline.pvsm`. A state stores the *absolute* input path;
in the dialog choose "Search files under specified directory" and point it at the folder containing the data file.
