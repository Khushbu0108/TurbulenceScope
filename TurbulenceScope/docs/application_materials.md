# Application materials (Handshake AI: ParaView Specialist)

**Read this first: what you can truthfully say today.**

| Claim | Status |
|---|---|
| Built a reproducible Python analysis of JHTDB summary files, with 41 passing tests | **True** (executed) |
| Validated vorticity/Q/λ₂ code against closed-form flows | **True** (executed, synthetic fields) |
| Wrote a ParaView Python pipeline | **True**, but it is **unexecuted** until you run it (docs/reproducibility.md §4) |
| Produced ParaView figures / a `.pvsm` / analysed a real 3-D JHTDB cutout | **Not yet.** Do not say this until it is done |
| Prior research exposure | Your résumé lists *Undergraduate Researcher, UR2PhD Program (Feb–May 2026)*: Python/Pandas/NumPy analysis and scientific visualization of turbulent-flow patterns. You told me this involved ParaView through the DaViM Lab; your résumé does not name ParaView or the lab. **Describe it in your own words and only as far as it is true.** |
| "2+ years of ParaView experience" | The posting requires **at least 2 years of work experience using ParaView**. Your résumé shows about 4 months of research (Feb–May 2026). Do not state or imply 2 years. If the application form asks, answer exactly. A project like this is evidence of skill, not a substitute for the stated requirement; the recruiter decides. |

Bracketed text `[…]` must be replaced with facts only you know.

## 30-second explanation

"TurbulenceScope is a reproducible workflow for the Johns Hopkins forced isotropic turbulence dataset. I analysed the supplied
energy-history and spectrum files with tested Python, found data-quality issues such as duplicated records, and showed that the
spectrum does not cleanly establish a −5/3 range. I also wrote vortex-identification code, validated it against analytic flows,
and a ParaView pipeline for velocity slices, streamlines and Q-isosurfaces. [Once run: and I applied it to a real cutout.] I built
it to be checkable, including a set of tasks for testing whether AI can do this kind of work."

## Two-minute technical walkthrough

1. **Data and honesty about scope (20 s).** Two JHTDB summary files: a time history of energy and Re_λ, and a radial spectrum. They are global
   summaries, not velocity fields, so they can answer questions about time evolution and scale distribution, not about structures.
2. **What the data showed (35 s).** The history file has 507 rows but 504 unique records: three exact duplicates and a short irregular stretch near t = 2.05; I report them and exclude duplicates only from statistics. Energy rises from 0.682 to
   a maximum of 0.730 at t = 4.26 and falls back to 0.680. Re_λ ranges 405.8–438.2; it averages 433.0 for t ≤ 2.04, matching JHTDB's quoted 433, but 418.0 over the whole file. Summing E(k) gives 0.7060 against the mean energy 0.7053, a 0.1 % consistency check.
3. **The −5/3 question (25 s).** Fitted slopes run from −1.44 to −1.73 depending on the range, and the compensated spectrum is only roughly flat (it varies by a factor of about 1.4–1.7 over k ≈ 3–60). So I report "roughly −1.4 to −1.7 over k ≈ 3–60, consistent with but not establishing −5/3", not a Kolmogorov law.
4. **Vortex tooling (25 s).** Vorticity, Q and λ₂ from the velocity-gradient tensor, tested against Taylor–Green, ABC flow, solid-body rotation and pure strain, with measured second-order convergence. Thresholds are swept as multiples of Q_rms because the isosurface depends strongly on the threshold.
5. **ParaView and reproducibility (15 s).** A pipeline with Slice, Stream Tracer, Glyph, Gradient and Contour at fixed colour ranges and camera, plus a launcher; one-command reproduction of the tested analysis. The real-cutout stage is [pending / done: …].

## Five technical interview questions and accurate answers

**1. Why not just isosurface vorticity magnitude to find vortices?**
Vorticity magnitude is also large in shear layers, where fluid does not swirl. Q = ½(‖Ω‖² − ‖S‖²) is positive where rotation dominates strain; λ₂ < 0 flags a pressure-minimum-like core. Both are Galilean-invariant but threshold-dependent and can disagree, so I sweep thresholds and report the reference statistic.

**2. What does your code do for a non-periodic cutout, and what's the accuracy?**
Second-order central differences inside, one-sided second-order stencils on faces, so the outer cells are less accurate; I verified exactness on linear fields and a measured convergence rate of about 2 on Taylor–Green. JHTDB's own README quotes ~7 % RMS error for its database-side gradients, so I treat vorticity/Q from a cutout as approximate and stride-dependent.

**3. Point data or cell data, and why does it matter in ParaView?**
Filters take an association: `['POINTS','velocity']` vs `['CELLS','velocity']`. Gradient, Contour and Stream Tracer need point data (Contour interpolates along edges; Stream Tracer interpolates the vector at points). If the array is cell data, apply Cell Data to Point Data, accepting a smoothing step. A wrong association is a classic reason for an empty output.

**4. How do you compare time steps fairly?**
Fix the colour-map range (`AutomaticRescaleRangeMode = 'Never'` plus an explicit range from the global min/max), fix the camera, and fix isovalues in absolute terms or as multiples of a *single* reference statistic. Auto-rescaling per step makes different fields look alike.

**5. Your spectrum fit gave slopes between −1.44 and −1.73. What do you conclude?**
That E(k) is approximately a power law over k ≈ 3–60 with an exponent that depends on the fitting range, and that this is roughly consistent with, but does not establish, a −5/3 inertial range: the compensated spectrum varies by a factor of ~1.4–1.7 rather than being constant, Re_λ ≈ 400–440 limits the range, and forcing acts at k ≤ 2. A spectrum also contains no information about flux or structures, so I wouldn't claim a cascade or vortex tubes from it.

## GitHub repository description (≤ 350 characters)

> Reproducible Python + ParaView workflow for JHU Turbulence Database forced isotropic turbulence: tested analysis of energy/Re_λ/spectrum data, vorticity & Q/λ₂ tooling validated on analytic flows, a ParaView pipeline, and a proposed benchmark of AI scientific-visualization tasks.

## Paragraph for the Handshake AI application

"I'm a Computer Science student at the University of Houston with undergraduate research experience analysing and visualising turbulent-flow data [describe, in your own words, what you did and with which tools, and for how long]. To go further, I built TurbulenceScope, an open repository around the Johns Hopkins forced isotropic turbulence dataset: tested Python analysis of its energy, Reynolds-number and spectrum data (including reporting duplicated records and declining to claim a −5/3 range the data do not support), vortex-identification code validated against analytic flows, a ParaView Python pipeline for slices, streamlines and Q-isosurfaces [state accurately whether it has been run and on what data], and five proposed tasks, with rubrics and failure modes, for evaluating whether AI systems can perform scientific-visualization work. I'm careful about separating what an analysis supports from what it doesn't, which I think is the core skill in reviewing AI-generated technical responses."

## Demonstration plan (about 8 minutes)

0. **Before you start:** run the ParaView smoke test once and fix any API errors; have a real cutout loaded if you have one.
1. **Repo and honesty slide (30 s):** README → "Verification status" table. Say what is run and what is not.
2. **Data (60 s):** `python scripts/inspect_dataset.py data/ener_Re_time.txt` → show the duplicate-time line. Open `results/figures/fig_energy_Re_history.png`.
3. **Spectrum (90 s):** open `fig_energy_spectrum_loglog.png`, then `fig_spectrum_compensated.png`; show `spectrum_fit_sensitivity.csv`. Make the "no −5/3 claim" point.
4. **Tests (30 s):** `python -m unittest discover -s tests -v` → 41 OK. Name two physics checks (solid-body rotation Q = Ω²; convergence order).
5. **ParaView (3 min):** File → Open the cutout (or `fields.vtk`); Information panel: bounds, arrays, **point** data, velocity 3 components. Then: colour by `velocity` Magnitude (fixed range) → Slice (z-normal) → Glyph on the slice → Stream Tracer (point cloud) → Gradient (vorticity, Q) → Contour at 1×, 2×, 4× Q_rms. Say why each filter is there, and show how the isosurface changes with threshold. Show `q_crosscheck.json`.
6. **Reproducibility (60 s):** show `python scripts/run_paraview.py …` and the `.pvsm`; mention the state-file path caveat.
7. **AI evaluation (60 s):** open `docs/ai_evaluation_tasks.md`; walk through Task 4's two bugs; say plainly that no model has been run against it.
