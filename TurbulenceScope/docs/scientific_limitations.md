# Scientific limitations

## What is supported today
* Descriptive statements about the **supplied global summaries** (time history of energy and Re_λ; radial spectrum).

## What is *not* supported, and why
1. **No spatial result exists yet.** No real JHTDB cutout was obtained: the build sandbox could not reach JHTDB (policy
   403), and cutouts need a personal token. Questions 4 and 5 are *pending*. No spatial figure of turbulence is presented.
2. **The ParaView pipeline has not been executed** by the author of this repository (ParaView was unavailable where it was
   written). Its numerics are covered by NumPy tests; its ParaView calls are untested.
3. **A −5/3 inertial range is not established.** Fitted slopes depend on the range chosen: −1.57 (k = 3–10), −1.59 (3–23),
   −1.73 (5–15), −1.44 (10–30), −1.45 (10–50), −1.53 (20–60), −1.60 (30–57), −1.96 (30–100). In most ranges the slope is
   *shallower* than −5/3 (−1.44 to −1.60), in one it is steeper (−1.73 for k = 5–15), and it steepens beyond k ≈ 60. The compensated
   spectrum k^{5/3}E(k) is *roughly* flat but not constant: over k = 3–60 it varies by a factor of 1.66 (0.288 at k = 3 to 0.477 at k = 45),
   over k = 10–60 by a factor of 1.41, with a broad maximum near k = 45 and a fall-off beyond. Whether that counts as a plateau is a
   judgement, which is why no −5/3 law is asserted. With Re_λ ≈ 405–438 an inertial range of limited extent is expected, and a non-flat compensated spectrum can arise from
   several effects (e.g. finite Reynolds number, forcing at k ≤ 2, a spectral bottleneck); **the supplied files cannot discriminate
   between them**. The defensible statement is: "E(k) is approximately a power law with exponent between about −1.4 and −1.7 over
   k ≈ 3–60, depending on the range, and does not settle −5/3."
4. **The high-k upturn.** E(k) reaches its minimum at k = 447 and rises for the last ~35 shells. The files do not say why; data of this
   kind can show such behaviour from truncation/aliasing treatment, but that is not confirmed here. Do not interpret k ≳ 300 physically.
5. **Duplicates and irregular sampling.** The history file has three exact duplicate records (t = 2.04, 3.06, 5.06) and two 0.01
   steps (2.04→2.05→2.06) in an otherwise 0.02-spaced series; it also ends at 10.04, not at the 10.056 / 10.048 quoted by the
   JHTDB pages (which themselves disagree). The cause is unknown; the data are kept as supplied, duplicates excluded from statistics only.
6. **Re_λ "around 433".** True for t ≲ 2 (mean 433.0), but over the whole file Re_λ ranges 405.8–438.2 (mean 418.0).
7. **Energy and Re_λ are essentially uncorrelated over the full record** (Pearson r = −0.05): energy rises then falls with a maximum at
   t = 4.26, whereas Re_λ peaks at t = 1.70 and is lowest at t = 6.84. This is an observation; Re_λ also depends on dissipation, which
   these files do not contain, so the files cannot explain it.
8. **Unit and normalisation ambiguity.** No units are given. Σₖ E(k) matches the mean energy to 0.1 %, supporting a common
   normalisation, but this is inferred.

## Limitations that will apply once a cutout is analysed
* **Gradient accuracy.** JHTDB reports ~7 % RMS error for database-side gradients; this project's own cutout gradients are
  second-order finite differences, which under-resolve the smallest scales (a smoothed gradient biases Q low in magnitude and
  fragments or merges structures depending on stride). Report the stride and dx with every result.
* **Cutout size.** Structures larger than the cutout are truncated; boundary cells use one-sided stencils. With the default
  token (≤ 4,096 cells ≈ 16³) a cutout is far too small to show vortex tubes, which span many grid points.
* **Threshold dependence.** Q isosurfaces change in volume, topology (merging/splitting) and apparent "thickness" with the
  threshold. Report several thresholds and the reference statistic.
* **High vorticity ≠ vortex.** Shear layers have high vorticity without swirling motion. Q and λ₂ address this partly but
  are Galilean-invariant, *threshold-dependent* indicators, not proofs. They also disagree in places.
* **No hairpin claim.** Hairpin vortices are a feature associated with wall-bounded flows; nothing here provides evidence
  for them in forced isotropic turbulence, and the project makes no such claim.
* **Single or few snapshots.** Statistics from one small cutout at one time are not converged ensemble statistics.
* **Streamlines** are instantaneous-field trajectories, not particle paths in an unsteady flow.
