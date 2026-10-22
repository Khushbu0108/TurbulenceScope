# Methodology

Every quantity below lists its definition, required inputs, assumptions/units and how it was validated.
"Validated" means checked by an automated test in `tests/` unless stated otherwise.

## 1. What each input can and cannot support

| Input | Is | Supports | Does **not** support |
|---|---|---|---|
| `ener_Re_time.txt` | global scalars vs. time | Q1, Q2 | any spatial statement; why the values change |
| `spectrum.txt` | 1-D shell spectrum E(k) | Q3 | spatial structures, phases, vortices, the time-dependence of the spectrum (JHTDB describes it as time-averaged; the file itself does not say) |
| a 3-D velocity cutout (**not yet obtained**) | u(x,y,z) at one or more times | Q4, Q5 | statements about scales larger than the cutout |

Units: JHTDB describes the data in simulation units (domain [0, 2π]³, ν = 0.000185). The supplied files carry no unit
labels, so figures say "code units, as supplied".

## 2. Summary-file analysis (executed and tested)

**Parsing** (`turbulencescope/summary_io.py`). Header up to the dashed separator is ignored; every data line must have the
expected number of float fields or `SummaryFormatError` is raised (no row is silently dropped). Source files are only read.

**QC.** Non-finite count; monotonic time; duplicates by exact `time` match with a check whether the other columns
conflict; step sizes versus nominal 0.02. Duplicates are *reported* (`results/tables/history_duplicate_records.csv`).
For statistics only, exact duplicate rows are excluded so they do not bias averages; the file itself is untouched.

**Time-weighted mean.** ⟨f⟩ = (1/T)∫ f dt by the trapezoid rule on the actual times (correct for the non-uniform spacing).
Input: `time`, `energy` or `Re_lambda`. Validated against a hand-computed irregular series.

**dE/dt.** `numpy.gradient(energy, time)` on the de-duplicated record. Sign changes located by linear interpolation. This
is a diagnostic of the supplied curve; it says nothing about forcing, dissipation or transfer.

**Spectrum consistency check.** Total energy implied by E(k): Σₖ E(k)Δk with Δk = 1. This presumes E(k) is the
shell-integrated spectrum in the same normalisation as `energy`; that is *tested, not assumed*: the sum is 0.70605 versus
the time-weighted mean energy 0.70531 (difference +0.10 %). Agreement is evidence the two files share a normalisation and
that the spectrum is time-averaged over a comparable period; it is not proof of either.

**Re_λ check.** JHTDB's README quotes Re_λ = 433 averaged over t = 0 … 2.048. The file's time-weighted mean over the
records up to t = 2.04 is 433.02. (Over the whole file the mean is 418.05.)

**Local slope and power-law fits.** Local exponent d ln E/d ln k by sliding least squares in log-log space (odd window,
default 21). `fit_power_law` fits ln E = ln A + s ln k over [k_min, k_max] and reports slope, standard error, R², and the
number of decades. `candidate_scaling_ranges` is a *screening* step (local slope within ±0.25 of −5/3 for ≥ 10 consecutive
wavenumbers). Validated: recovers an exact −5/3 and −2.5 power law to machine precision and finds the plateau in a synthetic spectrum.
Limits: a high R² on a short range is not evidence of a power law; the standard error ignores correlation between
neighbouring shells and model error, so it understates uncertainty. The fit-sensitivity table (`results/tables/spectrum_fit_sensitivity.csv`)
is the honest summary.

## 3. Spatial analysis (code validated on analytic flows; **no real cutout analysed yet**)

Array convention in code: `u[nx, ny, nz, 3]`, components (u_x, u_y, u_z).

* **Velocity-gradient tensor** A_ij = ∂u_i/∂x_j. Non-periodic (cutout): central differences inside, second-order one-sided
  differences on faces (`numpy.gradient(edge_order=2)`). Periodic (full-box synthetic fields only): central with wrap-around.
  Cross-check: a Fourier-spectral gradient. *Validated:* exact on linear fields including boundaries; second-order
  convergence measured on Taylor–Green (observed rates 1.8–2.2); spectral gradient matches closed form to 1e-10.
* **Vorticity** ω = ∇×u. Validated against closed-form Taylor–Green vorticity, ABC flow (Beltrami: ω = u), and solid-body
  rotation (ω_z = 2Ω).
* **Strain S = (A+Aᵀ)/2, rotation Ω = (A−Aᵀ)/2.**
* **Q-criterion** Q = ½(‖Ω‖²_F − ‖S‖²_F); regions Q > 0 are rotation-dominated. Validated: solid-body rotation gives Q = Ω²;
  pure strain gives Q = −a²; for divergence-free fields Q = −½ A_ij A_ji; Galilean invariance.
* **λ₂** = middle eigenvalue of S² + Ω²; vortex core where λ₂ < 0. Validated: solid-body rotation gives −Ω²; pure strain gives +a².
* **Threshold sensitivity.** For thresholds f·Q_ref (default Q_ref = RMS of Q over the volume, recorded in the output),
  report volume fraction of {Q > threshold} and the number of 6-connected components. The choice of Q_ref is a modelling
  decision; thresholds expressed in multiples of a field statistic make cutouts of different size comparable but not equivalent.
* **Interior statistics.** For non-periodic cutouts the threshold table excludes a 2-voxel margin where one-sided stencils apply.

**Synthetic fields** (`taylor_green`, `abc_flow`) exist only to test the above and to smoke-test ParaView. They are not
turbulence, and outputs derived from them live under `results/synthetic_validation/` with "SYNTHETIC" in the figure titles.

## 4. Reading a real JHTDB cutout

The official cutout page does not state output formats, variable names, index base or axis order. Therefore:

1. Run `python scripts/inspect_dataset.py <file>`; record dataset names, shapes, attributes.
2. Pass `--dataset` and `--axis-order` explicitly to `analyze_spatial_fields.py`; pass `--dx` from the metadata
   (the README implies 2π/1024 ≈ 0.00614 for the full grid, but it does not say that for strided or sub-sampled requests: a stride changes dx).
3. Verify the velocity array: shape (…, 3); mean-flow and RMS plausible; **divergence_rms / vorticity_rms small** (reported
   in `spatial_summary.json`). A wrong axis order typically produces a large divergence.
4. Write what you found here.

## 5. ParaView stage

`paraview/turbulence_pipeline.py` (status: written, **not executed**; see `docs/reproducibility.md`). Filters used and why:
Slice (interior cross-section), Stream Tracer (instantaneous trajectories; *not* material paths), Glyph (direction on a plane
whose colour already encodes magnitude), Gradient (vorticity and Q via VTK's operator), Contour (Q iso-surfaces swept over
several thresholds). The ParaView Q is cross-checked against the NumPy Q on the interior when the input file carries `Q_numpy`;
agreement is expected to O(dx²) and a large mismatch indicates an array-order or spacing error.
