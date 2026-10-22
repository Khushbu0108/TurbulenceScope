# Literature review

**How this was done.** Abstracts and metadata for Papers A and B and the arXiv version of Paper C were read through a
web-fetch tool on 2026-10-21; JHTDB's dataset, cutout, citing and README pages were read the same way. **Full texts of
all three papers were not read**, so statements below are limited to what the abstracts/pages support. Where a section
says "not available to me", the detail would be in the full paper. Nothing here is a reproduction of a published result,
and no published figure is reused.

## Paper A: McKeown, Pumir, Rubinstein, Brenner & Ostilla-Mónico

*Energy transfer and vortex structures: Visualizing the incompressible turbulent energy cascade.* arXiv:2204.09023
(v2, 2023). https://arxiv.org/abs/2204.09023

* **Question.** How does kinetic energy move between scales, and does a simple model configuration (two antiparallel
  vortex tubes) capture the energy transfer seen in statistically steady homogeneous isotropic turbulence (HIT)?
* **Methods (per abstract/summary).** Velocity is band-filtered in logarithmically spaced wavenumber bands; energy
  transfer is analysed in Fourier space and mapped back to physical space, then compared with vorticity and strain filtered to
  the same bands. A steady HIT simulation is compared with a transient vortex-tube reconnection flow.
* **Findings that can inform this project.** In HIT, intense energy-transfer regions were reported to show no significant
  correlation with filtered vorticity or strain, and no simple self-similar mechanism was found; the vortex-tube flow was
  more orderly; the authors caution that transient flows may not represent HIT.
* **How it informs TurbulenceScope.** It is a direct reason for the project's central caution: *high vorticity or a
  Q-isosurface is not, by itself, an explanation of energy transfer.* The project therefore reports where structures are
  and how threshold-dependent they are, and makes no claim about the cascade mechanism.
* **Differences / limitations.** This project does **no** band filtering, computes **no** energy-transfer fields, and uses
  a small cutout of a different simulation. It does not reproduce or test any claim of Paper A.

## Paper B: Bürger, Treib, Westermann, Werner, Lalescu, Szalay, Meneveau & Eyink

*Vortices within vortices: hierarchical nature of vortex tubes in turbulence.* arXiv:1210.3325 (v3, 2023).
https://arxiv.org/abs/1210.3325

* **Question.** Are small-scale vortical structures simply advected by larger ones, or deformed by them?
* **Methods (per abstract).** The JHU turbulence database was combined with a GPU-oriented interactive visualization tool to
  produce high-quality videos for examining structures at several scales. The vortex-identification criterion, specific
  queries and quantitative details are *not available to me* (abstract only).
* **Findings that can inform this project.** The abstract reports that small structures are sheared and twisted by
  larger-scale motions rather than passively carried, and relates this to the cascade.
* **How it informs TurbulenceScope.** It motivates visual inspection of nested/tube-like structure and shows JHTDB is a
  legitimate source for vortex visualization. It also motivates scale-aware thinking: a single Q threshold selects a
  *range* of scales, so the project sweeps thresholds.
* **Differences / limitations.** Paper B uses a dedicated interactive renderer and, per the abstract, time-resolved
  video; this project uses ParaView on a small cutout, static or few-step. The project cannot confirm hierarchical
  structure and does not claim to.

## Paper C: Li, Perlman, Wan, Yang, Meneveau, Burns, Chen, Szalay & Eyink

*A public turbulence database cluster and applications to study Lagrangian evolution of velocity increments in
turbulence.* Journal of Turbulence 9, No. 31 (2008). DOI 10.1080/14685240802376389. (arXiv:0804.1703 read for the
abstract; the publisher page could not be fetched.)

* **Question.** How can a very large DNS be made publicly usable, and what does it reveal about the Lagrangian
  evolution of velocity increments (pressure vs. viscous effects on intermittency)?
* **Methods (per abstract).** A database cluster storing forced isotropic 1024³ DNS data, accessed through web services;
  users run code locally that retrieves only needed data and uses built-in differentiation and interpolation; accuracy
  tests of those operations; conditional averages of velocity-increment evolution.
* **Findings that can inform this project.** Practical: the dataset is designed to be queried remotely without
  downloading it, and the built-in derivative operators have finite accuracy. JHTDB's README for this dataset states
  about 7 % RMS error for database-side gradients (finite-difference/spline rather than the DNS's spectral operators).
* **How it informs TurbulenceScope.** It is the dataset's reference publication and a required citation. The gradient-error
  statement is why `docs/scientific_limitations.md` treats vorticity/Q from a cutout as *approximate* and why this project
  tests its own finite-difference gradients against closed-form flows instead of assuming accuracy.
* **Differences / limitations.** The project does not study Lagrangian evolution or velocity increments.

## JHTDB documentation used

* Dataset page (https://turbulence.pha.jhu.edu/Forced_isotropic_turbulence.aspx): 1024³ nodes; pseudo-spectral DNS; velocity
  (3 components) and pressure; coarse set of 5,028 frames, t = 0 … 10.056, saved every 10 DNS steps; fine set for a very
  short range (t ≈ 0.0002 … 0.0198); Re_λ "around 433"; forcing holds energy in |k| ≤ 2 constant; downloadable
  energy/Re_λ time history and time-averaged radial spectrum; cutout and Python/MATLAB/C/.NET access.
* Dataset README (https://turbulence.idies.jhu.edu/docs/isotropic/README-isotropic1024.pdf): domain [0, 2π]³, ν = 0.000185,
  DNS time step 0.0002, stored dt = 0.002. **Discrepancy:** the README states 5,024 frames to t = 10.048, whereas the
  dataset page states 5,028 frames to t = 10.056. The supplied `ener_Re_time.txt` ends at t = 10.04 with 0.02 spacing.
  Resolve against the returned data before relying on either figure.
* Cutout page (https://turbulence.idies.jhu.edu/database/cutout): a personal turbulence token is entered in the form;
  with the default token a cutout cannot exceed 4,096 cells; results go to private SciServer scratch storage and are
  deleted after at least 3 days (currently 180); cutouts above ~2 GB may be unstable. **Output formats, dataset/variable
  names, index base and axis order are not stated on that page**, so this project does not assume them.

## Bibliography

1. McKeown, R., Pumir, A., Rubinstein, S. M., Brenner, M. P., & Ostilla-Mónico, R. (2022; rev. 2023). Energy transfer and
   vortex structures: Visualizing the incompressible turbulent energy cascade. arXiv:2204.09023.
2. Bürger, K., Treib, M., Westermann, R., Werner, S., Lalescu, C. C., Szalay, A., Meneveau, C., & Eyink, G. L. (2012; rev.
   2023). Vortices within vortices: hierarchical nature of vortex tubes in turbulence. arXiv:1210.3325.
3. Li, Y., Perlman, E., Wan, M., Yang, Y., Meneveau, C., Burns, R., Chen, S., Szalay, A., & Eyink, G. (2008). A public
   turbulence database cluster and applications to study Lagrangian evolution of velocity increments in turbulence.
   *Journal of Turbulence*, 9(31). https://doi.org/10.1080/14685240802376389
4. Perlman, E., Burns, R., Li, Y., & Meneveau, C. (2007). Data exploration of turbulence simulations using a database
   cluster. *Supercomputing SC07*, ACM/IEEE.
5. Hunt, J. C. R., Wray, A. A., & Moin, P. (1988). Eddies, streams, and convergence zones in turbulent flows. *Center for
   Turbulence Research Report CTR-S88*. (Q-criterion; cited from general knowledge, not re-read for this project.)
6. Jeong, J., & Hussain, F. (1995). On the identification of a vortex. *Journal of Fluid Mechanics*, 285, 69–94.
   (λ₂ criterion; cited from general knowledge, not re-read for this project.)
7. JHTDB, Isotropic Turbulence on 1024³ Grid (isotropic1024coarse). https://doi.org/10.7281/T1KK98XB
