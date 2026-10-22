# Data provenance

## Files in this directory

| File | Content | Rows | Notes |
|---|---|---|---|
| `ener_Re_time.txt` | columns `time energy Re_lambda` | 507 data rows (504 unique) | t = 0 … 10.04; contains 3 exact duplicate records and an irregular 0.01 step near t = 2.05 (see `results/tables/`) |
| `spectrum.txt` | columns `k E(k)` | 482 | integer shells k = 1 … 482 |

Both are summary files **supplied by the project owner**. The JHTDB dataset page for Forced Isotropic Turbulence lists
a time history of total kinetic energy and Re_λ and a radial energy spectrum E(k) as downloadable files; the owner's
copies are assumed to be those files. The exact download date/URL was not recorded: *please add it here*. SHA-256 of the
copies used for every result in this repository are written to `results/tables/summary_results.json`.

These files are **not** velocity fields and cannot be used to reconstruct one.

## Source dataset

Johns Hopkins Turbulence Databases (JHTDB), *Forced Isotropic Turbulence* (listed by JHTDB as
"Isotropic Turbulence on 1024³ Grid, isotropic1024coarse").

* Dataset page: https://turbulence.pha.jhu.edu/Forced_isotropic_turbulence.aspx
* Database: https://turbulence.idies.jhu.edu/database
* Dataset DOI listed on JHTDB's citing page: https://doi.org/10.7281/T1KK98XB
* License: Open Data Commons Attribution License (ODC-By) per JHTDB's citing page. Attribution is required.

## Required citations (from JHTDB's citing page, https://turbulence.idies.jhu.edu/citing)

* Y. Li, E. Perlman, M. Wan, Y. Yang, C. Meneveau, R. Burns, S. Chen, A. Szalay & G. Eyink. "A public turbulence database
  cluster and applications to study Lagrangian evolution of velocity increments in turbulence." *Journal of Turbulence* 9,
  No. 31, 2008.
* E. Perlman, R. Burns, Y. Li, and C. Meneveau. "Data Exploration of Turbulence Simulations using a Database Cluster."
  *Supercomputing SC07*, ACM/IEEE, 2007.

Suggested acknowledgement for talks/posters: "Data obtained from the JHTDB at http://turbulence.idies.jhu.edu."

> Verified by the author of this repository on 2026-10-21 by reading the pages above through a web-fetch tool. JHTDB may
> change its requirements; check the citing page before publishing.

## Spatial cutouts (not included)

No spatial velocity data are committed. See `docs/reproducibility.md` for how to obtain a cutout with your own JHTDB
token. Cutouts, HDF5/XDMF files and tokens are excluded by `.gitignore`.
