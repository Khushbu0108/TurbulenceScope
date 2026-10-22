"""Spatial-field analysis for a 3-D velocity cutout (NumPy only).

Array convention used throughout this module: ``u`` has shape ``(nx, ny, nz, 3)`` with
components ordered (u_x, u_y, u_z) and axes ordered (x, y, z). ``dx`` is the grid spacing
per axis. **The axis order and spacing of a real JHTDB file must be confirmed from the
returned file's own metadata** (see docs/methodology.md); ``load_velocity`` does not guess.

Definitions (A_ij = du_i/dx_j is the velocity-gradient tensor)
  vorticity           omega = curl u,  omega_i = eps_ijk A_kj
  strain-rate tensor  S = (A + A^T)/2
  rotation tensor     Omega = (A - A^T)/2
  Q-criterion         Q = (1/2)(||Omega||_F^2 - ||S||_F^2)   (Hunt, Wray & Moin 1988)
  lambda_2            second-largest eigenvalue of S^2 + Omega^2 (Jeong & Hussain 1995)
For incompressible flow, Q = -(1/2) A_ij A_ji.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- gradients


def velocity_gradient(u: np.ndarray, dx, periodic: bool = False) -> np.ndarray:
    """Return A with shape (nx, ny, nz, 3, 3), A[..., i, j] = d u_i / d x_j.

    periodic=False: second-order central differences in the interior and second-order
        one-sided differences on the boundary faces (``numpy.gradient(edge_order=2)``).
        Use this for JHTDB *cutouts*, which are sub-volumes and are not periodic.
    periodic=True: second-order central differences with wrap-around. Only valid when the
        array is a whole periodic box.
    """
    if u.ndim != 4 or u.shape[-1] != 3:
        raise ValueError("u must have shape (nx, ny, nz, 3)")
    dx = np.broadcast_to(np.asarray(dx, float), (3,))
    A = np.empty(u.shape[:3] + (3, 3), dtype=np.result_type(u.dtype, np.float32))
    for i in range(3):
        for j in range(3):
            if periodic:
                A[..., i, j] = (np.roll(u[..., i], -1, axis=j) - np.roll(u[..., i], 1, axis=j)) / (2 * dx[j])
            else:
                A[..., i, j] = np.gradient(u[..., i], dx[j], axis=j, edge_order=2)
    return A


def spectral_velocity_gradient(u: np.ndarray, L) -> np.ndarray:
    """Fourier-spectral gradient for a whole periodic box of side(s) ``L``.

    Independent of the finite-difference route and used as a cross-check on periodic
    synthetic fields. The Nyquist mode is zeroed for the odd-derivative operator.
    """
    L = np.broadcast_to(np.asarray(L, float), (3,))
    n = u.shape[:3]
    A = np.empty(u.shape[:3] + (3, 3), dtype=float)
    ks = []
    for ax in range(3):
        k = 2 * np.pi * np.fft.fftfreq(n[ax], d=L[ax] / n[ax])
        if n[ax] % 2 == 0:
            k[n[ax] // 2] = 0.0
        shape = [1, 1, 1]
        shape[ax] = n[ax]
        ks.append(k.reshape(shape))
    for i in range(3):
        uh = np.fft.fftn(u[..., i])
        for j in range(3):
            A[..., i, j] = np.real(np.fft.ifftn(1j * ks[j] * uh))
    return A


# ------------------------------------------------------------------------- derived fields


def vorticity(A: np.ndarray) -> np.ndarray:
    """omega = (A_zy - A_yz, A_xz - A_zx, A_yx - A_xy)."""
    return np.stack([A[..., 2, 1] - A[..., 1, 2],
                     A[..., 0, 2] - A[..., 2, 0],
                     A[..., 1, 0] - A[..., 0, 1]], axis=-1)


def strain_rotation(A: np.ndarray):
    At = np.swapaxes(A, -1, -2)
    return 0.5 * (A + At), 0.5 * (A - At)


def q_criterion(A: np.ndarray) -> np.ndarray:
    S, W = strain_rotation(A)
    return 0.5 * (np.sum(W * W, axis=(-1, -2)) - np.sum(S * S, axis=(-1, -2)))


def lambda2(A: np.ndarray, chunk: int = 1 << 18) -> np.ndarray:
    """Second-largest eigenvalue of S^2 + Omega^2 (symmetric 3x3 -> eigvalsh, ascending)."""
    S, W = strain_rotation(A)
    M = S @ S + W @ W
    flat = M.reshape(-1, 3, 3)
    out = np.empty(flat.shape[0])
    for a in range(0, flat.shape[0], chunk):
        ev = np.linalg.eigvalsh(flat[a:a + chunk])
        out[a:a + chunk] = ev[:, 1]
    return out.reshape(A.shape[:3])


def divergence(A: np.ndarray) -> np.ndarray:
    return A[..., 0, 0] + A[..., 1, 1] + A[..., 2, 2]


def speed(u: np.ndarray) -> np.ndarray:
    return np.sqrt(np.sum(u.astype(float) ** 2, axis=-1))


# --------------------------------------------------------------------------- statistics


def field_statistics(a: np.ndarray) -> dict:
    a = np.asarray(a, float).ravel()
    q = np.percentile(a, [1, 5, 25, 50, 75, 95, 99])
    return {"n": int(a.size), "min": float(a.min()), "max": float(a.max()), "mean": float(a.mean()),
            "std": float(a.std()), "rms": float(np.sqrt(np.mean(a**2))),
            "p01": float(q[0]), "p05": float(q[1]), "p25": float(q[2]), "median": float(q[3]),
            "p75": float(q[4]), "p95": float(q[5]), "p99": float(q[6])}


def q_threshold_sensitivity(Q: np.ndarray, factors=(0.5, 1.0, 2.0, 4.0), reference: float | None = None,
                            min_voxels: int = 1) -> list[dict]:
    """Volume fraction and connected-component count of {Q > f * Q_ref}.

    ``reference`` defaults to the RMS of Q over the volume. Thresholds that are fixed
    multiples of a field statistic make cutouts of different size comparable, but the
    choice of reference is itself a modelling decision and is recorded in the output.
    Connectivity is 6-connected (face neighbours).
    """
    from scipy import ndimage

    ref = float(np.sqrt(np.mean(np.asarray(Q, float) ** 2))) if reference is None else float(reference)
    rows = []
    for f in factors:
        thr = f * ref
        mask = Q > thr
        lab, n = ndimage.label(mask)
        if n:
            sizes = np.bincount(lab.ravel())[1:]
            kept = int((sizes >= min_voxels).sum())
            largest = int(sizes.max())
        else:
            kept, largest = 0, 0
        rows.append({"factor": float(f), "threshold": thr, "reference_Q_rms": ref,
                     "volume_fraction": float(mask.mean()), "n_components": kept,
                     "largest_component_voxels": largest,
                     "largest_component_fraction_of_marked": float(largest / mask.sum()) if mask.any() else 0.0})
    return rows


# ------------------------------------------------------------------------------ I/O


def load_velocity(path: str | Path, dataset: str | None = None, time_index: int = 0,
                  axis_order: str = "zyxc") -> np.ndarray:
    """Load a velocity cutout and return it as (nx, ny, nz, 3).

    Supported: ``.npy`` (already (nx,ny,nz,3)), and HDF5 (``.h5``/``.hdf5``, requires h5py).
    For HDF5 you must say how the stored array is laid out via ``axis_order``:
    ``"zyxc"`` means stored as (z, y, x, component); ``"xyzc"`` means (x, y, z, component).
    The default is a documented assumption to **verify against the file's metadata**;
    the loader prints the dataset names/shapes it finds when ``dataset`` is None and
    refuses to guess among several candidates.
    """
    p = Path(path)
    if p.suffix == ".npy":
        u = np.load(p)
        if u.ndim != 4 or u.shape[-1] != 3:
            raise ValueError(f"{p}: expected shape (nx,ny,nz,3), got {u.shape}")
        return u
    if p.suffix in {".h5", ".hdf5"}:
        try:
            import h5py
        except ImportError as e:  # pragma: no cover
            raise ImportError("h5py is required to read HDF5 files: pip install h5py") from e
        with h5py.File(p, "r") as f:
            cands = []
            f.visititems(lambda n, o: cands.append((n, o.shape)) if isinstance(o, h5py.Dataset) else None)
            if dataset is None:
                vec = [(n, s) for n, s in cands if len(s) == 4 and s[-1] == 3]
                if len(vec) != 1:
                    raise ValueError(f"Cannot choose the velocity dataset automatically. Datasets: {cands}. "
                                     f"Pass dataset=<name>.")
                dataset = vec[0][0]
            arr = np.asarray(f[dataset])
        if arr.ndim != 4 or arr.shape[-1] != 3:
            raise ValueError(f"{dataset}: expected 4-D (.., 3) array, got {arr.shape}")
        if axis_order == "zyxc":
            arr = np.transpose(arr, (2, 1, 0, 3))
        elif axis_order != "xyzc":
            raise ValueError("axis_order must be 'zyxc' or 'xyzc'")
        return np.ascontiguousarray(arr)
    raise ValueError(f"Unsupported velocity file type: {p.suffix}")


def write_legacy_vtk(path: str | Path, fields: dict[str, np.ndarray], spacing=(1.0, 1.0, 1.0),
                     origin=(0.0, 0.0, 0.0), title: str = "TurbulenceScope field") -> Path:
    """Write point-data fields to a legacy binary VTK STRUCTURED_POINTS file (ParaView-readable).

    ``fields`` maps names to arrays shaped (nx,ny,nz) [scalar] or (nx,ny,nz,3) [vector]. All
    must share (nx,ny,nz). Data are float32, big-endian, x-fastest as VTK requires. Dependency-free
    so that the file can be generated without VTK installed.
    """
    p = Path(path)
    shape = None
    for name, a in fields.items():
        s = a.shape[:3]
        if shape is not None and s != shape:
            raise ValueError(f"field {name} has shape {s}, expected {shape}")
        shape = s
    nx, ny, nz = shape
    with open(p, "wb") as f:
        f.write(b"# vtk DataFile Version 3.0\n")
        f.write(title.replace("\n", " ")[:250].encode("ascii", "replace") + b"\n")
        f.write(b"BINARY\nDATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {nx} {ny} {nz}\n".encode())
        f.write(f"ORIGIN {origin[0]:.9g} {origin[1]:.9g} {origin[2]:.9g}\n".encode())
        f.write(f"SPACING {spacing[0]:.9g} {spacing[1]:.9g} {spacing[2]:.9g}\n".encode())
        f.write(f"POINT_DATA {nx * ny * nz}\n".encode())
        for name, a in fields.items():
            ncomp = 1 if a.ndim == 3 else a.shape[-1]
            kind = "SCALARS" if ncomp == 1 else "VECTORS"
            # reorder to x-fastest: transpose (x,y,z[,c]) -> (z,y,x[,c]) then C-order flatten
            order = (2, 1, 0) if ncomp == 1 else (2, 1, 0, 3)
            buf = np.ascontiguousarray(np.transpose(a, order)).astype(">f4").tobytes()
            if kind == "SCALARS":
                f.write(f"SCALARS {name} float 1\nLOOKUP_TABLE default\n".encode())
            else:
                f.write(f"VECTORS {name} float\n".encode())
            f.write(buf)
            f.write(b"\n")
    return p


# ------------------------------------------------------- synthetic validation fields
# These exist ONLY to validate the code against closed-form answers and to smoke-test the
# ParaView pipeline when no real cutout is available. They are not turbulence.


def grid(n: int, L: float = 2 * np.pi):
    x = np.arange(n) * (L / n)
    return np.meshgrid(x, x, x, indexing="ij"), L / n


def taylor_green(n: int, L: float = 2 * np.pi):
    """Taylor-Green vortex: u=(sin x cos y cos z, -cos x sin y cos z, 0). Divergence-free.

    Analytic gradient entries are used in tests (see tests/test_spatial_analysis.py).
    """
    (X, Y, Z), dx = grid(n, L)
    u = np.stack([np.sin(X) * np.cos(Y) * np.cos(Z), -np.cos(X) * np.sin(Y) * np.cos(Z),
                  np.zeros_like(X)], axis=-1)
    return u, dx


def abc_flow(n: int, A=1.0, B=np.sqrt(2 / 3), C=np.sqrt(1 / 3), L: float = 2 * np.pi):
    """Arnold-Beltrami-Childress flow: Beltrami (omega = u), divergence-free, periodic.

    u = (A sin z + C cos y, B sin x + A cos z, C sin y + B cos x). Analytically curl u = u.
    """
    (X, Y, Z), dx = grid(n, L)
    u = np.stack([A * np.sin(Z) + C * np.cos(Y), B * np.sin(X) + A * np.cos(Z),
                  C * np.sin(Y) + B * np.cos(X)], axis=-1)
    return u, dx
