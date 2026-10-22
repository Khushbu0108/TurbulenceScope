"""Tests for turbulencescope.spatial against closed-form results (no real cutout needed)."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from turbulencescope import spatial as sp  # noqa: E402


def linear_field(n, M, dx=0.1):
    """u = M @ x on an n^3 grid; returns u (n,n,n,3)."""
    x = np.arange(n) * dx
    X = np.stack(np.meshgrid(x, x, x, indexing="ij"), axis=-1)
    return X @ np.asarray(M, float).T, dx


class TestGradientsExactCases(unittest.TestCase):
    def test_gradient_of_linear_field_is_exact_including_boundaries(self):
        M = np.array([[1.0, 2.0, -1.0], [0.5, -3.0, 4.0], [2.0, 1.0, 0.25]])
        u, dx = linear_field(8, M)
        A = sp.velocity_gradient(u, dx)
        self.assertTrue(np.allclose(A, M, atol=1e-12))

    def test_solid_body_rotation(self):
        w = 1.7
        M = [[0, -w, 0], [w, 0, 0], [0, 0, 0]]
        u, dx = linear_field(8, M)
        A = sp.velocity_gradient(u, dx)
        omega = sp.vorticity(A)
        self.assertTrue(np.allclose(omega[..., 2], 2 * w))
        self.assertTrue(np.allclose(omega[..., :2], 0))
        self.assertTrue(np.allclose(sp.q_criterion(A), w**2))
        self.assertTrue(np.allclose(sp.lambda2(A), -(w**2)))

    def test_pure_strain_has_negative_Q_and_positive_lambda2(self):
        a = 0.8
        u, dx = linear_field(8, np.diag([a, -a, 0.0]))
        A = sp.velocity_gradient(u, dx)
        self.assertTrue(np.allclose(sp.vorticity(A), 0))
        self.assertTrue(np.allclose(sp.q_criterion(A), -(a**2)))
        self.assertTrue(np.allclose(sp.lambda2(A), a**2))
        self.assertTrue(np.allclose(sp.divergence(A), 0))

    def test_Q_equals_minus_half_AijAji_when_divergence_free(self):
        rng = np.random.default_rng(0)
        B = rng.normal(size=(3, 3))
        B -= np.trace(B) / 3 * np.eye(3)  # trace-free => divergence-free linear field
        u, dx = linear_field(6, B)
        A = sp.velocity_gradient(u, dx)
        q_ref = -0.5 * np.einsum("...ij,...ji->...", A, A)
        self.assertTrue(np.allclose(sp.q_criterion(A), q_ref, atol=1e-10))

    def test_galilean_invariance(self):
        u, dx = sp.taylor_green(16)
        A1 = sp.velocity_gradient(u, dx, periodic=True)
        A2 = sp.velocity_gradient(u + np.array([3.0, -2.0, 1.0]), dx, periodic=True)
        self.assertTrue(np.allclose(A1, A2, atol=1e-10))

    def test_shape_validation(self):
        with self.assertRaises(ValueError):
            sp.velocity_gradient(np.zeros((4, 4, 4)), 1.0)


class TestAnalyticFlows(unittest.TestCase):
    def test_taylor_green_vorticity_matches_closed_form(self):
        n = 48
        u, dx = sp.taylor_green(n)
        (X, Y, Z), _ = sp.grid(n)
        exact = np.stack([-np.cos(X) * np.sin(Y) * np.sin(Z), -np.sin(X) * np.cos(Y) * np.sin(Z),
                          2 * np.sin(X) * np.sin(Y) * np.cos(Z)], axis=-1)
        w_fd = sp.vorticity(sp.velocity_gradient(u, dx, periodic=True))
        w_sp = sp.vorticity(sp.spectral_velocity_gradient(u, 2 * np.pi))
        self.assertLess(np.abs(w_sp - exact).max(), 1e-10)           # spectral: machine precision
        err_fd = np.abs(w_fd - exact).max()
        self.assertLess(err_fd, 0.01)                                # 2nd-order FD at n=48
        self.assertGreater(err_fd, 1e-6)                             # i.e. FD is genuinely approximate

    def test_finite_difference_converges_second_order(self):
        errs = []
        for n in (16, 32, 64):
            u, dx = sp.taylor_green(n)
            (X, Y, Z), _ = sp.grid(n)
            exact = 2 * np.sin(X) * np.sin(Y) * np.cos(Z)
            wz = sp.vorticity(sp.velocity_gradient(u, dx, periodic=True))[..., 2]
            errs.append(np.abs(wz - exact).max())
        rates = [np.log2(errs[i] / errs[i + 1]) for i in range(2)]
        self.assertTrue(all(1.8 < r < 2.2 for r in rates), rates)

    def test_taylor_green_is_divergence_free(self):
        u, dx = sp.taylor_green(32)
        d = sp.divergence(sp.spectral_velocity_gradient(u, 2 * np.pi))
        self.assertLess(np.abs(d).max(), 1e-10)

    def test_abc_flow_is_beltrami(self):
        u, dx = sp.abc_flow(32)
        omega = sp.vorticity(sp.spectral_velocity_gradient(u, 2 * np.pi))
        self.assertLess(np.abs(omega - u).max(), 1e-10)
        self.assertLess(np.abs(sp.divergence(sp.spectral_velocity_gradient(u, 2 * np.pi))).max(), 1e-10)

    def test_nonperiodic_gradient_on_subvolume_close_to_periodic(self):
        u, dx = sp.taylor_green(40)
        Ap = sp.velocity_gradient(u, dx, periodic=True)
        An = sp.velocity_gradient(u, dx, periodic=False)
        interior = (slice(2, -2),) * 3
        self.assertTrue(np.allclose(Ap[interior], An[interior], atol=1e-12))


class TestThresholdSensitivity(unittest.TestCase):
    def test_volume_fraction_monotone_and_components(self):
        Q = np.zeros((20, 20, 20))
        Q[2:5, 2:5, 2:5] = 10.0
        Q[12:14, 12:14, 12:14] = 3.0
        rows = sp.q_threshold_sensitivity(Q, factors=(0.1, 0.5, 1.5), reference=2.0)
        vf = [r["volume_fraction"] for r in rows]
        self.assertTrue(all(vf[i] >= vf[i + 1] for i in range(2)))
        self.assertEqual([r["n_components"] for r in rows], [2, 2, 1])
        self.assertEqual(rows[0]["largest_component_voxels"], 27)

    def test_reference_defaults_to_rms(self):
        Q = np.random.default_rng(1).normal(size=(6, 6, 6))
        r = sp.q_threshold_sensitivity(Q, factors=(1.0,))[0]
        self.assertAlmostEqual(r["reference_Q_rms"], float(np.sqrt(np.mean(Q**2))))

    def test_field_statistics(self):
        s = sp.field_statistics(np.arange(101.0))
        self.assertEqual((s["min"], s["max"], s["median"]), (0.0, 100.0, 50.0))


class TestIO(unittest.TestCase):
    def test_npy_roundtrip_and_shape_check(self):
        u, _ = sp.taylor_green(8)
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "u.npy"
            np.save(p, u)
            self.assertTrue(np.array_equal(sp.load_velocity(p), u))
            np.save(p, u[..., :2])
            with self.assertRaises(ValueError):
                sp.load_velocity(p)

    def test_unsupported_extension(self):
        with self.assertRaises(ValueError):
            sp.load_velocity("nothing.csv")

    def test_legacy_vtk_layout_is_x_fastest(self):
        n = (3, 4, 5)
        s = np.arange(np.prod(n), dtype=float).reshape(n)  # s[x,y,z]
        v = np.stack([s, 2 * s, 3 * s], axis=-1)
        with tempfile.TemporaryDirectory() as d:
            p = sp.write_legacy_vtk(Path(d) / "t.vtk", {"s": s, "v": v}, spacing=(0.5, 0.5, 0.5))
            raw = p.read_bytes()
        head = raw.split(b"POINT_DATA")[0].decode()
        self.assertIn("DIMENSIONS 3 4 5", head)
        self.assertIn("SPACING 0.5 0.5 0.5", head)
        i = raw.index(b"LOOKUP_TABLE default\n") + len(b"LOOKUP_TABLE default\n")
        npts = int(np.prod(n))
        got = np.frombuffer(raw[i:i + 4 * npts], dtype=">f4")
        # x varies fastest: element [x=1,y=0,z=0] is index 1; [x=0,y=1,z=0] is index nx
        self.assertEqual(got[1], s[1, 0, 0])
        self.assertEqual(got[n[0]], s[0, 1, 0])
        self.assertEqual(got[n[0] * n[1]], s[0, 0, 1])


if __name__ == "__main__":
    unittest.main()
