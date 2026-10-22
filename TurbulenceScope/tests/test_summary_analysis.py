"""Tests for summary-file parsing, QC and spectrum analysis. Run: python -m unittest discover -s tests -v
(also collected by pytest if installed)."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from turbulencescope import spectrum_analysis as sa  # noqa: E402
from turbulencescope import summary_io as sio  # noqa: E402

DATA = REPO / "data"


class TestParsers(unittest.TestCase):
    def test_real_history_file(self):
        h = sio.read_energy_history(DATA / "ener_Re_time.txt")
        self.assertEqual(list(h.columns), ["time", "energy", "Re_lambda"])
        self.assertEqual(len(h), 507)
        self.assertAlmostEqual(h["time"].iloc[0], 0.0)
        self.assertAlmostEqual(h["time"].iloc[-1], 10.04)

    def test_real_spectrum_file(self):
        s = sio.read_spectrum(DATA / "spectrum.txt")
        self.assertEqual(len(s), 482)
        self.assertEqual(int(s["k"].iloc[0]), 1)
        self.assertEqual(int(s["k"].iloc[-1]), 482)
        self.assertAlmostEqual(s["E"].iloc[0], 0.299999995142)

    def test_trailing_whitespace_is_tolerated(self):
        txt = " time energy Re_lambda\n-----\n0.0\t1.0\t2.0 \n0.1\t1.1\t2.1 \n"
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.txt"
            p.write_text(txt)
            self.assertEqual(len(sio.read_energy_history(p)), 2)

    def test_malformed_line_raises_not_dropped(self):
        txt = "k E\n-----\n1 0.5\n2 abc\n"
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "s.txt"
            p.write_text(txt)
            with self.assertRaises(sio.SummaryFormatError):
                sio.read_spectrum(p)

    def test_missing_separator_raises(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "s.txt"
            p.write_text("k E\n1 0.5\n")
            with self.assertRaises(sio.SummaryFormatError):
                sio.read_spectrum(p)

    def test_source_files_not_modified(self):
        import hashlib
        before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.glob("*.txt")}
        sio.read_energy_history(DATA / "ener_Re_time.txt")
        sio.read_spectrum(DATA / "spectrum.txt")
        after = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.glob("*.txt")}
        self.assertEqual(before, after)


class TestQC(unittest.TestCase):
    def setUp(self):
        self.h = sio.read_energy_history(DATA / "ener_Re_time.txt")

    def test_duplicates_reported_not_removed(self):
        q = sio.qc_energy_history(self.h)
        self.assertEqual(sorted(q.duplicate_rows["time"].unique().tolist()), [2.04, 3.06, 5.06])
        self.assertEqual(len(q.duplicate_rows), 6)
        self.assertEqual(q.n_duplicate_times_conflicting, 0)
        self.assertEqual(len(self.h), 507)  # input frame untouched

    def test_irregular_steps_detected(self):
        q = sio.qc_energy_history(self.h)
        self.assertEqual(len(q.irregular_steps), 2)
        self.assertTrue(np.allclose(q.irregular_steps["dt"], 0.01))

    def test_no_nonfinite_and_monotonic(self):
        q = sio.qc_energy_history(self.h)
        self.assertEqual(q.n_nonfinite, 0)
        self.assertTrue(q.is_monotonic_nondecreasing)

    def test_conflicting_duplicate_is_flagged(self):
        import pandas as pd
        df = pd.DataFrame({"time": [0.0, 0.02, 0.02, 0.04], "energy": [1, 2, 3, 4], "Re_lambda": [5, 6, 6, 8]})
        self.assertEqual(sio.qc_energy_history(df).n_duplicate_times_conflicting, 1)

    def test_spectrum_qc(self):
        q = sio.qc_spectrum(sio.read_spectrum(DATA / "spectrum.txt"))
        self.assertTrue(q.k_contiguous)
        self.assertEqual((q.k_min, q.k_max, q.n_nonfinite, q.n_nonpositive, q.duplicated_k), (1, 482, 0, 0, 0))


class TestSpectrumMath(unittest.TestCase):
    def test_power_law_fit_recovers_exponent(self):
        k = np.arange(1, 200)
        E = 3.0 * k ** (-5 / 3)
        r = sa.fit_power_law(k, E, 5, 100)
        self.assertAlmostEqual(r["slope"], -5 / 3, places=10)
        self.assertAlmostEqual(r["amplitude_A"], 3.0, places=8)
        self.assertAlmostEqual(r["r_squared"], 1.0, places=10)

    def test_local_slope_of_pure_power_law(self):
        k = np.arange(1, 300)
        s = sa.local_loglog_slope(k, k ** -2.5, window=11)
        self.assertTrue(np.allclose(s[10:-10], -2.5, atol=1e-9))

    def test_local_slope_rejects_even_window(self):
        with self.assertRaises(ValueError):
            sa.local_loglog_slope(np.arange(1, 50), np.arange(1, 50) ** -1.0, window=10)

    def test_candidate_ranges_on_synthetic(self):
        k = np.arange(1, 400)
        E = np.where(k < 100, k ** (-5 / 3), 100 ** (-5 / 3) * np.exp(-(k - 100) / 15))
        df = sa.candidate_scaling_ranges(k, E)
        self.assertGreaterEqual(len(df), 1)
        self.assertLess(df["kmin"].iloc[0], 30)
        self.assertTrue(np.allclose(df["slope"].iloc[0], -5 / 3, atol=0.05))

    def test_fit_needs_points(self):
        with self.assertRaises(ValueError):
            sa.fit_power_law(np.arange(1, 10), np.ones(9), 5, 6)

    def test_compensated_flat_for_five_thirds(self):
        k = np.arange(1, 50)
        self.assertTrue(np.allclose(sa.compensated(k, 2.0 * k ** (-5 / 3)), 2.0))

    def test_integrate_spectrum(self):
        k = np.arange(1, 11)
        self.assertAlmostEqual(sa.integrate_spectrum(k, np.ones(10)), 10.0)

    def test_upturn_detection(self):
        k = np.arange(1, 11)
        E = np.array([5, 4, 3, 2, 1, 0.5, 0.2, 0.1, 0.3, 0.6])
        self.assertEqual(sa.spectral_upturn_onset(k, E), (8, True))

    def test_real_spectrum_sum_matches_mean_energy_within_half_percent(self):
        s = sio.read_spectrum(DATA / "spectrum.txt")
        h = sio.read_energy_history(DATA / "ener_Re_time.txt")
        e_spec = sa.integrate_spectrum(s["k"].to_numpy(), s["E"].to_numpy())
        e_mean = sa.history_summary(h)["energy"]["mean_time_weighted"]
        self.assertLess(abs(e_spec - e_mean) / e_mean, 0.005)


class TestHistoryMath(unittest.TestCase):
    def test_summary_values_on_real_data(self):
        h = sio.read_energy_history(DATA / "ener_Re_time.txt")
        s = sa.history_summary(h)
        self.assertEqual(s["n_records_used"], 504)
        self.assertAlmostEqual(s["energy"]["max"], 0.7298512)
        self.assertAlmostEqual(s["energy"]["t_at_max"], 4.26)
        self.assertAlmostEqual(s["Re_lambda"]["min"], 405.7751)
        self.assertAlmostEqual(s["Re_lambda"]["t_at_min"], 6.84)
        self.assertAlmostEqual(s["Re_lambda"]["max"], 438.2334)

    def test_time_weighted_mean_handles_irregular_spacing(self):
        import pandas as pd
        df = pd.DataFrame({"time": [0.0, 1.0, 3.0], "energy": [0.0, 1.0, 1.0], "Re_lambda": [1.0, 1.0, 1.0]})
        s = sa.history_summary(df, dedupe=False)
        self.assertAlmostEqual(s["energy"]["mean_time_weighted"], (0.5 + 2.0) / 3.0)

    def test_early_time_mean_Re_matches_published_433(self):
        # JHTDB README quotes Re_lambda = 433 averaged over t = 0..2.048 (cited in docs/literature_review.md)
        h = sio.read_energy_history(DATA / "ener_Re_time.txt").drop_duplicates()
        e = h[h["time"] <= 2.048]
        t = e["time"].to_numpy()
        mean = np.trapezoid(e["Re_lambda"], t) / (t[-1] - t[0])
        self.assertAlmostEqual(mean, 433.0, delta=0.1)

    def test_sign_change_interpolation(self):
        self.assertAlmostEqual(sa.sign_change_times([0, 1, 2], [-1, 1, 1])[0], 0.5)


if __name__ == "__main__":
    unittest.main()
