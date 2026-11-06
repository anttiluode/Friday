"""Regression checks for mathematical claims, not just printed output."""
import contextlib
import importlib
import importlib.util
import io
from pathlib import Path
import runpy
import unittest

import numpy as np


class AuditTests(unittest.TestCase):
    def audit(self):
        self.assertIsNotNone(importlib.util.find_spec("audit"), "audit module is missing")
        return importlib.import_module("audit")

    def test_matches_original_matrices(self):
        a = self.audit()
        with contextlib.redirect_stdout(io.StringIO()):
            original = runpy.run_path(str(Path(__file__).parent / "original/magnetic.py"))
        for mode, golden in [("rational", False), ("golden", True)]:
            expected = original["build_magnetic_laplacian"](golden_angle_wiring=golden)
            np.testing.assert_allclose(a.laplacian(a.bridge_phases(mode)), expected, atol=1e-14)

    def test_analytic_dark_modes_survive_arbitrary_outer_phases(self):
        a = self.audit()
        D, values = a.dark_basis()
        np.testing.assert_allclose(D.conj().T @ D, np.eye(15), atol=1e-14)
        np.testing.assert_allclose(D[a.PORTS], 0, atol=1e-14)
        rng = np.random.default_rng(42)
        for _ in range(5):
            L = a.laplacian(rng.uniform(-10, 10, 5))
            np.testing.assert_allclose(L @ D, D * values, atol=1e-13)
            self.assertAlmostEqual(np.linalg.eigvalsh(L)[0], 2-np.sqrt(3), places=12)

    def test_same_outer_flux_preserves_spectrum(self):
        a = self.audit()
        p = a.bridge_phases("golden")
        left = np.linalg.eigvalsh(a.laplacian(p))
        right = np.linalg.eigvalsh(a.laplacian([0, 0, 0, 0, sum(p)]))
        np.testing.assert_allclose(left, right, atol=1e-13)

    def test_sensor_placement_and_rational_control(self):
        a = self.audit()
        for mode, expected in [("rational", (15, 27, 30)),
                               ("golden", (15, 30, 30)),
                               ("rational_11", (15, 30, 30))]:
            L = a.laplacian(a.bridge_phases(mode))
            sensors = [a.PORTS, a.PORTS+1, np.sort(np.r_[a.PORTS, a.PORTS+1])]
            self.assertEqual(tuple(a.observable_dimension(L, s) for s in sensors), expected)

    def test_internal_opening_exposes_modes_but_weakly(self):
        a = self.audit()
        closed = a.laplacian(a.bridge_phases("golden"))
        weak = a.laplacian(a.bridge_phases("golden"), delta=1e-4)
        opened = a.laplacian(a.bridge_phases("golden"), delta=0.1)
        self.assertEqual(a.observable_dimension(weak, a.PORTS), 30)
        self.assertLess(a.finite_singular_values(closed, a.PORTS)[-1], 1e-12)
        self.assertTrue(7e-6 < a.finite_singular_values(weak, a.PORTS)[-1] < 8e-6)
        self.assertTrue(0.007 < a.finite_singular_values(opened, a.PORTS)[-1] < 0.008)


if __name__ == "__main__":
    unittest.main()
