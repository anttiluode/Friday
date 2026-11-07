"""Regression checks for the mirror-lock claims (fast subset of mirror_lock.py)."""
import unittest

import numpy as np

import audit
import mirror_lock as M


class MirrorLockTests(unittest.TestCase):
    def test_dark_modes_are_the_odd_sector_of_a_port_fixing_mirror(self):
        L = M.laplacian(M.base_magnitudes())
        np.testing.assert_allclose(L, audit.laplacian(audit.bridge_phases("golden")), atol=1e-14)
        np.testing.assert_allclose(M.S @ L, L @ M.S, atol=1e-13)
        np.testing.assert_allclose(M.S @ M.S, np.eye(M.N), atol=1e-13)
        np.testing.assert_allclose(M.S[np.ix_(M.PORTS, M.PORTS)], np.eye(5), atol=1e-14)
        D, _ = audit.dark_basis()
        np.testing.assert_allclose(M.S @ D, -D, atol=1e-13)

    def test_stored_words_are_invisible_at_ports_and_readable_by_key(self):
        words = [(0, 0, 0, 0, 0), (1, 0, 1, 1, 0), (1, 1, 1, 1, 1), (0, 1, 0, 0, 1)]
        ref = M.public_reading(M.stored(words[0]))
        for w in words:
            self.assertLess(M.rel(M.public_reading(M.stored(w)), ref), 1e-13)
            bits, margin = M.decode(M.stored(w))
            self.assertEqual(tuple(bits), w)
            self.assertGreater(margin.min(), 0.05)

    def test_public_traffic_cannot_write_and_its_future_is_word_independent(self):
        hist = M.random_public_history(np.random.default_rng(3), 60)
        m = M.run(M.base_magnitudes(), hist, "renormalize")
        self.assertLess(np.abs(M.antisym(m)).max(), 1e-12)
        a, b = M.stored((0, 0, 0, 0, 0)), M.stored((1, 1, 0, 1, 0))
        for drive, f in hist[:30]:
            ra, rb = M.response(a, drive, [f])[0][M.PORTS], M.response(b, drive, [f])[0][M.PORTS]
            self.assertLess(M.rel(ra, rb), 1e-12)
            a, b = M.write(a, drive, f), M.write(b, drive, f)

    def test_key_write_is_local_and_signed(self):
        up = M.antisym(M.key_write(M.base_magnitudes(), 2, +1))
        down = M.antisym(M.key_write(M.base_magnitudes(), 2, -1))
        self.assertGreater(up[2], 1e-3)
        np.testing.assert_allclose(down, -up, atol=1e-13)
        np.testing.assert_allclose(np.delete(up, 2), 0, atol=1e-13)

    def test_flipflop_retains_a_word_under_public_traffic(self):
        w = (1, 0, 0, 1, 1)
        m = M.base_magnitudes()
        for c, bit in enumerate(w):
            m = M.key_write(m, c, +1 if bit else -1)
        for drive, f in M.random_public_history(np.random.default_rng(11), 150):
            m = M.flipflop_step(m, drive, f)
        self.assertEqual(tuple(M.decode(m)[0]), w)


if __name__ == "__main__":
    unittest.main()
