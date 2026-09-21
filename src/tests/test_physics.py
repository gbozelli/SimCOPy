"""Testes de fisica: contra solucao analitica, nao contra outro simulador."""
import sys, unittest
import numpy as np
sys.path.insert(0, ".")
from simcopy.kernels.ssfm import propagate, FiberKernelParams
from simcopy.kernels import birefringence as bir
from simcopy.core.units import beta2_from_dispersion

class TestDispersion(unittest.TestCase):
    def test_gaussian_broadening_matches_closed_form(self):
        """Pulso gaussiano, gamma=0, alpha=0: T(z) = T0*sqrt(1+(z/L_D)^2)."""
        N, fs = 1 << 14, 1e12
        t = (np.arange(N) - N // 2) / fs
        T0 = 20e-12
        A = np.exp(-0.5 * (t / T0) ** 2).astype(complex)[None, :]
        omega = 2 * np.pi * np.fft.fftfreq(N, 1 / fs)
        b2 = beta2_from_dispersion(16e-6, 1550e-9)
        LD = T0**2 / abs(b2)
        for L in (0.5 * LD, 1.0 * LD, 2.0 * LD):
            out = propagate(A, omega, FiberKernelParams(L, b2, 0.0, 0.0),
                            max_step=L / 200)[0]
            I = np.abs(out) ** 2
            Tm = np.sqrt(np.sum(I * t**2) / np.sum(I)) * np.sqrt(2)
            Tt = T0 * np.sqrt(1 + (L / LD) ** 2)
            self.assertLess(abs(Tm - Tt) / Tt, 1e-3,
                            f"L/LD={L/LD:.1f}: medido {Tm*1e12:.3f} ps, "
                            f"teorico {Tt*1e12:.3f} ps")

    def test_energy_conserved_without_attenuation(self):
        N, fs = 1 << 12, 1e12
        rng = np.random.default_rng(0)
        A = (rng.standard_normal((2, N)) + 1j * rng.standard_normal((2, N))) * 1e-2
        omega = 2 * np.pi * np.fft.fftfreq(N, 1 / fs)
        p = FiberKernelParams(80e3, -2e-26, 1.3e-3, 0.0, pmd_coefficient=0.0)
        out = propagate(A, omega, p, rng=np.random.default_rng(1))
        e0, e1 = np.sum(np.abs(A) ** 2), np.sum(np.abs(out) ** 2)
        self.assertLess(abs(e1 - e0) / e0, 1e-10)

    def test_attenuation_matches_beer_lambert(self):
        N = 1 << 10
        A = np.ones((1, N), dtype=complex)
        omega = 2 * np.pi * np.fft.fftfreq(N, 1e-12)
        a = 0.2e-3 / 4.343
        out = propagate(A, omega, FiberKernelParams(100e3, 0.0, 0.0, a))
        self.assertAlmostEqual(np.mean(np.abs(out) ** 2), np.exp(-a * 100e3), places=12)

class TestBirefringence(unittest.TestCase):
    def test_su2_is_haar_uniform_on_poincare_sphere(self):
        rng = np.random.default_rng(11)
        S = []
        for _ in range(20000):
            U = bir.random_su2(rng)
            ax, ay = U[0, 0], U[1, 0]          # entrada linear em X
            s0 = abs(ax) ** 2 + abs(ay) ** 2
            S.append([(abs(ax) ** 2 - abs(ay) ** 2) / s0,
                      2 * np.real(ax * np.conj(ay)) / s0,
                      2 * np.imag(ax * np.conj(ay)) / s0])
        S = np.array(S)
        for k, lab in enumerate("123"):
            h = np.histogram(S[:, k], bins=8, range=(-1, 1))[0] / (len(S) / 8)
            self.assertLess(np.abs(h - 1).max(), 0.06, f"S{lab} nao uniforme: {h}")

    def test_su2_is_unitary(self):
        rng = np.random.default_rng(5)
        for _ in range(200):
            U = bir.random_su2(rng)
            np.testing.assert_allclose(U @ U.conj().T, np.eye(2), atol=1e-12)

    def test_dgd_scales_with_sqrt_length(self):
        """DGD total ~ D_PMD*sqrt(L): a assinatura do modelo de passo grosso."""
        N, fs = 1 << 12, 2e12
        omega = 2 * np.pi * np.fft.fftfreq(N, 1 / fs)
        dpmd = 0.1e-12 / 31.62
        dw = 2 * np.pi * 5e9
        def mean_dgd(L, n_real=60):
            out = []
            for r in range(n_real):
                J = np.empty((2, 2, 2), dtype=complex)
                for iw, w0 in enumerate((0.0, dw)):
                    for k, v in enumerate(((1, 0), (0, 1))):
                        A = np.zeros((2, N), dtype=complex)
                        A[0] = v[0] * np.exp(1j * w0 * np.arange(N) / fs)
                        A[1] = v[1] * np.exp(1j * w0 * np.arange(N) / fs)
                        o = propagate(A, omega,
                            FiberKernelParams(L, 0.0, 0.0, 0.0, dpmd, 50.0),
                            rng=np.random.default_rng(1000 + r), max_step=1e3)
                        J[iw, :, k] = o[:, N // 2] * np.exp(-1j * w0 * (N // 2) / fs)
                M = np.linalg.solve(J[0], J[1])
                ev = np.linalg.eigvals(M)
                out.append(abs(np.angle(ev[0] / ev[1]) / dw))
            return np.mean(out)
        d80, d320 = mean_dgd(80e3), mean_dgd(320e3)
        self.assertGreater(d80, 0.0)
        self.assertLess(abs(d320 / d80 - 2.0), 0.35,
                        f"escala em sqrt(L) falhou: {d320/d80:.3f} (esperado ~2)")

if __name__ == "__main__":
    unittest.main(verbosity=2)
