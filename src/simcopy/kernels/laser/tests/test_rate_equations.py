"""Equacoes de taxa do laser contra formas fechadas e estatistica de Langevin."""

import dataclasses

import numpy as np
import pytest

from simcopy.kernels.laser.conversions import facet_power
from simcopy.kernels.laser.rate_equations import (
    RateEquationParams,
    derivatives,
    gain,
    integrate,
    steady_state,
)
from simcopy.testing import approx

# Parametros do professor em densidade. O R_sp pequeno e mantido nestes cenarios
# deterministas para reproduzir as referencias historicas do notebook.
PROFESSOR = RateEquationParams(
    gain_coefficient=3.526970094117647e-12,
    transparency_density=1.25e24,
    photon_lifetime=2.8931579685554005e-12,
    carrier_lifetime=2.1e-9,
    volume=1.8e-16,
    spontaneous_emission_rate=3.968253968253968e-16,
    linewidth_enhancement=3.0,
)
N_TH_PROFESSOR = 1.348e24
I_TH_AGRAWAL = 0.018612835881208975
S_30MA = 1.1524800754727474e21
S_31MA = 1.2528004738871203e21


def agrawal(volume, r_sp_number=0.0):
    """Exercicios 3.7--3.10 de Agrawal convertidos de numero para densidade."""
    return RateEquationParams(
        gain_coefficient=6e3 * volume,
        transparency_density=1e8 / volume,
        photon_lifetime=1.2593408321736882e-12,
        carrier_lifetime=2e-9,
        volume=volume,
        spontaneous_emission_rate=r_sp_number / volume,
    )


def test_gain():
    assert gain(1.3e24, 1e20, PROFESSOR) == approx(176348504705.88168, rel=1e-12)
    compressed = dataclasses.replace(PROFESSOR, gain_compression=1e-22)
    assert gain(1.3e24, 1e20, compressed) == approx(174585019658.82288, rel=1e-12)


def test_derivatives():
    # Referencia calculada diretamente das Eqs. 3.3.1, 3.3.2 e 3.3.16.
    dn, ds, dphi = derivatives(1.3e24, 1e20, 30e-3, PROFESSOR)
    assert dn == approx(4.0356904289192e32, rel=1e-9)
    assert ds == approx(-1.6929456451764774e31, rel=1e-9)
    assert dphi == approx(-253941846776.47162, rel=1e-9)


def test_derivatives_fase_sem_compressao():
    # A fase usa o ganho sem o fator (1 - eps*S), conforme a Eq. 3.3.16.
    compressed = dataclasses.replace(PROFESSOR, gain_compression=1e-22)
    assert derivatives(1.3e24, 1e20, 30e-3, compressed)[2] == approx(
        -253941846776.47162, rel=1e-9
    )


@pytest.mark.parametrize(
    ("current", "expected_s"),
    [
        (20e-3, 1.4927609132902528e20),
        (30e-3, 1.1524800754727474e21),
    ],
)
def test_steady_state_acima_do_limiar(current, expected_s):
    n, s = steady_state(current, PROFESSOR)
    assert s == approx(expected_s, rel=1e-12)
    assert n == approx(N_TH_PROFESSOR, rel=1e-12)


def test_steady_state_abaixo_do_limiar():
    n, s = steady_state(10e-3, PROFESSOR)
    assert n == approx(7.281760586870891e23, rel=1e-12)
    assert s == approx(1.815220279146239e-28, rel=1e-9)


def test_steady_state_sem_emissao_espontanea_abaixo_do_limiar():
    without_spontaneous_emission = dataclasses.replace(
        PROFESSOR, spontaneous_emission_rate=0.0
    )
    n, s = steady_state(10e-3, without_spontaneous_emission)
    assert s == 0.0
    assert n == approx(7.281760586870891e23, rel=1e-12)


def test_steady_state_agrawal_exercicio_3_8():
    params = agrawal(volume=1e-16)
    n, s = steady_state(2 * I_TH_AGRAWAL, params)
    assert s * params.volume == approx(146300.37494201775, rel=1e-12)
    assert n * params.volume == approx(232344368.10802943, rel=1e-12)


def test_steady_state_nao_depende_do_volume_escolhido():
    for volume in (1e-16, 1.8e-16, 1.0):
        params = agrawal(volume)
        photons = steady_state(2 * I_TH_AGRAWAL, params)[1] * volume
        assert photons == approx(146300.37494201775, rel=1e-12)


def test_steady_state_agrawal_exercicio_3_10_com_emissao_espontanea():
    params = agrawal(volume=1e-16, r_sp_number=1588132417296.3528)
    n, s = steady_state(2 * I_TH_AGRAWAL, params)
    assert s * params.volume == approx(146303.5141263201, rel=1e-10)
    assert n * params.volume == approx(232342558.93250957, rel=1e-10)


def test_steady_state_zera_as_derivadas_com_compressao():
    compressed = dataclasses.replace(PROFESSOR, gain_compression=1e-22)
    n, s = steady_state(30e-3, compressed)
    dn, ds, _ = derivatives(n, s, 30e-3, compressed)
    assert abs(dn) * compressed.carrier_lifetime / n < 1e-12
    assert abs(ds) * compressed.photon_lifetime / s < 1e-12
    assert n > N_TH_PROFESSOR
    assert s < S_30MA


def test_curva_potencia_corrente():
    def power(current):
        photon_density = steady_state(current, PROFESSOR)[1]
        return facet_power(
            photon_density,
            1.8e-16,
            88174252.35294119,
            2410.0,
            1.2815779723541474e-19,
        )

    assert power(20e-3) == approx(0.0003658782784898758, rel=1e-9)
    assert power(30e-3) == approx(0.0028247485732891897, rel=1e-9)
    assert (power(30e-3) - power(20e-3)) / 10e-3 == approx(0.2458870294799314, rel=1e-9)


def test_integrate_estacionario_e_ponto_fixo():
    n0, s0 = steady_state(30e-3, PROFESSOR)
    n, s, phi = integrate(np.full(2000, 30e-3), 1e-12, PROFESSOR, n0, s0)
    assert n[-1] == approx(N_TH_PROFESSOR, rel=1e-9)
    assert s[-1] == approx(S_30MA, rel=1e-9)
    assert abs(phi[-1]) < 1e-9


def test_integrate_degrau_converge_para_o_estacionario():
    n0, s0 = steady_state(30e-3, PROFESSOR)
    n, s, _ = integrate(np.full(8000, 31e-3), 1e-12, PROFESSOR, n0, s0)
    assert s[-1] == approx(S_31MA, rel=1e-6)
    assert n[-1] == approx(N_TH_PROFESSOR, rel=1e-6)


def test_integrate_oscilacao_de_relaxacao():
    n0, s0 = steady_state(30e-3, PROFESSOR)
    dt = 1e-12
    _, s, _ = integrate(np.full(8000, 31e-3), dt, PROFESSOR, n0, s0)
    deviation = s - S_31MA
    time = np.arange(deviation.size) * dt
    crossings = np.where(np.diff(np.sign(deviation)) != 0)[0][:9]
    crossing_times = time[crossings] - deviation[crossings] * dt / (
        deviation[crossings + 1] - deviation[crossings]
    )
    omega = 2 * np.pi / (2 * np.mean(np.diff(crossing_times)))
    assert omega == approx(39003400152.31235, rel=0.01)

    peaks = [
        k
        for k in range(1, deviation.size - 1)
        if deviation[k] > deviation[k - 1] and deviation[k] > deviation[k + 1]
    ][:2]
    gamma = np.log(deviation[peaks[0]] / deviation[peaks[1]]) / ((peaks[1] - peaks[0]) * dt)
    assert gamma == approx(2447390140.743383, rel=0.05)


def test_integrate_subpassos_equivalem_a_grade_fina():
    n0, s0 = steady_state(20e-3, PROFESSOR)
    coarse_current = np.linspace(20e-3, 30e-3, 501)
    fine_current = np.linspace(20e-3, 30e-3, 5001)
    coarse = integrate(coarse_current, 1e-12, PROFESSOR, n0, s0, max_step=0.1e-12)
    fine = integrate(fine_current, 0.1e-12, PROFESSOR, n0, s0, max_step=0.1e-12)
    assert np.max(np.abs(coarse[1] / fine[1][::10] - 1)) < 1e-9
    assert np.max(np.abs(coarse[0] / fine[0][::10] - 1)) < 1e-9
    assert np.max(np.abs(coarse[2] - fine[2][::10])) < 1e-9


NOISY = RateEquationParams(
    gain_coefficient=1e-12,
    transparency_density=1e24,
    photon_lifetime=2e-12,
    carrier_lifetime=2e-9,
    volume=1e-16,
    spontaneous_emission_rate=5e27,
)
NOISE_N0 = 1.5e24
NOISE_S0 = 2e21
NOISE_DT = 1e-13
NOISE_REALIZATIONS = 100_000


def one_noise_step(seed=20):
    n0 = np.full(NOISE_REALIZATIONS, NOISE_N0)
    s0 = np.full(NOISE_REALIZATIONS, NOISE_S0)
    result = integrate(
        np.array([30e-3, 30e-3]),
        NOISE_DT,
        NOISY,
        n0,
        s0,
        max_step=NOISE_DT,
        rng=np.random.default_rng(seed),
    )
    dn, ds, dphi = derivatives(NOISE_N0, NOISE_S0, 30e-3, NOISY)
    deterministic = (
        NOISE_N0 + dn * NOISE_DT,
        NOISE_S0 + ds * NOISE_DT,
        dphi * NOISE_DT,
    )
    return result, deterministic


def test_integrate_variancia_do_ruido_de_fotons():
    (_, s, _), (_, deterministic_s, _) = one_noise_step()
    increments = s[:, 1] - deterministic_s
    expected_variance = 2 * NOISY.spontaneous_emission_rate * NOISE_S0 * NOISE_DT
    assert np.var(increments) == approx(expected_variance, rel=0.025)


def test_integrate_variancia_do_ruido_de_fase():
    (_, _, phi), (_, _, deterministic_phi) = one_noise_step()
    increments = phi[:, 1] - deterministic_phi
    expected_variance = NOISY.spontaneous_emission_rate * NOISE_DT / (2 * NOISE_S0)
    assert np.var(increments) == approx(expected_variance, rel=0.025)


def test_integrate_forcas_de_langevin_independentes():
    (_, s, phi), (_, deterministic_s, deterministic_phi) = one_noise_step()
    correlation = np.corrcoef(s[:, 1] - deterministic_s, phi[:, 1] - deterministic_phi)[
        0, 1
    ]
    assert abs(correlation) < 0.02


def test_integrate_ruido_reproduzido_pela_seed():
    current = np.full(20, 30e-3)
    args = (current, 0.1e-12, NOISY, NOISE_N0, NOISE_S0)
    a = integrate(*args, rng=np.random.default_rng(31))
    b = integrate(*args, rng=np.random.default_rng(31))
    c = integrate(*args, rng=np.random.default_rng(32))
    assert all(np.array_equal(x, y) for x, y in zip(a, b))
    assert not np.array_equal(a[1], c[1])


def test_integrate_rsp_zero_desliga_ruido():
    without_spontaneous_emission = dataclasses.replace(NOISY, spontaneous_emission_rate=0.0)
    current = np.full(20, 30e-3)
    deterministic = integrate(
        current, 0.1e-12, without_spontaneous_emission, NOISE_N0, NOISE_S0
    )
    with_rng = integrate(
        current,
        0.1e-12,
        without_spontaneous_emission,
        NOISE_N0,
        NOISE_S0,
        rng=np.random.default_rng(33),
    )
    assert all(np.array_equal(x, y) for x, y in zip(deterministic, with_rng))


def test_integrate_nao_deixa_densidade_de_fotons_negativa():
    clipping_case = RateEquationParams(
        gain_coefficient=0.0,
        transparency_density=0.0,
        photon_lifetime=1e9,
        carrier_lifetime=1e9,
        volume=1.0,
        spontaneous_emission_rate=1e30,
    )
    _, s, _ = integrate(
        np.array([0.0, 0.0]),
        1e-10,
        clipping_case,
        np.zeros(20_000),
        np.full(20_000, 1e20),
        max_step=1e-10,
        rng=np.random.default_rng(34),
    )
    assert np.min(s[:, 1]) == 0.0
    assert np.all(s >= 0.0)


def test_integrate_conjunto_produz_realizacoes_independentes():
    n_realizations = 32
    _, s, phi = integrate(
        np.full(40, 30e-3),
        NOISE_DT,
        NOISY,
        np.full(n_realizations, NOISE_N0),
        np.full(n_realizations, NOISE_S0),
        max_step=NOISE_DT,
        rng=np.random.default_rng(35),
    )
    assert s.shape == (n_realizations, 40)
    assert phi.shape == (n_realizations, 40)
    assert np.unique(s[:, -1]).size == n_realizations
    assert np.unique(phi[:, -1]).size == n_realizations


def test_integrate_largura_de_linha_sem_alpha_h():
    # Para alpha_H = 0, Eq. 3.3.37: Delta_nu = R_sp/(4*pi*S). Um passo basta
    # para medir a taxa de difusao: Delta_nu = Var(dphi)/(2*pi*dt).
    (_, _, phi), (_, _, deterministic_phi) = one_noise_step(seed=36)
    phase_increment = phi[:, 1] - deterministic_phi
    measured_linewidth = np.var(phase_increment) / (2 * np.pi * NOISE_DT)
    expected_linewidth = NOISY.spontaneous_emission_rate / (4 * np.pi * NOISE_S0)
    assert measured_linewidth == approx(expected_linewidth, rel=0.025)
