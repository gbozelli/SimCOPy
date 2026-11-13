"""Resposta de pequeno sinal do laser (Agrawal, Eqs. 3.3.17--3.3.23)."""

import dataclasses

import numpy as np

from simcopy.kernels.laser.rate_equations import (
    RateEquationParams,
    derivatives,
    integrate,
    steady_state,
)
from simcopy.kernels.laser.small_signal import (
    bandwidth_3db,
    linearize,
    modulation_response,
    normalized_modulation_response,
)
from simcopy.testing import approx

PARAMS = RateEquationParams(
    gain_coefficient=3.526970094117647e-12,
    transparency_density=1.25e24,
    photon_lifetime=2.8931579685554005e-12,
    carrier_lifetime=2.1e-9,
    volume=1.8e-16,
    spontaneous_emission_rate=3.968253968253968e-16,
    linewidth_enhancement=3.0,
)
CURRENT = 30e-3


def test_linearize_sem_compressao():
    small_signal = linearize(CURRENT, PARAMS)
    assert small_signal.carrier_density == approx(1.348e24, rel=1e-12)
    assert small_signal.photon_density == approx(1.1524800754727476e21, rel=1e-12)
    assert small_signal.a_sn == approx(4064762760.2588296, rel=1e-12)
    assert small_signal.a_ns == approx(345643069223.5289, rel=1e-12)
    assert small_signal.gamma_n == approx(4540953236.449306, rel=1e-12)
    assert small_signal.gamma_r == approx(2270476618.224653, rel=1e-12)
    assert small_signal.omega_r == approx(37413928048.88923, rel=1e-12)


def test_linearize_compressao_usa_jacobiano_completo():
    compressed = dataclasses.replace(PARAMS, gain_compression=1e-22)
    n, s = steady_state(CURRENT, compressed)
    small_signal = linearize(CURRENT, compressed)

    # Jacobiano numerico das equacoes originais, independente das formulas de linearize.
    dn_step, ds_step = n * 1e-4, s * 1e-4
    dn_plus, ds_plus, _ = derivatives(n + dn_step, s, CURRENT, compressed)
    dn_minus, ds_minus, _ = derivatives(n - dn_step, s, CURRENT, compressed)
    dn_s_plus, ds_s_plus, _ = derivatives(n, s + ds_step, CURRENT, compressed)
    dn_s_minus, ds_s_minus, _ = derivatives(n, s - ds_step, CURRENT, compressed)

    assert (ds_plus - ds_minus) / (2 * dn_step) == approx(small_signal.a_sn, rel=1e-9)
    assert (ds_s_plus - ds_s_minus) / (2 * ds_step) == approx(
        -small_signal.gamma_p, rel=1e-9
    )
    assert (dn_s_plus - dn_s_minus) / (2 * ds_step) == approx(-small_signal.a_ns, rel=1e-9)
    assert (dn_plus - dn_minus) / (2 * dn_step) == approx(-small_signal.gamma_n, rel=1e-9)


def test_modulation_response_em_frequencias_de_referencia():
    frequency = np.array([0.0, 1e9, 5e9, 10e9])
    expected = np.array(
        [
            1.0032039841437237e23 + 0j,
            1.03175793e23 - 2.15585624e21j,
            3.02015266e23 - 1.03074921e23j,
            -5.47384276e22 - 6.14175588e21j,
        ]
    )
    assert modulation_response(frequency, CURRENT, PARAMS) == approx(expected, rel=2e-8)


def test_resposta_normalizada_e_banda_de_3db():
    assert normalized_modulation_response(0.0, CURRENT, PARAMS) == 1.0 + 0.0j
    bandwidth = bandwidth_3db(CURRENT, PARAMS)
    assert bandwidth == approx(10313707264.829136, rel=1e-12)
    # O projeto define a banda pela amplitude da modulacao de potencia cair a 1/2.
    assert abs(normalized_modulation_response(bandwidth, CURRENT, PARAMS)) == approx(
        0.5, rel=1e-12
    )


def test_resposta_de_pequeno_sinal_confere_com_integracao_no_tempo():
    frequency = 5e9
    current_amplitude = 10e-6
    dt = 1e-12
    time = np.arange(10_000) * dt
    current = CURRENT + current_amplitude * np.cos(2 * np.pi * frequency * time)
    n0, s0 = steady_state(current[0], PARAMS)
    _, photon_density, _ = integrate(current, dt, PARAMS, n0, s0, max_step=0.1e-12)

    # Lock-in complexo depois de 2 ns de acomodacao. Para entrada cos(wt), o fator 2
    # recupera a amplitude complexa da componente positiva.
    steady = slice(2000, None)
    centered = photon_density[steady] - np.mean(photon_density[steady])
    measured = (
        2
        * np.mean(centered * np.exp(-1j * 2 * np.pi * frequency * time[steady]))
        / current_amplitude
    )
    expected = modulation_response(frequency, CURRENT, PARAMS)
    assert abs(measured) == approx(abs(expected), rel=0.015)
    assert abs(np.angle(measured / expected)) < 0.02
