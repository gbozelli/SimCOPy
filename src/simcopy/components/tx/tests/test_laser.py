"""Integracao dos kernels de laser com Waveform e SimContext."""

import numpy as np
import pytest

from simcopy.components.tx.laser import (
    CWLaser,
    DirectModulatedLaser,
    rate_equation_params,
)
from simcopy.config.resolve import resolve_rate_equation_laser
from simcopy.config.system import CWLaserConfig, RateEquationLaserConfig
from simcopy.core.context import SimContext
from simcopy.core.grid import TimeGrid
from simcopy.core.signals import Domain, SignalError, SignalGeometry, Waveform
from simcopy.core.units import C_LIGHT
from simcopy.kernels.laser.conversions import facet_power
from simcopy.kernels.laser.rate_equations import steady_state
from simcopy.testing import approx

FREQUENCY = C_LIGHT / 1550e-9
GRID = TimeGrid(sampling_rate=1e12, n_samples=64)
GEOMETRY = SignalGeometry(n_symbols=16, bits_per_symbol=1, samples_per_symbol=4, n_pol=1)
RATE_CONFIG = RateEquationLaserConfig(
    cavity_length=300e-6,
    cavity_width=3e-6,
    cavity_height=0.2e-6,
    group_index=3.4,
    mirror_loss=2410.0,
    internal_loss=1510.0,
    differential_gain=4e-20,
    transparency_density=1.25e24,
    confinement=1.0,
    carrier_lifetime=2.1e-9,
    spontaneous_emission_factor=1.5,
    gain_compression=0.0,
    linewidth_enhancement=3.0,
)


def waveform(values=None, *, domain=Domain.ELECTRICAL, channel_id=7):
    if values is None:
        values = np.zeros(GRID.n_samples)
    samples = np.asarray(values, dtype=complex).reshape(1, GRID.n_samples)
    return Waveform(
        samples=samples,
        grid=GRID,
        domain=domain,
        geometry=GEOMETRY,
        channel_id=channel_id,
    )


def test_rate_equation_params_mapeia_configuracao_resolvida():
    resolved = resolve_rate_equation_laser(RATE_CONFIG, FREQUENCY)
    params = rate_equation_params(RATE_CONFIG, resolved)
    assert params.gain_coefficient == resolved.gain_coefficient
    assert params.transparency_density == RATE_CONFIG.transparency_density
    assert params.photon_lifetime == resolved.photon_lifetime
    assert params.carrier_lifetime == RATE_CONFIG.carrier_lifetime
    assert params.volume == resolved.volume
    assert params.spontaneous_emission_rate == resolved.spontaneous_emission_rate
    assert params.gain_compression == RATE_CONFIG.gain_compression
    assert params.linewidth_enhancement == RATE_CONFIG.linewidth_enhancement


def test_cw_laser_devolve_waveform_optico_alinhado_em_x():
    source = CWLaser(CWLaserConfig(power=1e-3, linewidth=500e3), FREQUENCY)
    output = source(waveform(), SimContext(GRID, root_seed=10))
    assert output.samples.shape == (2, GRID.n_samples)
    assert np.all(output.samples[1] == 0.0)
    assert output.power == approx(1e-3, rel=1e-12)
    assert output.domain is Domain.OPTICAL
    assert output.grid is GRID
    assert output.geometry is GEOMETRY
    assert output.f0 == FREQUENCY
    assert output.channel_id == 7


def test_cw_laser_reproduz_o_stream_do_contexto():
    source = CWLaser(CWLaserConfig(power=1e-3, linewidth=500e3), FREQUENCY)
    x = waveform()
    first = source(x, SimContext(GRID, root_seed=11))
    repeated = source(x, SimContext(GRID, root_seed=11))
    different = source(x, SimContext(GRID, root_seed=12))
    assert np.array_equal(first.samples, repeated.samples)
    assert not np.array_equal(first.samples, different.samples)


def test_cw_laser_separa_streams_por_canal_e_papel():
    config = CWLaserConfig(power=1e-3, linewidth=500e3)
    context = SimContext(GRID, root_seed=13)
    channel_7 = CWLaser(config, FREQUENCY)(waveform(channel_id=7), context)
    channel_8 = CWLaser(config, FREQUENCY)(waveform(channel_id=8), context)
    local_oscillator = CWLaser(config, FREQUENCY, rng_role="rx_lo_phase")(
        waveform(channel_id=7), context
    )
    assert not np.array_equal(channel_7.samples, channel_8.samples)
    assert not np.array_equal(channel_7.samples, local_oscillator.samples)


def test_cw_laser_exige_channel_id_para_sortear_ruido():
    source = CWLaser(CWLaserConfig(), FREQUENCY)
    with pytest.raises(SignalError, match="channel_id"):
        source(waveform(channel_id=None), SimContext(GRID, root_seed=14))


def test_direct_modulated_laser_converte_sinal_em_corrente():
    pattern = np.resize(np.array([-1.0, 0.0, 1.0]), GRID.n_samples)
    source = DirectModulatedLaser(
        RATE_CONFIG,
        FREQUENCY,
        bias_ratio=2.0,
        modulation_ratio=0.5,
        noise=False,
    )
    threshold = resolve_rate_equation_laser(RATE_CONFIG, FREQUENCY).threshold_current
    assert source.current(waveform(pattern)) == approx(
        threshold * (2.0 + 0.5 * pattern), rel=1e-12
    )


def test_direct_modulated_laser_exige_sinal_eletrico():
    source = DirectModulatedLaser(
        RATE_CONFIG, FREQUENCY, bias_ratio=2.0, modulation_ratio=0.0, noise=False
    )
    with pytest.raises(SignalError, match="espera entrada electrical"):
        source(
            waveform(domain=Domain.OPTICAL),
            SimContext(GRID, root_seed=15),
        )


def test_direct_modulated_laser_deterministico_em_corrente_constante():
    source = DirectModulatedLaser(
        RATE_CONFIG, FREQUENCY, bias_ratio=2.0, modulation_ratio=0.0, noise=False
    )
    x = waveform()
    output = source(x, SimContext(GRID, root_seed=16))

    resolved = resolve_rate_equation_laser(RATE_CONFIG, FREQUENCY)
    params = rate_equation_params(RATE_CONFIG, resolved)
    _, photon_density = steady_state(2 * resolved.threshold_current, params)
    expected_power = facet_power(
        photon_density,
        resolved.volume,
        resolved.group_velocity,
        RATE_CONFIG.mirror_loss,
        resolved.photon_energy,
    )
    assert output.power == approx(expected_power, rel=1e-10)
    assert np.max(np.abs(np.abs(output.samples[0]) ** 2 - expected_power)) < 1e-14
    assert np.all(output.samples[1] == 0.0)
    assert output.domain is Domain.OPTICAL
    assert output.f0 == FREQUENCY
    assert output.channel_id == x.channel_id


def test_direct_modulated_laser_ruido_e_reproduzivel():
    source = DirectModulatedLaser(
        RATE_CONFIG,
        FREQUENCY,
        bias_ratio=2.0,
        modulation_ratio=0.2,
        noise=True,
    )
    values = 0.5 * np.sin(2 * np.pi * np.arange(GRID.n_samples) / 16)
    x = waveform(values)
    first = source(x, SimContext(GRID, root_seed=17))
    repeated = source(x, SimContext(GRID, root_seed=17))
    different = source(x, SimContext(GRID, root_seed=18))
    assert np.array_equal(first.samples, repeated.samples)
    assert not np.array_equal(first.samples, different.samples)
    assert np.all(first.samples[1] == 0.0)


def test_direct_modulated_laser_so_exige_canal_quando_ha_ruido():
    x = waveform(channel_id=None)
    context = SimContext(GRID, root_seed=19)
    noisy = DirectModulatedLaser(
        RATE_CONFIG, FREQUENCY, bias_ratio=2.0, modulation_ratio=0.0, noise=True
    )
    deterministic = DirectModulatedLaser(
        RATE_CONFIG, FREQUENCY, bias_ratio=2.0, modulation_ratio=0.0, noise=False
    )
    with pytest.raises(SignalError, match="channel_id"):
        noisy(x, context)
    assert deterministic(x, context).channel_id is None
