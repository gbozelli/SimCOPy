"""Grandezas derivadas: o que deixa de ser configuravel e passa a ser calculado."""
import pytest
from tests._approx import approx
from simcopy.config.io import load
from simcopy.config.resolve import resolve

@pytest.fixture(scope="module")
def ref():
    s, sim = load("examples/400zr_4ch_75ghz.yaml")
    return s, sim, resolve(s, sim)

def test_taxa_de_amostragem_unica(ref):
    s, sim, r = ref
    assert sim.sampling_rate == approx(59.84e9 * 16)        # 957,44 GHz

def test_sps_e_derivado_nao_configurado(ref):
    s, sim, r = ref
    assert all(c.samples_per_symbol == 16 for c in r.channels)

def test_janela_bate_com_o_vpi(ref):
    """TimeWindow = NumberOfSymbols / SymbolRate = 2^18 / 59,84 GBd."""
    s, sim, r = ref
    assert r.window == approx(2**18 / 59.84e9, rel=1e-12)
    assert r.n_samples == 2**22                                    # potencia de 2

def test_gamma_bate_com_o_vpi(ref):
    s, sim, r = ref
    assert r.fibers[0].gamma * 1e3 == approx(1.3, rel=0.02)   # 1/W/km

def test_beta2_negativo_para_D_positivo(ref):
    s, sim, r = ref
    assert r.fibers[0].beta2 < 0


# ------------------------------------------------------ contrato: lambda = c/f
# Estes testes existem para pegar deriva na integracao. Se algum ponto do
# resolve voltar a usar lambda fixo (por exemplo `... else 1550e-9`), a 191 THz
# o erro em beta2 e de ~2,6 %, grande demais para passar.
from simcopy.config.system import (SystemConfig, LinkConfig, SpanConfig,
                                   FiberConfig, CoherentChannel)
from simcopy.config.simulation import SimulationConfig
from simcopy.core.units import C_LIGHT, beta2_from_dispersion, gamma_from_n2

@pytest.mark.parametrize("f_THz", [191.0, 193.7, 196.0])
def test_resolve_usa_lambda_derivado_da_frequencia(f_THz):
    s = SystemConfig(channels=(CoherentChannel(index=0, symbol_rate=59.84e9, seed=0),),
                     link=LinkConfig(spans=(SpanConfig(FiberConfig(80e3)),)),
                     center_frequency=f_THz * 1e12)
    r = resolve(s, SimulationConfig.from_symbols(59.84e9, 16, 1024))
    lam = C_LIGHT / (f_THz * 1e12)
    assert r.fibers[0].beta2 == approx(beta2_from_dispersion(16e-6, lam), rel=1e-12)
    assert r.fibers[0].gamma == approx(gamma_from_n2(2.6e-20, 80e-12, lam), rel=1e-12)
