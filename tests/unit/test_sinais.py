"""Tipos de sinal: as relacoes de dimensao e as recusas de entrada invalida."""
import numpy as np
import pytest
from tests._approx import approx
from simcopy.core.grid import TimeGrid
from simcopy.core.signals import SignalGeometry, Waveform, Domain, SignalError

def test_geometria_calcula_bits_e_amostras():
    g = SignalGeometry(n_symbols=1000, bits_per_symbol=4, samples_per_symbol=16, n_pol=2)
    assert g.n_bits == 1000 * 4 * 2          # simbolos x log2(M) x polarizacoes
    assert g.n_samples == 1000 * 16

def test_geometria_recusa_tres_polarizacoes():
    with pytest.raises(SignalError):
        SignalGeometry(10, 4, 16, n_pol=3)

def _wf(samples, n=128):
    return Waveform(samples, TimeGrid(1e12, n), Domain.OPTICAL,
                    SignalGeometry(8, 4, 16, samples.shape[0]))

def test_waveform_recusa_array_real():
    with pytest.raises(SignalError, match="complexo"):
        _wf(np.ones((2, 128)))

def test_waveform_recusa_numero_de_amostras_errado():
    with pytest.raises(SignalError, match="nao batem com a grade"):
        _wf(np.ones((2, 100), complex))

def test_potencia_soma_as_duas_polarizacoes():
    """Regressao do bug C-5: o codigo antigo somava a polarizacao X duas vezes."""
    s = np.vstack([np.full(128, 1.0 + 0j), np.full(128, 2.0 + 0j)])
    assert _wf(s).power == approx(1.0 + 4.0)

def test_precondicao_de_dominio():
    w = _wf(np.ones((2, 128), complex))
    with pytest.raises(SignalError, match="espera entrada electrical"):
        w.require(Domain.ELECTRICAL, "fotodiodo")


# ------------------------------------------------ quarta representacao: digital
from simcopy.core.signals import DigitalSignal

def _dig(n_sym=100, sps=2, n_pol=2, **kw):
    geo = SignalGeometry(n_sym, 4, sps, n_pol)
    return DigitalSignal(np.zeros((n_pol, n_sym * sps), complex), 59.84e9, geo, **kw)

def test_digital_tem_taxa_propria_fora_da_grade():
    """2 amostras/simbolo a 59,84 GBd = 119,68 GSa/s, nao os 957,44 da grade."""
    d = _dig(sps=2)
    assert d.sampling_rate == approx(119.68e9)

def test_digital_recusa_numero_de_amostras_errado():
    geo = SignalGeometry(100, 4, 2, 2)
    with pytest.raises(SignalError, match="DigitalSignal"):
        DigitalSignal(np.zeros((2, 150), complex), 59.84e9, geo)

def test_digital_recusa_resolucao_invalida():
    with pytest.raises(SignalError, match="resolution_bits"):
        _dig(resolution_bits=0)

def test_digital_sem_quantizacao_e_o_padrao():
    assert _dig().resolution_bits is None

def test_digital_extrai_uma_amostra_por_simbolo():
    geo = SignalGeometry(4, 4, 2, 1)
    s = np.array([[0, 1, 2, 3, 4, 5, 6, 7]], dtype=complex)
    d = DigitalSignal(s, 59.84e9, geo)
    np.testing.assert_array_equal(d.symbol_samples(0), [[0, 2, 4, 6]])
    np.testing.assert_array_equal(d.symbol_samples(1), [[1, 3, 5, 7]])
