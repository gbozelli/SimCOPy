"""Cascata de configuracao: o que cada decisao libera e o que bloqueia."""
import pytest
from tests._approx import approx
from simcopy.config.system import (IMDDChannel, CoherentChannel, SystemConfig,
    LinkConfig, SpanConfig, FiberConfig, Pol, QAM16, PAM4, ConfigError)
from simcopy.core.units import C_LIGHT

def _sys(*chans):
    return SystemConfig(channels=chans,
                        link=LinkConfig(spans=(SpanConfig(FiberConfig(125e3)),)))

def test_imdd_bloqueia_qam():
    with pytest.raises(ConfigError, match="kind = IMDD"):
        IMDDChannel(index=0, symbol_rate=50e9, seed=0, modulation=QAM16)

def test_dcs_bloqueia_pam4():
    with pytest.raises(ConfigError, match="kind = DCS"):
        CoherentChannel(index=0, symbol_rate=50e9, seed=0, modulation=PAM4)

def test_imdd_e_sempre_sp():
    assert IMDDChannel(index=0, symbol_rate=50e9, seed=0).pol is Pol.SP

def test_enlace_vetorial_se_houver_um_canal_dp():
    """[D-14]: SP e DP coexistem; basta um DP para o solver ser vetorial."""
    sp = CoherentChannel(index=0, symbol_rate=50e9, seed=0, pol=Pol.SP)
    dp = CoherentChannel(index=1, symbol_rate=50e9, seed=1, pol=Pol.DP)
    assert not _sys(sp).link_is_vectorial
    assert _sys(sp, dp).link_is_vectorial

def test_lambda_e_derivado_da_frequencia():
    """lambda e f tem um unico dono: lambda = c/f. Nao existe 1550e-9 fixo."""
    s = _sys(CoherentChannel(index=0, symbol_rate=50e9, seed=0))
    assert s.center_wavelength * s.center_frequency == approx(C_LIGHT)
    assert s.center_wavelength == approx(1547.715e-9, rel=1e-6)

def test_grade_wdm_simetrica_e_espacada():
    s = SystemConfig(channels=tuple(CoherentChannel(index=i, symbol_rate=59.84e9, seed=i)
                                    for i in range(4)),
                     link=LinkConfig(spans=(SpanConfig(FiberConfig(125e3)),)),
                     channel_spacing=75e9)
    f = s.channel_offsets()
    assert f == approx([-112.5e9, -37.5e9, 37.5e9, 112.5e9])
