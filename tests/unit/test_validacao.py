"""Invariantes checados antes de alocar memoria."""
import dataclasses
from simcopy.config.io import load
from simcopy.config.simulation import SimulationConfig
from simcopy.validation.rules import validate
from simcopy.validation.base import Severity

S, SIM = load("examples/400zr_4ch_75ghz.yaml")

def _erros(sim):
    return [v for v in validate(S, sim) if v.severity is Severity.ERROR]

def test_cenario_de_referencia_e_valido():
    assert _erros(SIM) == []

def test_nyquist_violado_com_sps_4():
    sim = SimulationConfig.from_symbols(59.84e9, 4, 2**14)
    assert any(v.validator == "RequiresNyquistSatisfied" for v in _erros(sim))

def test_sps_nao_inteiro_e_recusado():
    sim = dataclasses.replace(SIM, sampling_rate=900e9)
    assert any(v.validator == "IntegerSamplesPerSymbol" for v in _erros(sim))

def test_sps_do_dsp_tem_que_dividir_o_da_grade():
    """400ZR: grade a 16, DSP a 2 -> dizimacao por 8. DSP a 3 nao divide 16."""
    import dataclasses as dc
    ch = S.channels[0]
    rx3 = dc.replace(ch.rx, adc=dc.replace(ch.rx.adc, samples_per_symbol=3))
    s3 = dc.replace(S, channels=(dc.replace(ch, rx=rx3),) + S.channels[1:])
    erros = [v for v in validate(s3, SIM) if v.severity is Severity.ERROR]
    assert any(v.validator == "DSPRateDividesSimulationRate" for v in erros)
