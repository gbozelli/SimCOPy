"""Leitura de YAML: tx e rx dentro de cada canal, ancoras, coerencia f/lambda."""
import textwrap
import pytest
from tests._approx import approx
from simcopy.config.io import load
from simcopy.config.system import ConfigError

BASE = """
defaults:
  tx: &tx {laser: {power: 1.0e-3, linewidth: 500.0e3}, rolloff: 0.18}
system:
  {centro}
  channels:
    - {index: 0, kind: dcs, modulation: 16QAM, pol: DP, symbol_rate: 59.84e9,
       seed: 0, tx: *tx}
    - index: 1
      kind: dcs
      modulation: 16QAM
      pol: DP
      symbol_rate: 59.84e9
      seed: 1
      tx:
        <<: *tx
        laser: {power: 2.0e-3, linewidth: 100.0e3}
  link:
    spans: [{fiber: {length: 80.0e3}}]
simulation: {samples_per_symbol: 16, n_symbols: 1024}
"""

def _load(tmp_path, centro="center_frequency: 193.7e12"):
    p = tmp_path / "c.yaml"
    p.write_text(BASE.replace("{centro}", centro))
    return load(p)

def test_transmissor_fica_dentro_do_canal(tmp_path):
    s, _ = _load(tmp_path)
    assert s.channels[0].tx.laser.linewidth == 500e3

def test_sobrescrita_local_de_um_canal(tmp_path):
    s, _ = _load(tmp_path)
    assert s.channels[1].tx.laser.linewidth == 100e3       # sobrescrito
    assert s.channels[1].tx.rolloff == 0.18                # herdado da ancora

def test_aceita_comprimento_de_onda_no_lugar_da_frequencia(tmp_path):
    s, _ = _load(tmp_path, "center_wavelength: 1550.0e-9")
    assert s.center_wavelength == approx(1550e-9)

def test_recusa_frequencia_e_lambda_incoerentes(tmp_path):
    with pytest.raises(ConfigError, match="nao sao coerentes"):
        _load(tmp_path, "center_frequency: 193.7e12\n  center_wavelength: 1550.0e-9")
