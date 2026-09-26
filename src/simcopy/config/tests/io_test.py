"""Arquivo rascunho para alterar a estrutura de testes"""

from simcopy.config.io import _center_frequency, _num, _laser, _mzm, _pd, _adc, _coherent, _imdd, load

from simcopy.config.simulation import SimulationConfig
from simcopy.config.system import CoherentChannel, IMDDChannel, LaserConfig, SystemConfig
from simcopy.core.units import C_LIGHT
_INT_KEYS = {"index", "seed", "s21_order", "filter_order", "samples_per_symbol",
             "resolution_bits", "n_symbols", "root_seed"}

def test_num():
  assert (
    _num([59.84e9], list), 
    _num({59.84e9}, dict), 
    _num("59.84e9", str)
    ) == 59.84e9

def test_coherent():
  channel = {
    "index": 1,
    "symbol_rate": 32e9,
    "seed": 1234,
    "modulation": "DP-QPSK",
    "tx": {
      "laser": {"power": 0.0, "linewidth": 100e3},
      "mzm": {"vpi": 3.14, "vbias": 1.57},
      "iq": {"gain_imbalance": 0.1, "phase_imbalance": 5.0},
      "rolloff": 0.2
    },
    "rx": {
      "lo": {"power": 0.0, "linewidth": 100e3},
      "photodiode": {"responsivity": 1.0, "dark_current": 1e-9},
      "adc": {"resolution_bits": 8, "sampling_rate": 64e9},
      "optical_bandwidth": 70e9
    }
  }
  c = _coherent(channel)
  assert isinstance(c, CoherentChannel)
  assert c.index == 1
  assert c.tx.laser.power == 0.0
  assert c.rx.photodiode.responsivity == 1.0
  assert c.rx.adc.resolution_bits == 8
  assert c.rx.optical_bandwidth == 70e9
  assert c.tx.rolloff == 0.2
  assert c.tx.iq.gain_imbalance == 0.1
  assert c.tx.iq.phase_imbalance == 5.0
  assert c.tx.mzm.vpi == 3.14
  assert c.tx.mzm.vbias == 1.57
  assert c.rx.lo.power == 0.0
  assert c.rx.lo.linewidth == 100e3
  assert c.rx.lo.linewidth == 100e3
  # precisamos escrever outrps testes, mas vc pegou a ideia

def test_imdd():
  channel = {
    "index": 2,
    "symbol_rate": 10e9,
    "seed": 5678,
    "modulation": "OOK",
    "tx": {
      "laser": {"power": 0.0, "linewidth": 100e3},
      "mzm": {"vpi": 3.14, "vbias": 1.57}
    },
    "rx": {
      "photodiode": {"responsivity": 1.0, "dark_current": 1e-9},
      "adc": {"resolution_bits": 8, "sampling_rate": 20e9}
    }
  }
  c = _imdd(channel)
  assert isinstance(c, IMDDChannel)
  assert c.index == 2
  assert c.tx.laser.power == 0.0
  assert c.rx.photodiode.responsivity == 1.0
  assert c.rx.adc.resolution_bits == 8
  assert c.tx.mzm.vpi == 3.14
  assert c.tx.mzm.vbias == 1.57
    

def test_center_frequency():
    freq, lam = C_LIGHT / 1550e-9, 1550e-9
    assert _center_frequency({"center_frequency": freq, "center_wavelength": lam}) == freq
    assert _center_frequency({"center_wavelength": lam}) == freq
    assert _center_frequency({"center_frequency": freq}) == freq

def test_load():
    path = "tests_files_yaml/config_test.yaml"
    system, simulation = load(path)
    assert isinstance(system, SystemConfig)
    assert isinstance(simulation, SimulationConfig)