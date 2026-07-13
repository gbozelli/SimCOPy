import sys
import os
import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from dac_nyquist import dac_nyquist

@pytest.fixture
def sample_symbols():
    np.random.seed(42)
    return np.random.choice([1+1j, -1+1j, 1-1j, -1-1j], size=100)

@pytest.mark.parametrize("samples_per_symbol", [8, 16, 32])
def test_dac_nyquist_output_length(sample_symbols, samples_per_symbol):
    result = dac_nyquist(sample_symbols, samples_per_symbol=samples_per_symbol)
    
    expected_length = len(sample_symbols) * samples_per_symbol
    assert len(result) == expected_length

def test_dac_nyquist_output_type(sample_symbols):
    result = dac_nyquist(sample_symbols)
    
    assert isinstance(result, np.ndarray)
    assert np.iscomplexobj(result)

def test_dac_nyquist_zero_input():
    symbols = np.zeros(50, dtype=complex)
    result = dac_nyquist(symbols, samples_per_symbol=16)
    
    np.testing.assert_allclose(result, np.zeros(50 * 16, dtype=complex))