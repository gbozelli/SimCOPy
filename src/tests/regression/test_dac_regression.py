import sys
import os
import numpy as np
import pytest

# Caminhos
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, src_path)

from dac_nyquist import dac_nyquist
from deprecated.DAC_Nyquist import DAC_Nyquist as dac_old

def test_dac_regression():
    np.random.seed(42)
    symbols = np.random.choice([1+1j, -1+1j, 1-1j, -1-1j], size=50)
    
    # Compara nova e antiga
    res_new = dac_nyquist(symbols, samples_per_symbol=16, roll_off=0.2)
    res_old = dac_old(symbols, SpS=16, RollOff=0.2)
    
    # Tolerância necessária para operações de FFT
    np.testing.assert_allclose(res_new, res_old, atol=1e-10)