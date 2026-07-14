import sys
import os
import numpy as np
import pytest

# Caminhos
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, src_path)

from qam_modulation import qam_modulation
from deprecated.QAM_mod import QAM_mod as qam_old 

@pytest.mark.parametrize("M", [4, 8, 16, 32])
def test_qam_regression(M):
    num_symbols = 100
    N_bits = int(np.log2(M))
    bits = np.random.randint(0, 2, num_symbols * N_bits)
    
    # Compara a saída da nova função com a antiga
    res_new = qam_modulation(bits, M)
    res_old = qam_old(bits, M)
    
    np.testing.assert_array_equal(res_new, res_old)