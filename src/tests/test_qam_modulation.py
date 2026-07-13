import sys
import os
import numpy as np
import pytest

from qam_modulation import qam_modulation

@pytest.mark.parametrize("s, M, expected", [
    (np.array([0, 1, 0, 0, 1, 1, 1, 0]), 4, np.array([-1+1j, 1+1j, -1-1j, 1-1j])),
    (np.array([0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0]), 16, np.array([-3+3j, 1+1j, -3-3j])),
    (np.array([0, 0, 0, 0, 1, 1]), 8, np.array([-3+1j, -1-1j]))
])
def test_qam_modulation_values(s, M, expected):
    result = qam_modulation(s, M)
    np.testing.assert_array_equal(result, expected)

def test_qam_modulation_invalid_m():
    s = np.array([0, 1])
    
    # Testa se a função lança um erro caso um M inválido (ex: 3) seja passado
    with pytest.raises(KeyError):
        qam_modulation(s, 3)