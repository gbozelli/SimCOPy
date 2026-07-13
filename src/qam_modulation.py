import numpy as np

def qam_modulation(raw_bits: np.ndarray, modulation_order: int) -> np.ndarray:
    """
    Modulates a binary sequence into QAM symbols.
    """

    N_bits = int(np.log2(modulation_order))
    N_sym = len(raw_bits) // N_bits
    s_1 = np.reshape(raw_bits, (N_sym, N_bits))
    
    weights = 1 << np.arange(N_bits)[::-1]
    indices = s_1.dot(weights)
    
    mapping = {
        2: np.array([1j, -1j]),
        
        4: np.array([
            1+1j, -1+1j, 1-1j, -1-1j
        ]),
        
        8: np.array([
            -3+1j, -3-1j, -1+1j, -1-1j, 
            3+1j, 3-1j, 1+1j, 1-1j
        ]),
        
        16: np.array([
            -3-3j, -3-1j, -3+3j, -3+1j, 
            -1-3j, -1-1j, -1+3j, -1+1j, 
            3-3j, 3-1j, 3+3j, 3+1j, 
            1-3j, 1-1j, 1+3j, 1+1j
        ]),
        
        32: np.array([
            -1-1j, -3-1j, -1-3j, -3-3j, 
            -1+1j, -3+1j, -1+3j, -3+3j, 
            1-1j, 3-1j, 1-3j, 3-3j, 
            1+1j, 3+1j, 1+3j, 3+3j, 
            -3-5j, -5-1j, -1-5j, -5-3j, 
            -3+5j, -5+1j, -1+5j, -5+3j, 
            3-5j, 5-1j, 1-5j, 5-3j, 
            3+5j, 5+1j, 1+5j, 5+3j
        ])
    }
    
    return mapping[modulation_order][indices]