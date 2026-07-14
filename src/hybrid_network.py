import numpy as np

def HybridNetwork(A_sig,A_LO):
    A_1 = 1/2*(A_sig+A_LO)
    A_2 = 1/2*(A_sig-A_LO)
    A_3 = 1/2*(A_sig+1j*A_LO)
    A_4 = 1/2*(A_sig-1j*A_LO)
    return [A_1,A_2,A_3,A_4]

