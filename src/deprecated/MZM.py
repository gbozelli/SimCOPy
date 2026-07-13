import numpy as np

def MZM(A,Vrf_upper,Vrf_lower,Splitter_IL=1.0,Upper_IL=2.0,Lower_IL=3.0,Combiner_IL=1.1,Vpi = 5,Vbias=2.5):
    Aupper_1 = A/2*10**(-Splitter_IL/20)
    Alower_1 = A/2*10**(-Splitter_IL/20)
    Aupper_2 = Aupper_1*np.exp(1j*(-Vrf_upper+Vbias/2)/Vpi*np.pi)
    Alower_2 = Alower_1*np.exp(1j*(-Vrf_lower-Vbias/2)/Vpi*np.pi)
    Aupper_3 = Aupper_2*10**(-Upper_IL/20)
    Alower_3 = Alower_2*10**(-Lower_IL/20)
    Aupper_4 = Aupper_3*10**(-Combiner_IL/20)
    Alower_4 = Alower_3*10**(-Combiner_IL/20)
    Aout = Aupper_4+Alower_4
    return Aout
