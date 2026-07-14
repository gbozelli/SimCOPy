import numpy as np

def single_frequency_laser(t = None, P = 0.01, Delta_nu = 1e6, Freq_offset = 0, X_fraction = 1, phi_pol = 0):
    # Calculate the amplitude
    A = np.sqrt(P)
    # Calculate the phase noise
    ts = t[1]-t[0]
    freq_noise = np.sqrt(2*np.pi*Delta_nu*ts)*np.random.randn(len(t))
    phase_noise = np.cumsum(freq_noise)
    # Frequency offset
    phase = 2*np.pi*(Freq_offset)*t+phase_noise
    # Generation of the signal
    s = A*np.exp(1j*phase)
    Ax = np.sqrt(X_fraction)*np.exp(1j*phi_pol/2)*s
    Ay = np.sqrt(1-X_fraction)*np.exp(-1j*phi_pol/2)*s
    A = np.zeros((2,len(t)),dtype=complex)
    #print(np.shape(A))
    A[0,:] = Ax
    A[1,:] = Ay
    return(A)

