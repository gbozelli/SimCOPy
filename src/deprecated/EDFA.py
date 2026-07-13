import numpy as np

def EDFA(s_in,t,G_edfa_dB,NF_dB,lambda0=1550e-9):
    h = 6.626e-34     # Planck's constant in J/Hz
    c = 3e8           # Speed of light in the vacuum in m/s
    f0 = c/lambda0    # Operation frequency in Hz
    ts = t[1]-t[0]    # Sampling time in s
    fs = 1/ts         # Sampling frequency in Hz
    BW_sim = fs       # Simulation bandwidth in Hz

    # Gain
    G_edfa = 10**(G_edfa_dB/10)
    sx = s_in[0]
    sy = s_in[1]
    N = len(sx)

    sx = np.sqrt(G_edfa)*sx
    sy = np.sqrt(G_edfa)*sy

    # Noise
    F_n = 10**(NF_dB/10)
    n_sp = F_n/2
    S_ASE = n_sp*h*f0*(G_edfa-1)
    P_ASE = S_ASE*BW_sim

    s_ASEx = np.sqrt(P_ASE)/2*np.random.rand(N)+1j*np.sqrt(P_ASE)/2*np.random.rand(N)
    s_ASEy = np.sqrt(P_ASE)/2*np.random.rand(N)+1j*np.sqrt(P_ASE)/2*np.random.rand(N)
    s_outx = sx + s_ASEx
    s_outy = sy + s_ASEy

    s = np.zeros((2,len(t)),dtype = np.complex128)
    s[0,:] = s_outx
    s[1,:] = s_outy
    return s
