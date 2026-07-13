import numpy as np
from scipy import signal

def photodiode(s=None,t=None,R=0.9,BW=10e9,shot_noise=True,thermal_noise=True):
    ts = t[2]-t[1]
    fs = 1/ts
    BW_sim = 1/ts
    # Extract each polarization
    N = len(t)
    Ax = s[:]

    #print('The length of Ax is',np.shape(Ax))
    #print('The length of t is',np.shape(t))

    # Photodetection
    #I_opt = (np.abs(Ax))**2
    I_opt = np.multiply(abs(Ax),abs(Ax))
    i_pd = R*I_opt

    # Noise addition
    kB = 1.38e-23       # Boltzmann's constant in J/K
    Temp = 300          # Temperature in K
    Rl = 50             # Load resistance in Ohm
    q = 1.6e-19         # Electron charge in C
    if thermal_noise:
        # Thermal noise
        PSD_n_thermal = 4*kB*Temp/Rl
        P_n_thermal = PSD_n_thermal*BW_sim
        n_thermal = np.sqrt(P_n_thermal)*np.random.randn(1,N)
    else:
        n_thermal = 0

    if shot_noise:
        # Shot noise
        PSD_n_shot = 2*q*np.mean(np.abs(i_pd))
        P_n_shot = PSD_n_shot*BW_sim
        n_shot = np.sqrt(P_n_shot)*np.random.randn(1,N)
    else:
        n_shot = 0
    i_pd = i_pd + n_thermal+n_shot

    # Filtering
    wp = BW
    ws = 1.5*wp      # Begin of the stop band
    gpass = 3        # Maximum loss in the pass band in dB
    gstop = 20       # Minimum loss in the stop band in dB
    analog = False   #

    # Finding the order
    No, Wn = signal.ellipord(wp, ws, gpass, gstop, analog, fs)

    # Create the filter
    b, a = signal.ellip(No, gpass, gstop, Wn, btype='low', analog=False, fs=fs)

    i_pd = signal.lfilter(b,a,i_pd)

    return i_pd
