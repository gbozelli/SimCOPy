import numpy as np
import matplotlib.pyplot as plt

def spectrum_analysis(s,t):
    ts = t[1]-t[0]                       # Sampling period in s
    fs = 1/ts                            # Sampling frequancy in Hz
    Deltaf = 1/(ts*len(t))               # Frequency resolution in Hz
    f = np.arange(len(t))*Deltaf-fs/2    # Frequency vector in Hz

    plt.figure()
    plt.plot(f/1e9,30+10*np.log10(np.fft.fftshift(np.abs(np.fft.fft(s)))))
    plt.xlabel('Frequency [GHz]')
    plt.ylabel('Power [dBm]')

"""### 6.1.4 Optical transmitter and receiver"""

# DP-QAM transmitter
