import matplotlib.pyplot as plt
import numpy as np
# DAC_Nyquist: DAC considering Nyquist filtering
def DAC_Nyquist(s=None,SpS=16,RollOff=0.2,plot_flag=False):

    # Initial calculations
    N_sym = len(s)
    N_samples = N_sym*SpS

    # Conversion to the frequency domain
    S = np.fft.fftshift(np.fft.fft(s))
    S_os = np.zeros(N_samples,dtype=complex)

    # Ideal interpolation (Nyquist limit)
    ind_i = int(N_samples/2-N_sym/2)
    ind_f = ind_i+N_sym
    S_os[ind_i:ind_f] = S
    s_os = SpS*np.fft.ifft(np.fft.ifftshift(S_os))
    S_nyq = np.zeros(N_samples,dtype=complex)
    ind_i = int(N_samples/2-N_sym-N_sym/2)
    ind_f = ind_i+N_sym
    S_nyq[ind_i:ind_f] = S
    ind_i = int(N_samples/2-N_sym/2)
    ind_f = ind_i+N_sym
    S_nyq[ind_i:ind_f] = S
    ind_i = int(N_samples/2+N_sym-N_sym/2)
    ind_f = ind_i+N_sym
    S_nyq[ind_i:ind_f] = S

    ind = np.arange(N_samples)
    N_aux = int(RollOff*N_sym)
    ind_aux = np.arange(int(N_aux))
    s_aux_1 = 0.5*(1-np.cos(np.pi*ind_aux/N_aux))
    s_aux_2 = 0.5*(1+np.cos(np.pi*ind_aux/N_aux))
    window = np.zeros(N_samples)
    ind_i_1 = int(N_samples/2-N_sym/2-N_aux/2)
    ind_f_1 = ind_i_1+int(N_aux)
    window[ind_i_1:ind_f_1] = s_aux_1
    ind_i_2 = int(N_samples/2+N_sym/2-N_aux/2)
    ind_f_2 = ind_i_2+int(N_aux)
    window[ind_i_2:ind_f_2] = s_aux_2
    window[ind_f_1:ind_i_2] = np.ones(ind_i_2-ind_f_1)

    if plot_flag:
        plt.figure()
        plt.plot(np.abs(S_nyq)/max(np.abs(S_nyq)))
        plt.plot(np.abs(window)/max(window))

    S_nyq = window*S_nyq
    s_nyq = SpS*np.fft.ifft(np.fft.ifftshift(S_nyq))

    return s_nyq
