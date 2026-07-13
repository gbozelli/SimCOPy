import numpy as np
import matplotlib.pyplot as plt

def dac_nyquist(symbols: np.ndarray, samples_per_symbol: int = 16, roll_off: float = 0.2, plot_flag: bool = False) -> np.ndarray:
    """
    Generates a Nyquist-shaped signal from a sequence of symbols.
    """
    N_sym = len(symbols)
    N_samples = N_sym * samples_per_symbol

    S = np.fft.fftshift(np.fft.fft(symbols))
    
    S_nyq = np.zeros(N_samples, dtype=complex)
    
    ind_center = N_samples // 2
    half_sym = N_sym // 2
    
    S_nyq[ind_center - half_sym - N_sym : ind_center + half_sym - N_sym] = S
    S_nyq[ind_center - half_sym : ind_center + half_sym] = S
    S_nyq[ind_center - half_sym + N_sym : ind_center + half_sym + N_sym] = S

    N_aux = int(roll_off * N_sym)
    ind_aux = np.arange(N_aux)
    
    s_aux_1 = 0.5 * (1 - np.cos(np.pi * ind_aux / N_aux))
    s_aux_2 = 0.5 * (1 + np.cos(np.pi * ind_aux / N_aux))
    
    window = np.zeros(N_samples)
    
    ind_i_1 = ind_center - half_sym - N_aux // 2
    ind_f_1 = ind_i_1 + N_aux
    window[ind_i_1:ind_f_1] = s_aux_1
    
    ind_i_2 = ind_center + half_sym - N_aux // 2
    ind_f_2 = ind_i_2 + N_aux
    window[ind_i_2:ind_f_2] = s_aux_2
    
    window[ind_f_1:ind_i_2] = 1.0

    if plot_flag:
        plt.figure()
        plt.plot(np.abs(S_nyq) / np.max(np.abs(S_nyq)))
        plt.plot(np.abs(window) / np.max(window))
        plt.show()

    S_nyq_windowed = window * S_nyq
    s_nyq = samples_per_symbol * np.fft.ifft(np.fft.ifftshift(S_nyq_windowed))

    return s_nyq