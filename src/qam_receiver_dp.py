import numpy as np
import matplotlib.pyplot as plt
from qam_dem import qam_dem
from dac_nyquist import dac_nyquist

def QAM_receiver_DP(s,M,N_sync=128,N_inf=4096,SpS=16,RollOff=0.2,sync_seed=0,ts=1e-12,plot_flag=False,s_b=None):
    s_in = s
    s_x = s[0,:]
    s_y = s[1,:]

    # Generation of time vector
    N_samples = len(s)
    t = ts*(np.arange(N_samples))
    fs = 1/ts
    f = (np.arange(N_samples)/N_samples)*fs-fs/2

    # Coarse frequency offset compensation
    N_samples = len(s_x)
    N_sym = int(N_samples/SpS)
    S_nyq = np.zeros(N_samples,dtype=complex)
    ind_i = int(N_samples/2-N_sym/2)
    ind_f = int(ind_i+N_sym)
    window = np.zeros(N_samples)
    window[ind_i:ind_f] = np.ones(ind_f-ind_i)

    f_min, f_max, N_f = -500e6,500e6,101
    f_offset_array = np.linspace(f_min,f_max,N_f)
    P_array_x = np.linspace(f_min,f_max,N_f)
    P_array_y = np.linspace(f_min,f_max,N_f)
    for counter, offset in enumerate(f_offset_array):
        s_shift_x = np.exp(1j*2*np.pi*offset*t)*s_x
        S_shift_x = np.abs(np.fft.fftshift(np.fft.fft(s_shift_x)))
        S_f_x = S_shift_x*window
        P_array_x[counter] = sum(np.abs(S_f_x)**2)
        s_shift_y = np.exp(1j*2*np.pi*offset*t)*s_y
        S_shift_y = np.abs(np.fft.fftshift(np.fft.fft(s_shift_y)))
        S_f_y = S_shift_y*window
        P_array_y[counter] = sum(np.abs(S_f_y)**2)

    max_offset_opt = max(P_array)
    ind_offset_opt = np.where(P_array==max_offset_opt)
    f_offset_opt = f_offset_array[ind_offset_opt[0][0]]

    # Refining
    p = np.polyfit(f_offset_array,P_array,4)
    N_f = 2001
    f_offset_array_ext = np.linspace(f_min,f_max,N_f)
    P_array_ext = np.polyval(p,f_offset_array_ext)
    if plot_flag:
        plt.figure()
        plt.plot(f_offset_array/1e6,P_array,'o')
        plt.plot(f_offset_array_ext/1e6,P_array_ext)
        plt.ylabel('Power [a.u.]')
        plt.xlabel('Frequency offset [MHz]')

    max_offset_opt = max(P_array_ext)
    ind_offset_opt = np.where(P_array_ext==max_offset_opt)
    f_offset_opt = f_offset_array_ext[ind_offset_opt[0][0]]

    #f_offset_opt = 0
    s = np.exp(1j*2*np.pi*f_offset_opt*t)*s

    #s = s_in

    # Filter the signal using a Nyquist filter
    N_samples = len(s)
    N_sym = int(N_samples/SpS)
    S_nyq = np.zeros(N_samples,dtype=complex)
    ind_i = int(N_samples/2-N_sym-N_sym/2)
    ind_f = ind_i+N_sym

    ind = np.arange(N_samples)
    N_aux = int(RollOff*N_sym)
    ind_aux = np.arange(int(N_aux))
    s_aux_1 = 1*np.ones(N_aux)   #0.5*(1-np.cos(np.pi*ind_aux/N_aux))
    s_aux_2 = 1*np.ones(N_aux)   #0.5*(1+np.cos(np.pi*ind_aux/N_aux))
    window = np.zeros(N_samples)
    ind_i_1 = int(N_samples/2-N_sym/2-N_aux/2)
    ind_f_1 = ind_i_1+int(N_aux)
    window[ind_i_1:ind_f_1] = s_aux_1
    ind_i_2 = int(N_samples/2+N_sym/2-N_aux/2)
    ind_f_2 = ind_i_2+int(N_aux)
    window[ind_i_2:ind_f_2] = s_aux_2
    window[ind_f_1:ind_i_2] = np.ones(ind_i_2-ind_f_1)
    S = np.fft.fftshift(np.fft.fft(s))
    S_f = S*window
    s_f = SpS*np.fft.ifft(np.fft.ifftshift(S_f))
    s = s_f
    #s = s_in

    if plot_flag:
        plt.figure()
        plt.plot(f/1e9,20*np.log10(np.abs(S)/max(np.abs(S))),label='signal')
        plt.plot(f/1e9,20*np.log10(window/max(np.abs(window))),label='filter')
        plt.ylabel('Normalized power / Filter resp. [dB]')
        plt.xlabel('Frequency [GHz]')

    # Generation of the sync sequence
    # Initial zeros
    s_zeros_init = np.zeros(10,dtype=complex)
    # Sync sequence
    np.random.seed(sync_seed)
    s_sync_i = np.random.randint(0,2,N_sync,dtype=int)*2-1
    # Final zeros
    s_zeros_final = np.zeros(10,dtype=complex)
    s_b_i = np.concatenate([s_zeros_init,s_sync_i,s_zeros_final])
    # Digital-to-analogue conversion of the synchronization sequence
    s_sync_i = dac_nyquist(symbols=s_b_i, samples_per_symbol=SpS, roll_off=RollOff)

    # Correlation between signal and synchronization
    s_xx_i = np.correlate(s,s_sync_i)
    if plot_flag:
        plt.figure()
        plt.plot(abs(s_xx_i))
        plt.ylabel('Cross-correlation')
        plt.xlabel('Sample index')
        plt.xlim([0,len(s_xx_i)])

    # Generating the index of the maximum correlation
    max_i = np.max(np.abs(s_xx_i))
    ind_i = np.where(max_i==np.abs(s_xx_i))
    # Finding the inidices of the sync and information
    ind_sync = (np.arange(N_sync)+10)*SpS+ind_i[0][0]
    ind_inf = (np.arange(N_inf)+10+N_sync)*SpS+ind_i[0][0]
    # Sampling the signal
    s_sync = s[ind_sync]
    s_inf = s[ind_inf]

    if plot_flag:
        # In-phase component time domain
        plt.figure()
        plt.plot(t*1e9,np.real(s))
        plt.plot(ind_sync*ts*1e9,np.real(s_sync),'.')
        plt.plot(ind_inf*ts*1e9,np.real(s_inf),'.')
        plt.xlabel('Time [ns]')
        plt.ylabel('In-phase component [a.u.]')
        plt.title('Raw data')
        plt.xlim([0,np.max(t)*1e9])

        # Contellation
        plt.figure()
        plt.plot(np.real(s_inf),np.imag(s_inf),'.')
        plt.plot(np.real(s_sync),np.imag(s_sync),'.')
        plt.axis('square')
        plt.xlabel('In-phase [a.u.]')
        plt.ylabel('Quadrature [a.u.]')
        plt.title('Raw data')

    # Finding the rotation angle
    theta_hat = np.angle(s_xx_i[ind_i])
    s = s*np.exp(-1j*theta_hat)
    s_sync = s_sync*np.exp(-1j*theta_hat)
    s_inf = s_inf*np.exp(-1j*theta_hat)

    if plot_flag:
        # In-phase component time domain
        plt.figure()
        plt.plot(t*1e9,np.real(s))
        plt.plot(ind_sync*ts*1e9,np.real(s_sync),'.')
        plt.plot(ind_inf*ts*1e9,np.real(s_inf),'.')
        plt.xlabel('Time [ns]')
        plt.ylabel('In-phase component [a.u.]')
        plt.title('After phase rotation')
        plt.xlim([0,np.max(t)*1e9])

        # Contellation
        plt.figure()
        plt.plot(np.real(s_inf),np.imag(s_inf),'.')
        plt.plot(np.real(s_sync),np.imag(s_sync),'.')
        plt.axis('square')
        plt.xlabel('In-phase [a.u.]')
        plt.ylabel('Quadrature [a.u.]')
        plt.title('After phase rotation')

    # Finding the scaling factor
    K = max_i/1940
    s = s/K
    s_sync = s_sync/K
    s_inf = s_inf/K
    if plot_flag:
        # In-phase component time domain
        plt.figure()
        plt.plot(t*1e9,np.real(s))
        plt.plot(ind_sync*ts*1e9,np.real(s_sync),'.')
        plt.plot(ind_inf*ts*1e9,np.real(s_inf),'.')
        plt.xlabel('Time [ns]')
        plt.ylabel('In-phase component [a.u.]')
        plt.title('After scaling correction')
        plt.xlim([0,np.max(t)*1e9])

        # Contellation
        plt.figure()
        plt.plot(np.real(s_inf),np.imag(s_inf),'.')
        plt.plot(np.real(s_sync),np.imag(s_sync),'.')
        plt.axis('square')
        plt.xlabel('In-phase [a.u.]')
        plt.ylabel('Quadrature [a.u.]')
        plt.title('After scaling correction')

    # Desnormalized the amplitude of the modulated signal
    if M == 16:
        s_inf_dn = 3*s_inf
    if M == 32:
        s_inf_dn = 5*s_inf

    # Blind phase noise and residual frequency offset compensation
    N_block = 16
    N_overlap = 8

    ind_1 = 0
    ind_2 = ind_1+N_block
    phi = 0
    min_dphi = -np.pi/4
    max_dphi = np.pi/4
    N_dphi = 21
    dphi_array = np.linspace(min_dphi,max_dphi,N_dphi)

    phi_array = np.array([0])
    ind_array = np.array([0])
    counter1 = 0
    while ind_2<=len(s_inf_dn):
        s_block = s_inf_dn[ind_1:ind_2]
        d_array_block = np.zeros(N_dphi)
        for counter, dphi in enumerate(dphi_array):
            s_block_rot = s_block*np.exp(1j*(dphi+phi))
            s_hat, dist = qam_dem(s_block_rot,M)
            d_array_block[counter] = np.var(dist)

        #plt.figure()
        #plt.plot(d_array_block)

        d_min = min(d_array_block)
        ind_min = np.where(d_array_block==d_min)
        dphi_opt = dphi_array[ind_min[0][0]]
        phi = phi+dphi_opt
        counter1 += 1
        phi_array = np.append(phi_array, phi)
        ind_array = np.append(ind_array, int((ind_2+ind_1)/2))
        ind_1 = ind_1+N_block-N_overlap
        ind_2 = ind_1+N_block

    ind_array = np.concatenate(([0],ind_array,[len(s_inf_dn)]))
    phi_array = np.concatenate(([0],phi_array,[phi_array[-1]]))
    phi_interp = np.interp(np.arange(len(s_inf_dn)),ind_array,phi_array)

    if plot_flag:
        plt.figure()
        plt.plot(phi_interp)
        #print(len(phi_interp))
    s_inf_pc = s_inf_dn*np.exp(1j*phi_interp)

    if plot_flag:
        # Contellation
        plt.figure()
        plt.plot(np.real(s_inf_pc),np.imag(s_inf_pc),'.')
        plt.axis('square')
        plt.xlabel('In-phase [a.u.]')
        plt.ylabel('Quadrature [a.u.]')
        plt.title('After scaling correction')

    #f plot_flag:
    #   constellation_analysis(s_inf_pc)
    #constellation_analysis(s_inf_pc)

    # Demapping
    s_hat, dist = QAM_dem(s_inf_pc,M)

    if plot_flag:
        plt.figure()
        plt.stem(s_hat!=s_b)

    BER = np.mean(s_hat!=s_b)

    return [s_inf, s_sync, s_inf_dn, BER]

"""### 6.1.3 Optical and electro-optical modules"""


