import numpy as np
import matplotlib.pyplot as plt
from DAC_Nyquist import DAC_Nyquist
from QAM_dem import QAM_dem
from photodiode import photodiode
from HybridNetwork import HybridNetwork
from single_frequency_laser import single_frequency_laser
from spectrum_analysis import spectrum_analysis

def QAM_receiver(s,M,N_sync=128,N_inf=4096,SpS=16,RollOff=0.2,sync_seed=0,ts=1e-12,plot_flag=False,s_b=None):
    s_in = s

    # Generation of time vector
    N_samples = len(s)
    t = ts*(np.arange(N_samples))
    fs = 1/ts
    f = (np.arange(N_samples)/N_samples)*fs-fs/2

    # Coarse frequency offset compensation

    N_samples = len(s)
    N_sym = int(N_samples/SpS)
    S_nyq = np.zeros(N_samples,dtype=complex)
    ind_i = int(N_samples/2-N_sym/2)
    ind_f = int(ind_i+N_sym)
    window = np.zeros(N_samples)
    window[ind_i:ind_f] = np.ones(ind_f-ind_i)

    f_min, f_max, N_f = -500e6,500e6,101
    f_offset_array = np.linspace(f_min,f_max,N_f)
    P_array = np.linspace(f_min,f_max,N_f)

    for counter, offset in enumerate(f_offset_array):
        s_shift = np.exp(1j*2*np.pi*offset*t)*s
        S_shift = np.abs(np.fft.fftshift(np.fft.fft(s_shift)))
        S_f = S_shift*window
        P_array[counter] = sum(np.abs(S_f)**2)

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
    s_sync_i = DAC_Nyquist(s=s_b_i,SpS=SpS,RollOff=RollOff)

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
            s_hat, dist = QAM_dem(s_block_rot,M)
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
def QAM_receiver_DP(s,t,s_b,M,N_sync=128,N_MIMO=32,N_inf=4096,SpS=16,RollOff=0.2,sync_seed_X=0,sync_seed_Y=123,ts=1e-12,plot_flag=False,L=10e3,D=16e-6):
    #print('The length is',L)
    #print('The sipersion parameters is',D)
    s_in = s
    s_x = s[0,:]
    s_y = s[1,:]

    # Generation of time vector
    N_samples = len(s_x)
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
    P_array = (P_array_x+P_array_y)/2

    Npolord = 10
    max_offset_opt_x = max(P_array_x)
    ind_offset_opt_x = np.where(P_array_x==max_offset_opt_x)
    f_offset_opt_x = f_offset_array[ind_offset_opt_x[0][0]]
    p_x = np.polyfit(f_offset_array,P_array_x,Npolord)
    max_offset_opt_y = max(P_array_y)
    ind_offset_opt_y = np.where(P_array_y==max_offset_opt_y)
    f_offset_opt_y = f_offset_array[ind_offset_opt_y[0][0]]
    p_y = np.polyfit(f_offset_array,P_array_y,Npolord)
    max_offset_opt = max(P_array)
    ind_offset_opt = np.where(P_array==max_offset_opt)
    f_offset_opt = f_offset_array[ind_offset_opt[0][0]]
    p = np.polyfit(f_offset_array,P_array,Npolord)

    N_f = 2001
    f_offset_array_ext = np.linspace(f_min,f_max,N_f)
    P_array_ext_x = np.polyval(p_x,f_offset_array_ext)
    P_array_ext_y = np.polyval(p_y,f_offset_array_ext)
    P_array_ext = np.polyval(p,f_offset_array_ext)

    max_offset_opt_x = max(P_array_ext_x)
    ind_offset_opt_x = np.where(P_array_ext_x==max_offset_opt_x)
    f_offset_opt_x = f_offset_array_ext[ind_offset_opt_x[0][0]]

    max_offset_opt_y = max(P_array_ext_y)
    ind_offset_opt_y = np.where(P_array_ext_y==max_offset_opt_y)
    f_offset_opt_y = f_offset_array_ext[ind_offset_opt_y[0][0]]

    max_offset_opt = max(P_array_ext)
    ind_offset_opt = np.where(P_array_ext==max_offset_opt)
    f_offset_opt = f_offset_array_ext[ind_offset_opt[0][0]]

    if plot_flag:
        plt.figure()
        plt.plot(f_offset_array/1e6,P_array_x,'o')
        plt.plot(f_offset_array/1e6,P_array_y,'o')
        plt.plot(f_offset_array/1e6,P_array,'o')
        plt.plot(f_offset_array_ext/1e6,P_array_ext_x)
        plt.plot(f_offset_array_ext/1e6,P_array_ext_y)
        plt.plot(f_offset_array_ext/1e6,P_array_ext)
        plt.plot(f_offset_array_ext[ind_offset_opt_x[0][0]]/1e6,max_offset_opt_x,'^r')
        plt.plot(f_offset_array_ext[ind_offset_opt_y[0][0]]/1e6,max_offset_opt_y,'^r')
        plt.plot(f_offset_array_ext[ind_offset_opt[0][0]]/1e6,max_offset_opt,'^r')
        plt.ylabel('Power [a.u.]')
        plt.xlabel('Frequency offset [MHz]')

    #print('The frequency offset for x is:',f_offset_opt_x/1e6)
    #print('The frequency offset for y is:',f_offset_opt_y/1e6)
    #print('The frequency offset for x-y is:',f_offset_opt/1e6)
    s = np.exp(1j*2*np.pi*f_offset_opt*t)*s
    #s = s_in

    # Filter the signal using a square filter
    N_samples = len(s[0,:])
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

    sx = s[0,:]
    Sx = np.fft.fftshift(np.fft.fft(sx))
    S_f_x = Sx*window
    s_f_x = SpS*np.fft.ifft(np.fft.ifftshift(S_f_x))
    sy = s[1,:]
    Sy = np.fft.fftshift(np.fft.fft(sy))
    S_f_y = Sy*window
    s_f_y = SpS*np.fft.ifft(np.fft.ifftshift(S_f_y))

    s[0,:] = s_f_x
    s[1,:] = s_f_y
    #s = s_in

    #s = s_in

    if plot_flag:
        epsilon = 1e-12
        plt.figure()
        plt.plot(f/1e9,20*np.log10(np.abs(Sx)/max(np.abs(Sx))+epsilon),label='signal')
        plt.plot(f/1e9,20*np.log10(np.abs(Sy)/max(np.abs(Sy))+epsilon),label='signal')
        plt.plot(f/1e9,20*np.log10(np.abs(S_f_x)/max(np.abs(S_f_x))+epsilon),label='signal')
        plt.plot(f/1e9,20*np.log10(np.abs(S_f_y)/max(np.abs(S_f_y))+epsilon),label='signal')
        plt.plot(f/1e9,20*np.log10(window/max(np.abs(window))+epsilon),label='filter')
        plt.ylabel('Normalized power / Filter resp. [dB]')
        plt.xlabel('Frequency [GHz]')

    ####### Dispersion compensation #########

    # Conversion to frequency domain
    s_f_x_f = np.fft.fftshift(np.fft.fft(s_f_x))
    s_f_y_f = np.fft.fftshift(np.fft.fft(s_f_y))

    # Dispersion
    c = 3e8
    lambda0 = 1550e-9
    s_f_x_f *= np.exp(1j*(D*lambda0**2/(4*np.pi*c)*L*((2*np.pi)*f)**2))
    s_f_y_f *= np.exp(1j*(D*lambda0**2/(4*np.pi*c)*L*((2*np.pi)*f)**2))

    # Conversion to time domain
    s_f_x = np.fft.ifft(np.fft.ifftshift(s_f_x_f))
    s_f_y = np.fft.ifft(np.fft.ifftshift(s_f_y_f))
    s[0,:] = s_f_x
    s[1,:] = s_f_y

    # Generation of the sync sequence
    # Initial zeros
    s_zeros_init = np.zeros(10,dtype=complex)
    # Sync sequence
    np.random.seed(sync_seed_X)
    s_sync_i_x_d = np.random.randint(0,2,N_sync,dtype=int)*2-1
    np.random.seed(sync_seed_Y)
    s_sync_i_y_d = np.random.randint(0,2,N_sync,dtype=int)*2-1

    # Final zeros
    s_zeros_final = np.zeros(10,dtype=complex)
    s_b_i_x = np.concatenate([s_zeros_init,s_sync_i_x_d,s_zeros_final])
    s_b_i_y = np.concatenate([s_zeros_init,s_sync_i_y_d,s_zeros_final])

    # Digital-to-analogue conversion of the synchronization sequence
    s_sync_i_x = DAC_Nyquist(s=s_b_i_x,SpS=SpS,RollOff=RollOff)
    s_sync_i_y = DAC_Nyquist(s=s_b_i_y,SpS=SpS,RollOff=RollOff)

    # Correlation between signal and synchronization
    sx = s[0,:]
    sy = s[1,:]
    s_xx_i = np.correlate(sx,s_sync_i_x)
    s_xy_i = np.correlate(sx,s_sync_i_y)
    s_yy_i = np.correlate(sy,s_sync_i_y)
    s_yx_i = np.correlate(sy,s_sync_i_x)

    if plot_flag:
        plt.figure()
        plt.plot(abs(s_xx_i))
        plt.plot(abs(s_xy_i))
        plt.plot(abs(s_yy_i))
        plt.plot(abs(s_yx_i))
        plt.ylabel('Cross-correlation')
        plt.xlabel('Sample index')
        plt.xlim([0,len(s_xx_i)])
    s_corr_total = abs(s_xx_i)+abs(s_xy_i)+abs(s_yy_i)+abs(s_yx_i)

    #print(N_sync)
    #kaixo

    # Generating the index of the maximum correlation
    max_i = np.max(s_corr_total)
    ind_i = np.where(max_i==s_corr_total)
    # Finding the inidices of the sync and information
    ind_sync = (np.arange(N_sync)+10)*SpS+ind_i[0][0]
    ind_MIMO = (np.arange(N_MIMO)+10+N_sync)*SpS+ind_i[0][0]
    ind_inf = (np.arange(N_inf)+10+N_sync)*SpS+ind_i[0][0]
    #print(N_sync)

    #print(max(abs(s_xx_i[ind_i])))
    #print(max(abs(s_xy_i[ind_i])))
    #print(max(abs(s_yy_i[ind_i])))
    #print(max(abs(s_yx_i[ind_i])))
    #print(max(np.angle(s_xx_i[ind_i])))
    #print(max(np.angle(s_xy_i[ind_i])))
    #print(max(np.angle(s_yy_i[ind_i])))
    #print(max(np.angle(s_yx_i[ind_i])))

    H_hat = np.zeros((2,2),dtype=complex)
    H_hat[0][0] = s_xx_i[ind_i][0]
    H_hat[0][1] = s_xy_i[ind_i][0]
    H_hat[1][0] = s_yx_i[ind_i][0]
    H_hat[1][1] = s_yy_i[ind_i][0]
    H_inv = np.linalg.inv(H_hat)

    #H_hat = np.matrix([[s_xx_i[ind_i][0],s_xy_i[ind_i]][0],[s_yx_i[ind_i][0],s_yy_i[ind_i][0]]])
    #print(H_hat)
    #print(H_inv)

    s_sync = s[:,ind_sync]
    s_sync_x = sx[ind_sync]
    s_sync_y = sy[ind_sync]
    s_inf = s[:,ind_inf]
    s_inf_x = s_inf[0,:]
    s_inf_y = s_inf[1,:]

    s_sync_corr = np.matmul(H_inv,s_sync)
    s_sync_corr_x = s_sync_corr[0,:]
    s_sync_corr_y = s_sync_corr[1,:]

    s_inf_corr = np.matmul(H_inv,s_inf)
    s_inf_corr_x = s_inf_corr[0,:]
    s_inf_corr_y = s_inf_corr[1,:]

    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(np.real(s_sync_corr_x),'.')
    plt.plot(np.imag(s_sync_corr_x),'.')
    plt.plot(np.real(s_inf_corr_x),'.')
    plt.plot(np.imag(s_inf_corr_x),'.')

    plt.subplot(2,1,2)
    plt.plot(np.real(s_sync_corr_y),'.')
    plt.plot(np.imag(s_sync_corr_y),'.')
    plt.plot(np.real(s_inf_corr_y),'.')
    plt.plot(np.imag(s_inf_corr_y),'.')

    if plot_flag:
        # In-phase component time domain
        plt.figure()
        plt.plot(t*1e9,np.real(sx))
        plt.plot(ind_sync*ts*1e9,np.real(s_sync_x),'.')
        plt.plot(ind_inf*ts*1e9,np.real(s_inf_x),'.')
        plt.xlabel('Time [ns]')
        plt.ylabel('In-phase component [a.u.]')
        plt.title('Raw data')
        plt.xlim([0,np.max(t)*1e9])

        # Contellation
        plt.figure()
        plt.plot(np.real(s_inf_corr_x),np.imag(s_inf_corr_x),'.')
        plt.plot(np.real(s_sync_corr_x),np.imag(s_sync_corr_x),'.')
        plt.axis('square')
        plt.xlabel('In-phase [a.u.]')
        plt.ylabel('Quadrature [a.u.]')
        plt.title('Raw data')


    K_amp = np.mean(np.abs(s_sync_corr_x))
    #print(K_amp)

    s_sync_corr_x = s_sync_corr_x/K_amp
    s_sync_corr_y = s_sync_corr_y/K_amp
    s_inf_corr_x = s_inf_corr_x/K_amp
    s_inf_corr_y = s_inf_corr_y/K_amp

    # Contellation
    plt.figure()
    plt.plot(np.real(s_inf_corr_x),np.imag(s_inf_corr_x),'.')
    plt.plot(np.real(s_sync_corr_x),np.imag(s_sync_corr_x),'.')
    plt.axis('square')
    plt.xlabel('In-phase [a.u.]')
    plt.ylabel('Quadrature [a.u.]')
    plt.title('Raw data')

    # Desnormalized the amplitude of the modulated signal
    if M == 16:
        sx_inf_dn = 3*s_inf_corr_x
        sy_inf_dn = 3*s_inf_corr_y
    if M == 32:
        sx_inf_dn = 5*s_inf_corr_x
        sy_inf_dn = 5*s_inf_corr_y
    s_inf_dn = np.append(sx_inf_dn,sy_inf_dn)

    # Blind phase noise and residual frequency offset compensation
    #N_block = 8
    #N_overlap = 4
    N_block = 64
    N_overlap = 32
    #N_block = 32
    #N_overlap = 26

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
    while ind_2<=len(sx_inf_dn):
        s_block = sx_inf_dn[ind_1:ind_2]
        d_array_block = np.zeros(N_dphi)
        for counter, dphi in enumerate(dphi_array):
            s_block_rot = s_block*np.exp(1j*(dphi+phi))
            s_hat, dist = QAM_dem(s_block_rot,M)
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

    ind_array = np.concatenate(([0],ind_array,[len(sx_inf_dn)]))
    phi_array = np.concatenate(([0],phi_array,[phi_array[-1]]))
    phi_interp = np.interp(np.arange(len(sx_inf_dn)),ind_array,phi_array)

    if plot_flag:
        plt.figure()
        plt.plot(phi_interp)
        #print(len(phi_interp))
    sx_inf_pc = sx_inf_dn*np.exp(1j*(phi_interp))
    sy_inf_pc = sy_inf_dn*np.exp(1j*(phi_interp))

    if plot_flag:
        # Contellation
        Ni = 10
        plt.figure()
        plt.plot(np.real(sx_inf_pc[Ni:]),np.imag(sx_inf_pc[Ni:]),'.')
        plt.plot(np.real(sy_inf_pc[Ni:]),np.imag(sy_inf_pc[Ni:]),'.')
        plt.axis('square')
        plt.xlabel('In-phase [a.u.]')
        plt.ylabel('Quadrature [a.u.]')
        plt.title('After scaling correction')

    # Demapping
    s_hat_x, dist = QAM_dem(sx_inf_pc,M)
    s_hat_y, dist = QAM_dem(sy_inf_pc,M)

    s_b_x = s_b[0,:]
    s_b_y = s_b[1,:]
    #print(s_b_x)
    #print(s_b_y)

    if plot_flag:
        plt.figure()
        plt.stem(s_hat_x!=s_b_x)
        plt.stem(s_hat_y!=s_b_y)

    BER_X = np.mean(s_hat_x!=s_b_x)
    BER_Y = np.mean(s_hat_y!=s_b_y)
    BER = np.array([BER_X,BER_Y])

    return [s_inf, s_sync, s_inf_dn, BER]

# Optical receiver
def DPQAM_receiver(E_RX,E_TX,t,RollOff,ts,M,SpS,N_inf, s_b, P_laser = 0.01, Delta_nu = 100e3, N_sync= 128, sync_seed_X=0, sync_seed_Y=123,Freq_offset = 1e6, plot_flag = False,\
                  R=0.9, BW=4*20e9, shot_noise = True, thermal_noise = True,L=10e3,D=16e-6):

    #------------------------------------- Optical front end --------------------------------------------%

    spectrum_analysis(E_RX[0,:],t)

    # Separate polarizations
    E_RX_X = E_RX[0,:]
    E_RX_Y = E_RX[1,:]

    spectrum_analysis(E_TX[0,:],t)
    # Filter the desired channel



    # LO laser
    E_LO = single_frequency_laser(t = t, P = P_laser, Delta_nu = Delta_nu, Freq_offset = Freq_offset, X_fraction = 0.5, phi_pol = 0)

    # Separate the two polarizations of LO
    E_LO_X = E_LO[0,:]
    E_LO_Y = E_LO[1,:]

    # Process X polarization
    # Combine signal and LO
    E_1, E_2, E_3, E_4 = HybridNetwork(E_RX_X,E_LO_X)
    # Detect the outputs of 90-degree hybrid
    ipd_1 = photodiode(s=E_1,t=t,R=0.9,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_2 = photodiode(s=E_2,t=t,R=0.9,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_3 = photodiode(s=E_3,t=t,R=0.9,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_4 = photodiode(s=E_4,t=t,R=0.9,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    s_X_I = ipd_1-ipd_2   # Differential amplification
    s_X_Q = ipd_3-ipd_4   # Differential amplification

    spectrum_analysis(s_X_I,t)

    # Process Y polarization
    # Combine signal and LO
    E_1, E_2, E_3, E_4 = HybridNetwork(E_RX_Y,E_LO_Y)
    # Detect the outputs of 90-degree hybrid
    ipd_1 = photodiode(s=E_1,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_2 = photodiode(s=E_2,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_3 = photodiode(s=E_3,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_4 = photodiode(s=E_4,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    s_Y_I = ipd_1-ipd_2  # Differential amplification
    s_Y_Q = ipd_3-ipd_4  # Differential amplification

    # Cartesian to complex conversion
    s_rx_X = s_X_I+1j*s_X_Q
    s_rx_Y = s_Y_I+1j*s_Y_Q

    #print(np.shape(s_Y_I))

    m,n = np.shape(s_Y_I)
    s_rec = np.zeros((2,n),dtype=complex)

    #print(np.shape(s_Y_I))
    #print(np.shape(s_rec))
    #print(np.shape(s_rx_X))

    s_rec[0,:] = s_rx_X
    s_rec[1,:] = s_rx_Y


    #QAM_receiver_DP(s_rec, t = t, M =M, SpS = SpS, N_sync = N_sync, N_inf = N_inf,sync_seed_X=0,sync_seed_Y=123,\
    #                RollOff = RollOff,ts = ts,s_b=s_b[0,:],plot_flag = plot_flag)
    s_inf, s_sync, s_inf_dn, BER = QAM_receiver_DP(s_rec, t = t, M =M, SpS = SpS, N_sync = N_sync, N_inf = N_inf,sync_seed_X=0,sync_seed_Y=123,\
                    RollOff = RollOff,ts = ts,s_b=s_b,plot_flag = plot_flag,L=L,D=D)
    return BER