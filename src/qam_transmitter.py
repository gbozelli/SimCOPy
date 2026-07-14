import numpy as np
import matplotlib.pyplot as plt
from qam_modulation import qam_modulation
from dac_nyquist import dac_nyquist

def QAM_transmitter(N_sync=128,sync_seed=0,\
                    inf_seed=45,\
                    N_MIMO=32,MIMO_i=0,\
                    N_inf=1024,M=16,SpS=16,RollOff=0.2,ts=1e-12,
                    N_zeros_init=20,N_zeros_final=20,plot_flag=True):

    ## Sequence generation

    # Initial zeros
    s_zeros_init = np.zeros(N_zeros_init)
    # Syncronization sequence
    np.random.seed(sync_seed)
    s_sync_i = np.random.randint(0,2,N_sync,dtype=int)*2-1
    # MIMO sequence
    #print(MIMO_i)
    if MIMO_i == 1:
        s_MIMO_i = np.append(np.zeros(int(N_MIMO/2)),np.ones(int(N_MIMO/2)))
    else:
        s_MIMO_i = np.append(np.ones(int(N_MIMO/2)),np.zeros(int(N_MIMO/2)))
    s_MIMO_i = s_MIMO_i*np.mod(np.arange(N_MIMO)+MIMO_i,2)
    #print(s_MIMO_i)

    # Information sequence
    N_bits = int(np.log2(M))
    N_bits_total = int(N_bits*N_inf)
    np.random.seed(inf_seed)
    s_b = np.random.randint(0,2,N_bits_total,dtype=int)
    s_inf = qam_modulation(s_b,M)
    const_tx = s_inf.copy()
    s_inf = s_inf/np.max(np.real(s_inf))

    # Final zeros
    s_zeros_final = np.zeros(N_zeros_final)

    # Concatenate
    s = np.concatenate([s_zeros_init,s_sync_i,s_inf,s_zeros_final])

    ## Conversion to the continuous time
    s_cont = dac_nyquist(symbols=s, samples_per_symbol=SpS, roll_off=RollOff)

    ind = np.arange(len(s))*SpS

    t = ts*(np.arange(len(s_cont)))

    if plot_flag:
        plt.figure()
        plt.plot(t*1e9,np.real(s_cont))
        plt.plot(ind*ts*1e9,np.real(s),'.')
        plt.xlabel('Time [ns]')
        plt.ylabel('In-phase component [a.u.]')
        plt.xlim([0,np.max(t)*1e9])

        plt.figure()
        plt.plot(np.real(s_cont),np.imag(s_cont),alpha=0.2)
        plt.plot(np.real(s),np.imag(s),'.')
        plt.xlabel('In-phase')
        plt.ylabel('Quadrature')
        plt.axis('square')

    return [s_cont,t,s_b,const_tx]


