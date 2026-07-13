import matplotlib.pyplot as plt
import numpy as np
from photodiode import photodiode
from single_frequency_laser import single_frequency_laser
from HybridNetwork import HybridNetwork
from QAM_receiver_DP import QAM_receiver_DP

def DPQAM_receiver(E_RX,t,M,SpS,N_inf,RollOff,ts, s_b, channel_id=None, power_dBm=None,
                   P_laser=0.01, Delta_nu=0, N_sync=256,
                   sync_seed_X=0, sync_seed_Y=123,
                   Freq_offset=0,
                   plot_flag=False, BW=10e9, shot_noise=True,
                   thermal_noise=True, L=0, D=0, R=0.9):

    #------------------------------------- Optical front-end --------------------------------------------%

    E_RX_X = E_RX[0,:]
    E_RX_Y = E_RX[1,:]

    E_LO = single_frequency_laser(t=t, P=P_laser, Delta_nu=Delta_nu, Freq_offset=Freq_offset, X_fraction=0.5, phi_pol=0)
    E_LO_X, E_LO_Y = E_LO[0,:], E_LO[1,:]

    # Polarização X
    E_1, E_2, E_3, E_4 = HybridNetwork(E_RX_X, E_LO_X)
    ipd_1 = photodiode(s=E_1,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_2 = photodiode(s=E_2,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_3 = photodiode(s=E_3,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_4 = photodiode(s=E_4,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    s_X_I = ipd_1 - ipd_2
    s_X_Q = ipd_3 - ipd_4

    # Polarização Y
    E_1, E_2, E_3, E_4 = HybridNetwork(E_RX_Y, E_LO_Y)
    ipd_1 = photodiode(s=E_1,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_2 = photodiode(s=E_2,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_3 = photodiode(s=E_3,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_4 = photodiode(s=E_4,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    s_Y_I = ipd_1 - ipd_2
    s_Y_Q = ipd_3 - ipd_4

    # Conversão cartesiana → complexa
    s_rx_X = s_X_I + 1j*s_X_Q
    s_rx_Y = s_Y_I + 1j*s_Y_Q

    s_rec = np.vstack([s_rx_X, s_rx_Y])

    # DSP e decisão de símbolos
    s_inf, s_sync, s_inf_dn, BER, sx_inf_pc, sy_inf_pc = QAM_receiver_DP(
        s_rec, t=t, M=M, SpS=SpS, N_sync=N_sync, N_inf=N_inf,
        sync_seed_X=0, sync_seed_Y=123, RollOff=RollOff, ts=ts, s_b=s_b,
        plot_flag=False, L=L, D=D
    )

    const_rx = np.zeros((2,len(sx_inf_pc)),dtype=complex)
    const_rx[0,:] = sx_inf_pc
    const_rx[1,:] = sy_inf_pc

    #Antes estava apenas return BER
    return [BER, const_rx]

