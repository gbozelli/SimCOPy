import numpy as np

# QAM_mod: QAM modulator for values of M of 4, 8, 16, and 32
def QAM_mod(s,M):

    # Preliminary calculations
    N_bits = int(np.log2(M))
    N_sym = int(len(s)/N_bits)
    s_I = np.zeros(N_sym,dtype=int)
    s_Q = np.zeros(N_sym,dtype=int)

    # N_bits = 1, M = 2
    if N_bits == 1:
        s_1 = np.reshape(s,[N_sym,N_bits])
        for ind_sym in range(N_sym):
            if s_1[ind_sym,:]==0:
                s_I[ind_sym], s_Q[ind_sym] = 0,1
            elif s_1[ind_sym,:]==1:
                s_I[ind_sym], s_Q[ind_sym] = 0,-1

    # N_bits = 2, M = 4
    if N_bits == 2:
        s_1 = np.reshape(s,[N_sym,N_bits])
        for ind_sym in range(N_sym):
            if all(s_1[ind_sym,:]==[0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -1,1
            elif all(s_1[ind_sym,:]==[0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,1
            elif all(s_1[ind_sym,:]==[1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-1
            elif all(s_1[ind_sym,:]==[1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-1

    # N_bits = 3, M = 8
    if N_bits == 3:
        s_1 = np.reshape(s,[N_sym,N_bits])
        for ind_sym in range(N_sym):
            if all(s_1[ind_sym,:]==[0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -3,1
            elif all(s_1[ind_sym,:]==[0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,1
            elif all(s_1[ind_sym,:]==[1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,1
            elif all(s_1[ind_sym,:]==[1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 3,1
            elif all(s_1[ind_sym,:]==[0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,-1
            elif all(s_1[ind_sym,:]==[0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-1
            elif all(s_1[ind_sym,:]==[1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-1
            elif all(s_1[ind_sym,:]==[1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,-1

    # N_bits = 4, M = 16
    if N_bits == 4:
        s_1 = np.reshape(s,[N_sym,N_bits])
        for ind_sym in range(N_sym):
            if all(s_1[ind_sym,:]==[0,0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -3,3
            elif all(s_1[ind_sym,:]==[0,1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,3
            elif all(s_1[ind_sym,:]==[1,1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,3
            elif all(s_1[ind_sym,:]==[1,0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 3,3
            elif all(s_1[ind_sym,:]==[0,0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,1
            elif all(s_1[ind_sym,:]==[0,1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -1,1
            elif all(s_1[ind_sym,:]==[1,1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 1,1
            elif all(s_1[ind_sym,:]==[1,0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,1
            elif all(s_1[ind_sym,:]==[0,0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,-1
            elif all(s_1[ind_sym,:]==[0,1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-1
            elif all(s_1[ind_sym,:]==[1,1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-1,
            elif all(s_1[ind_sym,:]==[1,0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,-1
            elif all(s_1[ind_sym,:]==[0,0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -3,-3
            elif all(s_1[ind_sym,:]==[0,1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-3
            elif all(s_1[ind_sym,:]==[1,1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-3
            elif all(s_1[ind_sym,:]==[1,0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 3,-3

    # N_bits = 5, M = 32
    if N_bits == 5:
        s_1 = np.reshape(s,[N_sym,N_bits])
        for ind_sym in range(N_sym):
            if all(s_1[ind_sym,:]==[1,0,1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -3,5
            elif all(s_1[ind_sym,:]==[1,0,1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,5
            elif all(s_1[ind_sym,:]==[1,1,1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,5
            elif all(s_1[ind_sym,:]==[1,1,1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 3,5
            elif all(s_1[ind_sym,:]==[1,0,1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -5,3
            elif all(s_1[ind_sym,:]==[0,0,1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,3
            elif all(s_1[ind_sym,:]==[0,0,1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,3
            elif all(s_1[ind_sym,:]==[0,1,1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,3
            elif all(s_1[ind_sym,:]==[0,1,1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,3
            elif all(s_1[ind_sym,:]==[1,1,1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 5,3
            elif all(s_1[ind_sym,:]==[1,0,1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -5,1
            elif all(s_1[ind_sym,:]==[0,0,1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,1
            elif all(s_1[ind_sym,:]==[0,0,1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,1
            elif all(s_1[ind_sym,:]==[0,1,1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,1
            elif all(s_1[ind_sym,:]==[0,1,1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,1
            elif all(s_1[ind_sym,:]==[1,1,1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 5,1
            elif all(s_1[ind_sym,:]==[1,0,0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -5,-1
            elif all(s_1[ind_sym,:]==[0,0,0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,-1
            elif all(s_1[ind_sym,:]==[0,0,0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-1
            elif all(s_1[ind_sym,:]==[0,1,0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-1
            elif all(s_1[ind_sym,:]==[0,1,0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,-1,
            elif all(s_1[ind_sym,:]==[1,1,0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 5,-1
            elif all(s_1[ind_sym,:]==[1,0,0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -5,-3
            elif all(s_1[ind_sym,:]==[0,0,0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,-3,
            elif all(s_1[ind_sym,:]==[0,0,0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-3
            elif all(s_1[ind_sym,:]==[0,1,0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-3,
            elif all(s_1[ind_sym,:]==[0,1,0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,-3
            elif all(s_1[ind_sym,:]==[1,1,0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 5,-3
            elif all(s_1[ind_sym,:]==[1,0,0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -3,-5
            elif all(s_1[ind_sym,:]==[1,0,0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-5,
            elif all(s_1[ind_sym,:]==[1,1,0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-5
            elif all(s_1[ind_sym,:]==[1,1,0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 3,-5
    return s_I+1j*s_Q