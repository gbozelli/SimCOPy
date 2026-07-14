import numpy as np
import numpy.matlib as ml

def QAM_dem(s,M):
    # Preliminary calculations
    N_bits = int(np.log2(M))
    N_sym = len(s)
    N_bits_total = N_bits*N_sym
    sb = np.zeros(N_bits_total,dtype=int)
    sd = np.zeros(N_sym)

    if N_bits == 2:
        # Find the closest constellation point
        const_points = np.array([-1+1j,1+1j,\
                                -1-1j,1-1j])
        sc_ext = ml.repmat(s,M,1)

        const_points_ext = np.transpose(ml.repmat(const_points,N_sym,1))

        #print(np.shape(sc_ext))
        #print(np.shape(sc_ext))
        sym_dist = np.abs(sc_ext-const_points_ext)
        min_dist = np.min(sym_dist,axis=0)
        #print(sc_ext)
        #print(const_points_ext)
        #print(sym_dist)

        ind_array = np.zeros(N_sym)
        for ind_sym in range(N_sym):
            ind = np.where(min_dist[ind_sym]==sym_dist[:,ind_sym])
            if ind[0][0] == 0:
                sb_aux = np.array([0,1])
            elif ind[0][0] == 1:
                sb_aux = np.array([0,0])
            elif ind[0][0] == 2:
                sb_aux = np.array([1,1])
            elif ind[0][0] == 3:
                sb_aux = np.array([1,0])
            ind_i = ind_sym*N_bits
            ind_f = (ind_sym+1)*N_bits
            sb[ind_i:ind_f] = sb_aux
            sd[ind_sym]=min_dist[ind_sym]

    if N_bits == 3:
        # Find the closest constellation point
        const_points = np.array([-3+1j,-1+1j,1+1j,3+1j,\
                                -3-1j,-1-1j,1-1j,3-1j])
        sc_ext = ml.repmat(s,M,1)

        const_points_ext = np.transpose(ml.repmat(const_points,N_sym,1))

        #print(np.shape(sc_ext))
        #print(np.shape(sc_ext))
        sym_dist = np.abs(sc_ext-const_points_ext)
        min_dist = np.min(sym_dist,axis=0)
        #print(sc_ext)
        #print(const_points_ext)
        #print(sym_dist)

        ind_array = np.zeros(N_sym)
        for ind_sym in range(N_sym):
            ind = np.where(min_dist[ind_sym]==sym_dist[:,ind_sym])
            if ind[0][0] == 0:
                sb_aux = np.array([0,0,0])
            elif ind[0][0] == 1:
                sb_aux = np.array([0,1,0])
            elif ind[0][0] == 2:
                sb_aux = np.array([1,1,0])
            elif ind[0][0] == 3:
                sb_aux = np.array([1,0,0])
            elif ind[0][0] == 4:
                sb_aux = np.array([0,0,1])
            elif ind[0][0] == 5:
                sb_aux = np.array([0,1,1])
            elif ind[0][0] == 6:
                sb_aux = np.array([1,1,1])
            elif ind[0][0] == 7:
                sb_aux = np.array([1,0,1])

            ind_i = ind_sym*N_bits
            ind_f = (ind_sym+1)*N_bits
            sb[ind_i:ind_f] = sb_aux
            sd[ind_sym]=min_dist[ind_sym]

    if N_bits == 4:
        # Find the closest constellation point
        const_points = np.array([-3+3j,-1+3j,1+3j,3+3j,\
                                -3+1j,-1+1j,1+1j,3+1j,\
                                -3-1j,-1-1j,1-1j,3-1j,\
                                -3-3j,-1-3j,1-3j,3-3j])
        sc_ext = ml.repmat(s,M,1)
        const_points_ext = np.transpose(ml.repmat(const_points,N_sym,1))
        sym_dist = np.abs(sc_ext-const_points_ext)
        min_dist = np.min(sym_dist,axis=0)
        #print(sc_ext)
        #print(const_points_ext)
        #print(sym_dist)
        #print(min_dist)
        ind_array = np.zeros(N_sym)
        for ind_sym in range(N_sym):
            #print(min_dist[ind_sym])
            #print(sym_dist[:,ind_sym])
            ind = np.where(min_dist[ind_sym]==sym_dist[:,ind_sym])
            #print(ind[0][0])
            #ind_array[ind_sym] = ind[0][0]
            if ind[0][0] == 0:
                sb_aux = np.array([0,0,1,0])
            elif ind[0][0] == 1:
                sb_aux = np.array([0,1,1,0])
            elif ind[0][0] == 2:
                sb_aux = np.array([1,1,1,0])
            elif ind[0][0] == 3:
                sb_aux = np.array([1,0,1,0])
            elif ind[0][0] == 4:
                sb_aux = np.array([0,0,1,1])
            elif ind[0][0] == 5:
                sb_aux = np.array([0,1,1,1])
            elif ind[0][0] == 6:
                sb_aux = np.array([1,1,1,1])
            elif ind[0][0] == 7:
                sb_aux = np.array([1,0,1,1])
            elif ind[0][0] == 8:
                sb_aux = np.array([0,0,0,1])
            elif ind[0][0] == 9:
                sb_aux = np.array([0,1,0,1])
            elif ind[0][0] == 10:
                sb_aux = np.array([1,1,0,1])
            elif ind[0][0] == 11:
                sb_aux = np.array([1,0,0,1])
            elif ind[0][0] == 12:
                sb_aux = np.array([0,0,0,0])
            elif ind[0][0] == 13:
                sb_aux = np.array([0,1,0,0])
            elif ind[0][0] == 14:
                sb_aux = np.array([1,1,0,0])
            elif ind[0][0] == 15:
                sb_aux = np.array([1,0,0,0])

            ind_i = ind_sym*N_bits
            ind_f = (ind_sym+1)*N_bits
            sb[ind_i:ind_f] = sb_aux
            sd[ind_sym]=min_dist[ind_sym]

    if N_bits == 5:
        # Find the closest constellation point
        const_points = np.array([-3+5j,-1+5j,1+5j,3+5j,\
                                -5+3j,-3+3j,-1+3j,1+3j,3+3j,5+3j,\
                                -5+1j,-3+1j,-1+1j,1+1j,3+1j,5+1j,\
                                -5-1j,-3-1j,-1-1j,1-1j,3-1j,5-1j,\
                                -5-3j,-3-3j,-1-3j,1-3j,3-3j,5-3j,\
                                -3-5j,-1-5j,1-5j,3-5j])
        sc_ext = ml.repmat(s,M,1)
        const_points_ext = np.transpose(ml.repmat(const_points,N_sym,1))
        sym_dist = np.abs(sc_ext-const_points_ext)
        min_dist = np.min(sym_dist,axis=0)
        #print(sc_ext)
        #print(const_points_ext)
        #print(sym_dist)
        #print(min_dist)
        ind_array = np.zeros(N_sym)
        for ind_sym in range(N_sym):
            #print(min_dist[ind_sym])
            #print(sym_dist[:,ind_sym])
            ind = np.where(min_dist[ind_sym]==sym_dist[:,ind_sym])
            #print(ind[0][0])
            #ind_array[ind_sym] = ind[0][0]
            if ind[0][0] == 0:
                sb_aux = np.array([1,0,1,0,0])
            elif ind[0][0] == 1:
                sb_aux = np.array([1,0,1,1,0])
            elif ind[0][0] == 2:
                sb_aux = np.array([1,1,1,1,0])
            elif ind[0][0] == 3:
                sb_aux = np.array([1,1,1,0,0])
            elif ind[0][0] == 4:
                sb_aux = np.array([1,0,1,1,1])
            elif ind[0][0] == 5:
                sb_aux = np.array([0,0,1,1,1])
            elif ind[0][0] == 6:
                sb_aux = np.array([0,0,1,1,0])
            elif ind[0][0] == 7:
                sb_aux = np.array([0,1,1,1,0])
            elif ind[0][0] == 8:
                sb_aux = np.array([0,1,1,1,1])
            elif ind[0][0] == 9:
                sb_aux = np.array([1,1,1,1,1])
            elif ind[0][0] == 10:
                sb_aux = np.array([1,0,1,0,1])
            elif ind[0][0] == 11:
                sb_aux = np.array([0,0,1,0,1])
            elif ind[0][0] == 12:
                sb_aux = np.array([0,0,1,0,0])
            elif ind[0][0] == 13:
                sb_aux = np.array([0,1,1,0,0])
            elif ind[0][0] == 14:
                sb_aux = np.array([0,1,1,0,1])
            elif ind[0][0] == 15:
                sb_aux = np.array([1,1,1,0,1])
            elif ind[0][0] == 16:
                sb_aux = np.array([1,0,0,0,1])
            elif ind[0][0] == 17:
                sb_aux = np.array([0,0,0,0,1])
            elif ind[0][0] == 18:
                sb_aux = np.array([0,0,0,0,0])
            elif ind[0][0] == 19:
                sb_aux = np.array([0,1,0,0,0])
            elif ind[0][0] == 20:
                sb_aux = np.array([0,1,0,0,1])
            elif ind[0][0] == 21:
                sb_aux = np.array([1,1,0,0,1])
            elif ind[0][0] == 22:
                sb_aux = np.array([1,0,0,1,1])
            elif ind[0][0] == 23:
                sb_aux = np.array([0,0,0,1,1])
            elif ind[0][0] == 24:
                sb_aux = np.array([0,0,0,1,0])
            elif ind[0][0] == 25:
                sb_aux = np.array([0,1,0,1,0])
            elif ind[0][0] == 26:
                sb_aux = np.array([0,1,0,1,1])
            elif ind[0][0] == 27:
                sb_aux = np.array([1,1,0,1,1])
            elif ind[0][0] == 28:
                sb_aux = np.array([1,0,0,0,0])
            elif ind[0][0] == 29:
                sb_aux = np.array([1,0,0,1,0])
            elif ind[0][0] == 30:
                sb_aux = np.array([1,1,0,1,0])
            elif ind[0][0] == 31:
                sb_aux = np.array([1,1,0,0,0])

            ind_i = ind_sym*N_bits
            ind_f = (ind_sym+1)*N_bits
            sb[ind_i:ind_f] = sb_aux
            sd[ind_sym]=min_dist[ind_sym]

    return [sb,sd]

# DAC_Nyquist: DAC considering Nyquist filtering
