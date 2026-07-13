import numpy as np

def optical_filter(s,t,f0,BW):
    ts = t[1]-t[0]                       # Sampling period in s
    fs = 1/ts                            # Sampling frequancy in Hz
    Deltaf = 1/(ts*len(t))               # Frequency resolution in Hz
    f = np.arange(len(t))*Deltaf-fs/2    # Frequency vector in Hz
    ind1 = np.where(f<(f0-BW/2))
    ind2 = np.where(f<(f0+BW/2))
    #print(ind1[0][-1])
    #print(ind2[0][-1])
    #print(f[ind1[0][-1]]/1e9)
    #print(f[ind2[0][-1]]/1e9)

    S = np.fft.fftshift(np.fft.fft(s))

    H = np.zeros(len(f))
    H[ind1[0][-1]:ind2[0][-1]]=1
    Sf = S*H

    #plt.figure()
    #plt.plot(f,np.abs(S[0,:]))
    #plt.plot(f,np.abs(Sf[0,:]))

    sf = np.fft.ifft(np.fft.ifftshift(Sf))

    return sf
