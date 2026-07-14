import numpy as np

def fiber(A=0,t=0,L=100e3,DeltaL=1e3,D=16e-6,lambda0=1550e-9,alpha_dB=0.2e-3,n2=2.6e-20,Aeff=80e-12,SoP_rotation=None):
    # E: input electrical double-polarization field in V/m
    # t: time vector in [s]
    # L: fiberl length in [m]
    # DeltaL: the longitudinal step value in [m]
    # D: dispersion parameter in [s/m^2]
    # lambda0: operation wavelength
    # alpha_dB: attenuation coefficient in dB/m
    # n2: nonlinear refractive index in [W-1]
    # Aeff: effective modal area in [mˆ2]

    # Constants
    c = 3e8             # Speed of light in the vacuum in m/s
    pi = np.pi          # Pi number

    # Preliminary calculations
    Nt = len(t)                    # Number of time samples
    Deltat = t[1]-t[0]             # Time step in s
    fs = 1/Deltat                  # Sampling frequency in Hz
    f = np.arange(Nt)/(Nt/fs)-fs/2 # Frequency vector in Hz
    Nz = int(L/DeltaL)             # Number of segments
    z = np.linspace(0,L-DeltaL,Nz) # Propagation distance vector in m
    alpha_Np = alpha_dB/4.343      # Attenuation of in Np/m
    k = 2*pi/lambda0               # Wave number in m-1
    gamma = n2*k/Aeff              # Nonlinear coefficient in W-1m-1

    # Extraction of the polarizations
    Ax = A[0,:]                    # We extract the x polarization
    Ay = A[1,:]                    # We extract the y polarization

    # Propagate through the fiber
    for auxz in range(len(z)):

        # Scalar linear effects are simulated in the frequency domain

        # Conversion to frequency domain
        Ax_f = np.fft.fftshift(np.fft.fft(Ax))
        Ay_f = np.fft.fftshift(np.fft.fft(Ay))

        # Attenuation
        Ax_f = Ax_f*np.exp(-alpha_Np/2*DeltaL)
        Ay_f = Ay_f*np.exp(-alpha_Np/2*DeltaL)

        # Dispersion
        Ax_f *= np.exp(-1j*(D*lambda0**2/(4*np.pi*c)*DeltaL*((2*np.pi)*f)**2))
        Ay_f *= np.exp(-1j*(D*lambda0**2/(4*np.pi*c)*DeltaL*((2*np.pi)*f)**2))

        # Conversion to time domain
        Ax = np.fft.ifft(np.fft.ifftshift(Ax_f))
        Ay = np.fft.ifft(np.fft.ifftshift(Ay_f))

        # Nonlinear effects are simulated in the time domain
        Ixy = (abs(Ax))**2+(abs(Ay))**2
        Ax *= np.exp(1j*gamma*Ixy*DeltaL);
        Ay *= np.exp(1j*gamma*Ixy*DeltaL);

        # Polarization rotation can be performed either in time or frequency domain
        if SoP_rotation:
            theta_rand = np.random.rand(1)*np.pi
            phi_rand = np.random.rand(1)*2*np.pi
            Ax_aux = Ax * np.exp(1j*phi_rand/2)
            Ay_aux = Ay * np.exp(-1j*phi_rand/2)
            Ax = np.cos(theta_rand)*Ax_aux+np.sin(theta_rand)*Ay_aux
            Ay = -np.sin(theta_rand)*Ax_aux+np.cos(theta_rand)*Ay_aux

    # Create the output field from the two polarizations
    A = np.zeros((2,len(t)),dtype = np.complex128)
    A[0,:] = Ax
    A[1,:] = Ay
    return A

# Function that simulates the fiber using the SSFM
