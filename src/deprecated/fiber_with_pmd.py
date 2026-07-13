import numpy as np
from scipy.linalg import expm # Often used for matrix exponentials, though not strictly needed for this specific implementation

# Function that simulates the fiber using the SSFM, now with PMD
def fiber_with_pmd(A=0, t=0, L=100e3, DeltaL=1e3, D=16e-6, lambda0=1550e-9, alpha_dB=0.2e-3, n2=2.6e-20, Aeff=80e-12, SoP_rotation=None, pmd_coefficient=0):
    """
    Simulates optical fiber propagation using the Split-Step Fourier Method,
    including Chromatic Dispersion, Kerr Effect (SPM), and PMD.

    Args:
        A (np.ndarray): Input electric field [2xN] for X and Y polarizations.
        t (np.ndarray): Time vector.
        L (float): Total fiber length in meters.
        DeltaL (float): Length of each simulation step in meters.
        D (float): Chromatic dispersion parameter in s/m^2.
        lambda0 (float): Central wavelength in meters.
        alpha_dB (float): Attenuation in dB/m.
        n2 (float): Nonlinear refractive index in W^-1.
        Aeff (float): Effective modal area in m^2.
        SoP_rotation (bool): Legacy flag, functionality is now handled by the PMD model.
        pmd_coefficient (float): PMD coefficient in s/m^(1/2). A value of 0 disables PMD.

    Returns:
        np.ndarray: Output electric field [2xN].
    """
    # Constants
    c = 299792458.0
    pi = np.pi

    # Preliminary calculations
    Nt = len(t)
    Deltat = t[1] - t[0]
    fs = 1 / Deltat
    f = np.fft.fftshift(np.fft.fftfreq(Nt, d=Deltat)) # Use fftshift and fftfreq for correct frequency vector
    w = 2 * pi * f
    Nz = int(L / DeltaL)
    alpha_Np = alpha_dB / 4.343
    k = 2 * pi / lambda0
    gamma = n2 * k / Aeff if n2 > 0 else 0
    beta2 = -D * lambda0**2 / (2 * pi * c) # Dispersion in terms of beta2

    # --- Linear Operator (Dispersion and Attenuation) ---
    # This is applied in the frequency domain. It's more efficient to define it once.
    linear_operator = np.exp(-alpha_Np/2 * DeltaL + 0.5j * beta2 * w**2 * DeltaL)

    # Extraction of the polarizations
    Ax = A[0, :]
    Ay = A[1, :]

    # Propagate through the fiber
    for _ in range(Nz):
        # --- First half of linear step ---
        Ax_f = np.fft.fft(Ax)
        Ay_f = np.fft.fft(Ay)
        
        Ax_f *= linear_operator
        Ay_f *= linear_operator

        Ax = np.fft.ifft(Ax_f)
        Ay = np.fft.ifft(Ay_f)

        # ############################################################### #
        # ### NEW PMD SECTION ########################################### #
        # ############################################################### #
        if pmd_coefficient > 0:
            # 1. Calculate the mean DGD for this segment.
            mean_dgd_step = pmd_coefficient * np.sqrt(DeltaL)

            # 2. Draw a random DGD value for this step from a Maxwellian
            # distribution. We approximate this by taking the absolute value of a Gaussian.
            dgd_step = np.abs(np.random.normal(loc=0, scale=mean_dgd_step / np.sqrt(3)))
            
            # 3. Create a random rotation matrix (Jones Matrix) to orient the DGD.
            theta = np.random.uniform(0, np.pi)
            phi = np.random.uniform(0, 2 * np.pi)
            
            J11 = np.cos(theta/2)
            J12 = np.sin(theta/2) * np.exp(1j * phi)
            J21 = -np.sin(theta/2) * np.exp(-1j * phi)
            J22 = np.cos(theta/2)

            # 4. Apply the rotation in the time domain.
            Ax_rot = J11 * Ax + J12 * Ay
            Ay_rot = J21 * Ax + J22 * Ay
            
            # 5. Apply the DGD in the frequency domain (it's a phase shift).
            # This is the operator for the delay.
            dgd_operator_fast = np.exp(1j * w * dgd_step / 2)
            dgd_operator_slow = np.exp(-1j * w * dgd_step / 2)

            Ax_rot_f = np.fft.fft(Ax_rot) * dgd_operator_fast
            Ay_rot_f = np.fft.fft(Ay_rot) * dgd_operator_slow
            
            # 6. Convert back to time domain. The inverse rotation is not needed
            # because the next segment will have its own new random rotation.
            Ax = np.fft.ifft(Ax_rot_f)
            Ay = np.fft.ifft(Ay_rot_f)

        # --- Full Nonlinear Step (Kerr Effect) ---
        # This is applied in the time domain
        if gamma > 0:
            Ixy = (np.abs(Ax))**2 + (np.abs(Ay))**2
            phase_shift = np.exp(1j * gamma * Ixy * DeltaL)
            Ax *= phase_shift
            Ay *= phase_shift

        # --- Second half of linear step ---
        Ax_f = np.fft.fft(Ax)
        Ay_f = np.fft.fft(Ay)

        Ax_f *= linear_operator
        Ay_f *= linear_operator
        
        Ax = np.fft.ifft(Ax_f)
        Ay = np.fft.ifft(Ay_f)
    
    # Create the output field from the two polarizations
    A_out = np.zeros_like(A)
    A_out[0, :] = Ax
    A_out[1, :] = Ay
    return A_out