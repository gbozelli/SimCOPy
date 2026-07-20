import numpy as np

def edfa(
    input_signal: np.ndarray,
    t: np.ndarray,
    gain_db: float,
    noise_figure_db: float,
    wavelength: float = 1550e-9,
) -> np.ndarray:
    """
    Simulates an Erbium-Doped Fiber Amplifier (EDFA).
    """
    # Physical constants
    PLANCK_CONSTANT = 6.626e-34      # J·s
    SPEED_OF_LIGHT = 3e8             # m/s

    optical_frequency = SPEED_OF_LIGHT / wavelength

    sampling_period = t[1] - t[0]
    sampling_frequency = 1 / sampling_period
    simulation_bandwidth = sampling_frequency

    # Amplifier parameters
    linear_gain = 10 ** (gain_db / 10)
    noise_factor = 10 ** (noise_figure_db / 10)
    spontaneous_emission_factor = noise_factor / 2

    # Amplify signal
    output_signal = np.sqrt(linear_gain) * input_signal.copy()

    # ASE noise power
    ase_psd = (
        spontaneous_emission_factor
        * PLANCK_CONSTANT
        * optical_frequency
        * (linear_gain - 1)
    )

    ase_power = ase_psd * simulation_bandwidth

    num_samples = input_signal.shape[1]

    noise_std = np.sqrt(ase_power / 4)

    ase_noise = noise_std * (
        np.random.randn(2, num_samples)
        + 1j * np.random.randn(2, num_samples)
    )

    output_signal += ase_noise

    return output_signal