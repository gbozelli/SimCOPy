import importlib.util
import os
import numpy as np
import pytest

# Ensure pytest has a stable root location
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
src_dir = os.path.join(repo_root, 'src')
deprecated_dir = os.path.join(repo_root, 'src', 'deprecated')

import matplotlib
matplotlib.use('Agg')


def load_module_from_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

EDFA_old = load_module_from_path('EDFA_old', os.path.join(deprecated_dir, 'EDFA.py'))
EDFA_new = load_module_from_path('EDFA_new', os.path.join(src_dir, 'edfa.py'))
HybridNetwork_old = load_module_from_path('HybridNetwork_old', os.path.join(deprecated_dir, 'HybridNetwork.py'))
HybridNetwork_new = load_module_from_path('HybridNetwork_new', os.path.join(src_dir, 'hybrid_network.py'))
MZM_old = load_module_from_path('MZM_old', os.path.join(deprecated_dir, 'MZM.py'))
MZM_new = load_module_from_path('MZM_new', os.path.join(src_dir, 'mzm.py'))
optical_filter_old = load_module_from_path('optical_filter_old', os.path.join(deprecated_dir, 'optical_filter.py'))
optical_filter_new = load_module_from_path('optical_filter_new', os.path.join(src_dir, 'optical_filter.py'))
constellation_analysis_old = load_module_from_path('constellation_analysis_old', os.path.join(deprecated_dir, 'constellation_analysis.py'))
constellation_analysis_new = load_module_from_path('constellation_analysis_new', os.path.join(src_dir, 'constellation_analysis.py'))
spectrum_analysis_old = load_module_from_path('spectrum_analysis_old', os.path.join(deprecated_dir, 'spectrum_analysis.py'))
spectrum_analysis_new = load_module_from_path('spectrum_analysis_new', os.path.join(src_dir, 'spectrum_analysis.py'))
single_frequency_laser_old = load_module_from_path('single_frequency_laser_old', os.path.join(deprecated_dir, 'single_frequency_laser.py'))
single_frequency_laser_new = load_module_from_path('single_frequency_laser_new', os.path.join(src_dir, 'single_frequency_laser.py'))


def test_edfa_regression():
    np.random.seed(42)
    t = np.linspace(0, 1e-9, 128)
    s_in = np.zeros((2, len(t)), dtype=complex)

    out_new = EDFA_new.EDFA(s_in, t, G_edfa_dB=10, NF_dB=5)
    np.random.seed(42)
    out_old = EDFA_old.EDFA(s_in, t, G_edfa_dB=10, NF_dB=5)

    np.testing.assert_allclose(out_new, out_old, atol=1e-12)


def test_hybrid_network_regression():
    A_sig = np.array([1+1j, 2+0j])
    A_LO = np.array([0.5-0.5j, 1-1j])

    out_new = HybridNetwork_new.HybridNetwork(A_sig, A_LO)
    out_old = HybridNetwork_old.HybridNetwork(A_sig, A_LO)

    assert len(out_new) == len(out_old)
    for new_val, old_val in zip(out_new, out_old):
        np.testing.assert_allclose(new_val, old_val, atol=1e-12)


def test_mzm_regression():
    A = 1 + 1j
    Vrf_upper = 2.0
    Vrf_lower = -1.0

    out_new = MZM_new.MZM(A, Vrf_upper, Vrf_lower)
    out_old = MZM_old.MZM(A, Vrf_upper, Vrf_lower)

    np.testing.assert_allclose(out_new, out_old, atol=1e-12)


def test_optical_filter_regression():
    t = np.linspace(0, 1e-9, 256)
    s = np.exp(1j * 2 * np.pi * 1e9 * t)

    out_new = optical_filter_new.optical_filter(s, t, f0=1e9, BW=1e8)
    out_old = optical_filter_old.optical_filter(s, t, f0=1e9, BW=1e8)

    np.testing.assert_allclose(out_new, out_old, atol=1e-12)


def test_photodiode_regression():
    pytest.importorskip('scipy')

    photodiode_old = load_module_from_path('photodiode_old', os.path.join(deprecated_dir, 'photodiode.py'))
    photodiode_new = load_module_from_path('photodiode_new', os.path.join(src_dir, 'photodiode.py'))

    t = np.linspace(0, 1e-9, 256)
    s = np.exp(1j * 2 * np.pi * 1e9 * t)

    out_new = photodiode_new.photodiode(s, t, R=0.9, BW=1e9, shot_noise=False, thermal_noise=False)
    out_old = photodiode_old.photodiode(s, t, R=0.9, BW=1e9, shot_noise=False, thermal_noise=False)

    np.testing.assert_allclose(out_new, out_old, atol=1e-12)


def test_single_frequency_laser_regression():
    np.random.seed(101)
    t = np.linspace(0, 1e-9, 128)

    out_new = single_frequency_laser_new.single_frequency_laser(t, P=0.01, Delta_nu=1e6, Freq_offset=1e6, X_fraction=0.5, phi_pol=0.3)
    np.random.seed(101)
    out_old = single_frequency_laser_old.single_frequency_laser(t, P=0.01, Delta_nu=1e6, Freq_offset=1e6, X_fraction=0.5, phi_pol=0.3)

    np.testing.assert_allclose(out_new, out_old, atol=1e-12)


def test_plotter_smoke_regression():
    t = np.linspace(0, 1e-9, 128)
    s = np.exp(1j * 2 * np.pi * 1e9 * t)

    constellation_analysis_new.constellation_analysis(s)
    constellation_analysis_old.constellation_analysis(s)
    spectrum_analysis_new.spectrum_analysis(s, t)
    spectrum_analysis_old.spectrum_analysis(s, t)
