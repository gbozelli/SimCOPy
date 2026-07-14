import importlib.util
import os
import sys
import numpy as np

repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
src_dir = os.path.join(repo_root, 'src')
deprecated_dir = os.path.join(repo_root, 'src', 'deprecated')


def load_module_from_path(name, path, extra_paths=None):
    extra_paths = extra_paths or []
    old_sys_path = list(sys.path)
    try:
        for p in extra_paths:
            if p not in sys.path:
                sys.path.insert(0, p)
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path[:] = old_sys_path

qam_transmitter_old = load_module_from_path(
    'qam_transmitter_old', os.path.join(deprecated_dir, 'QAM_transmitter.py'),
    extra_paths=[deprecated_dir]
)
qam_transmitter_new = load_module_from_path(
    'qam_transmitter_new', os.path.join(src_dir, 'qam_transmitter.py'),
    extra_paths=[src_dir]
)


def test_qam_transmitter_regression():
    np.random.seed(123)
    M = 16
    N_inf = 256
    sync_seed = 7

    s_cont_old, t_old, s_b_old, const_tx_old = qam_transmitter_old.QAM_transmitter(
        N_inf=N_inf, M=M, sync_seed=sync_seed, plot_flag=False
    )

    np.random.seed(123)
    s_cont_new, t_new, s_b_new, const_tx_new = qam_transmitter_new.QAM_transmitter(
        N_inf=N_inf, M=M, sync_seed=sync_seed, plot_flag=False
    )

    np.testing.assert_array_equal(s_b_new, s_b_old)
    np.testing.assert_array_equal(const_tx_new, const_tx_old)
    np.testing.assert_array_equal(t_new, t_old)
    np.testing.assert_array_equal(s_cont_new, s_cont_old)
