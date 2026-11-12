"""Amplificador optico concentrado (EDFA): ganho + ruido de ASE (docs/amplificadores.md).

Agrawal, Sistemas de Comunicacao por Fibra Optica, secoes 7.2.3, 7.4 e 7.5:
    S_ASE = n_sp*h*nu_0*(G - 1)                   por polarizacao      (Eq. 7.2.11)
    F_n   = 2*n_sp*(1 - 1/G) + 1/G  ~  2*n_sp                         (Eq. 7.2.15)
    P_ASE = 2*S_ASE*delta_nu  (as duas polarizacoes: a ASE e despolarizada) (Eq. 7.4.1, 7.5.10)
Modelo do codigo do professor (notebook 4_EDFAsAndReceivers, draft.py `EDFA`):
    s_out = sqrt(G)*s_in + n_ASE(t), n_ASE ruido branco gaussiano na banda inteira da simulacao.
O ruido aqui e gaussiano circular (correcao C-2 da arquitetura; o professor usa np.random.rand).
"""
from __future__ import annotations
import numpy as np
from ..core.units import H_PLANCK

def ase_psd(n_sp, gain, frequency):
    """Densidade espectral da ASE por polarizacao, S_ASE [W/Hz] (Eq. 7.2.11)."""
    return n_sp * H_PLANCK * frequency * (gain - 1.0)

def nsp_from_noise_figure(noise_figure, gain):
    """n_sp que da a figura de ruido F_n (linear) no ganho G (linear): Eq. 7.2.15 invertida,
    n_sp = (G*F_n - 1)/(2*(G - 1)). Para G >> 1 tende a F_n/2, a aproximacao do professor."""
    return (gain * noise_figure - 1.0) / (2.0 * (gain - 1.0))

def amplify(field, gain, psd, sampling_rate, rng=None):
    """Campo (n_pol, N) amplificado: sqrt(G)*campo + ASE em cada polarizacao.

    A ASE de cada polarizacao e ruido gaussiano circular de densidade `psd` [W/Hz] na banda
    inteira da grade (largura `sampling_rate`): potencia psd*sampling_rate por polarizacao,
    metade em cada quadratura. Sem `rng`, so o ganho (amplificador sem ruido).

    Um campo de uma polarizacao (1, N) sai como (2, N), com Y = 0 antes do ruido: a ASE e
    despolarizada e entra nas duas (decisao de 2026-10-22).
    """
    field = np.asarray(field, dtype=complex)
    if field.shape[0] == 1:
        field = np.vstack([field, np.zeros_like(field)])
    out = np.sqrt(gain) * field
    if rng is not None and psd > 0.0:
        sigma = np.sqrt(psd * sampling_rate / 2.0)
        out = out + sigma * (rng.standard_normal(field.shape) + 1j * rng.standard_normal(field.shape))
    return out
