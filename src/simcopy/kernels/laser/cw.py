import numpy as np

def cw_field(n_samples, dt, power, linewidth, freq_offset, rng) -> np.ndarray:
    """
    Campo (N,) complexo: sqrt(P)·exp[i(-2π·Δf·t + φ(t))], φ passeio aleatório com
    incrementos N(0, 2π·Δν·dt). Agrawal Eq. 3.3.33–3.3.37 (docs/laser.md, seção 2).

    Δf é o desvio de frequência óptica (positivo = azul). Na convenção do Agrawal (portadora
    exp(-i·ω0·t), decisão de 2026-10-22), uma componente em +Δf tem envelope exp(-i·2π·Δf·t).

    A fase inicial φ(0) é sorteada uniforme em [0, 2π) (decisão de 2026-10-08); cada
    amostra seguinte soma um incremento gaussiano de variância 2π·Δν·dt.

    Args:
      n_samples: número de amostras
      dt: passo de tempo (s)
      power: potência média (W)
      linewidth: largura de linha (Hz)
      freq_offset: desvio de frequência óptica (Hz), positivo = azul
      rng: gerador numpy.random.Generator (no simulador, vem de ctx.rng_for(...))
    Outputs:
      cw_wave: campo contínuo (N,) complexo
    """
    phi0 = rng.uniform(0.0, 2 * np.pi)
    dphi = rng.normal(0, np.sqrt(2 * np.pi * linewidth * dt), n_samples - 1)
    phi = phi0 + np.concatenate(([0.0], np.cumsum(dphi)))
    cw_wave = np.sqrt(power) * np.exp(1j * (-2 * np.pi * freq_offset * dt * np.arange(n_samples) + phi))
    return cw_wave
