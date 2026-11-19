"""Processamento digital do receptor coerente (docs/sistema_coerente.md, secao 11).

Decisoes registradas no pacote (docs/decisoes.md do pacote): o receptor conhece o comprimento e a
dispersao verdadeiros do enlace (D-7), e nesta fase usa os simbolos transmitidos conhecidos
("data-aided"), deixando o DSP cego do VPI (MMA, Stokes, busca cega de fase) como marco proprio.
Com os simbolos conhecidos a medida ve a fisica do enlace, sem a convergencia de um algoritmo
adaptativo no meio.

Duas partes:
- sem os simbolos: compensacao de dispersao (Agrawal, secao 8.7.3: no receptor coerente a
  dispersao e desfeita no dominio eletrico, porque a corrente guarda a fase do campo) e, opcional,
  a estimativa do desvio de frequencia pela quarta potencia;
- com os simbolos: instante dos simbolos e desvio de frequencia, matriz de polarizacao, fase da
  portadora e equalizador de minimos quadrados.

Nada aqui corrige o desbalanco de quadratura do transmissor: o equalizador e linear nos complexos,
como o do VPI, e a imagem conj(a) fica como ruido.
"""
from __future__ import annotations
import numpy as np

def dispersion_compensation_response(freq, beta2_length: float, frequency_shift: float = 0.0,
                                     beta3_length: float = 0.0):
    """H(f) que desfaz a dispersao acumulada `beta2_length` = soma de beta2*L [s^2] das fibras.

    A fibra multiplica o espectro do campo por exp[+(i/2)*beta2*(2*pi*f_g)^2*L] (kernels/ssfm.py),
    com f_g a frequencia do eixo da FFT da grade optica. O receptor mistura o sinal com o
    oscilador local, cujo envelope e exp(-i*2*pi*F*t): a componente que estava em f_g aparece em
    f = f_g + F. `frequency_shift` e esse F [Hz] (desvio do canal mais o desvio do oscilador
    local). Com F, a compensacao tira tambem o atraso de grupo do canal em relacao ao centro da
    banda.
    """
    omega = 2.0 * np.pi * (np.asarray(freq, dtype=float) - frequency_shift)
    # com inclinacao da dispersao, a fibra aplica tambem -(i/6)*beta3*w^3*L (kernels/ssfm.py)
    return np.exp(-0.5j * beta2_length * omega ** 2 + (1j / 6.0) * beta3_length * omega ** 3)

def estimate_frequency_offset(signal, sampling_rate: float, max_offset: float) -> float:
    """Desvio de frequencia [Hz] de um sinal QAM (n_pol, N), pela raia da quarta potencia.

    Para constelacoes com simetria de 90 graus, E[a^4] != 0: o sinal elevado a quarta tem uma raia
    em 4 vezes o desvio. O pico do espectro de x^4 (somado nas polarizacoes) e procurado em
    |f| <= 4*max_offset e refinado por uma parabola nos tres pontos em torno do pico. O sinal do
    resultado segue o eixo da FFT: o sinal gira como exp(+i*2*pi*desvio*t).
    """
    x = np.atleast_2d(np.asarray(signal, dtype=complex))
    n = x.shape[1]
    spectrum = np.sum(np.abs(np.fft.fft(x ** 4, axis=1)) ** 2, axis=0)
    freq = np.fft.fftfreq(n, 1.0 / sampling_rate)
    inside = np.abs(freq) <= 4.0 * max_offset
    k = int(np.flatnonzero(inside)[np.argmax(spectrum[inside])])
    three = spectrum[[(k - 1) % n, k, (k + 1) % n]]
    if np.any(three <= 0.0):                      # raia sozinha (sinal sem ruido): sem refinar
        return float(freq[k] / 4.0)
    left, mid, right = np.log(three)
    denom = left - 2.0 * mid + right
    shift = 0.0 if denom == 0.0 else 0.5 * (left - right) / denom
    return float((freq[k] + shift * sampling_rate / n) / 4.0)

def timing_and_frequency(signal, symbols, samples_per_symbol: int, first_symbol: int = 0,
                         last_symbol: int | None = None) -> tuple[int, float]:
    """(atraso d [amostras], rotacao de fase por simbolo [rad]) entre o sinal (n_pol, N) e os
    simbolos transmitidos (n_pol, n_sym), pelo produto diferencial.

    O produto de cada amostra pelo conjugado da amostra um simbolo antes, somado nas polarizacoes,
        D(t) = sum_p x_p(t)*conj(x_p(t - T)),
    nao depende da matriz de polarizacao (para uma matriz unitaria, x^H(t - T) x(t) = a^H a), tira
    a fase da portadora e transforma um desvio de frequencia df numa rotacao constante
    exp(i*2*pi*df*T). A correlacao circular de D com o trem dos produtos transmitidos
    sum_q a_q[k]*conj(a_q[k-1]) tem o pico no instante dos simbolos, e a fase do pico e a rotacao
    por simbolo. Vale para |df| < Rs/2, e a correlacao e coerente na janela inteira, com qualquer
    ruido de fase dos lasers.

    A raia da quarta potencia (`estimate_frequency_offset`) nao serve de base para o 16QAM: ela e
    fraca e o ruido de fase a alarga (medido em 2026-11-15: de 7 a 84 vezes a mediana do espectro,
    conforme a semente, e errada com 2^13 simbolos).
    """
    x = np.atleast_2d(np.asarray(signal, dtype=complex))
    a = np.atleast_2d(np.asarray(symbols, dtype=complex))
    sps, n = int(samples_per_symbol), x.shape[1]
    last = a.shape[1] if last_symbol is None else last_symbol
    received = np.sum(x * np.conj(np.roll(x, sps, axis=1)), axis=0)
    sent = np.sum(a[:, first_symbol + 1:last] * np.conj(a[:, first_symbol:last - 1]), axis=0)
    impulses = np.zeros(n, dtype=complex)
    impulses[(first_symbol + 1) * sps:last * sps:sps] = sent
    corr = np.fft.ifft(np.fft.fft(received) * np.conj(np.fft.fft(impulses)))
    d = int(np.argmax(np.abs(corr)))
    rotation = float(np.angle(corr[d]))
    return (d - n if d > n // 2 else d), rotation

def _moving_average(z: np.ndarray, length: int) -> np.ndarray:
    """Media movel centrada de `length` pontos (encolhe nas bordas)."""
    if length <= 1:
        return z
    kernel = np.ones(length)
    return np.convolve(z, kernel, mode="same") / np.convolve(np.ones(z.size), kernel, mode="same")

def _least_squares_matrix(received: np.ndarray, sent: np.ndarray) -> np.ndarray:
    """H (n_pol, n_pol) de minimos quadrados em received ~ H @ sent."""
    return (received @ sent.conj().T) @ np.linalg.inv(sent @ sent.conj().T)

def carrier_and_polarization(received, sent, block_length: int, phase_window: int):
    """Fase da portadora phi[k] e matriz de polarizacao H, com os simbolos conhecidos.

    Modelo: received[:, k] = exp(i*phi[k]) * H @ sent[:, k] + ruido, com H constante na janela (a
    polarizacao da fibra nao muda durante a simulacao) e phi comum as duas polarizacoes (o mesmo
    laser e o mesmo oscilador local).

    1. Em blocos de `block_length` simbolos, H_b de minimos quadrados. A fase de cada bloco e o
       angulo entre H_b e uma referencia acumulada. Interpolada entre os centros dos blocos.
    2. Com essa fase tirada, um H so, de minimos quadrados na janela inteira.
    3. Fase fina, simbolo a simbolo: angulo da media movel, em `phase_window` simbolos, de
       sum_pol y*conj(a), com y = H^-1 @ received ja sem a fase dos blocos.
    Devolve (phi (n,), H (n_pol, n_pol)).
    """
    r = np.atleast_2d(np.asarray(received, dtype=complex))
    a = np.atleast_2d(np.asarray(sent, dtype=complex))
    n = r.shape[1]
    block_length = max(4 * r.shape[0], min(int(block_length), n))
    starts = np.arange(0, n - block_length + 1, block_length)
    centers = starts + 0.5 * (block_length - 1)
    theta = np.zeros(starts.size)
    reference = None
    for b, s in enumerate(starts):
        h = _least_squares_matrix(r[:, s:s + block_length], a[:, s:s + block_length])
        if reference is None:
            reference = h.copy()
            continue
        theta[b] = np.angle(np.vdot(reference, h))
        reference += h * np.exp(-1j * theta[b])
    theta = np.unwrap(theta)
    phi = np.interp(np.arange(n), centers, theta)
    h = _least_squares_matrix(r * np.exp(-1j * phi), a)
    for _ in range(2):                # fase fina e H de novo, com a fase melhor
        y = np.linalg.solve(h, r * np.exp(-1j * phi))
        z = _moving_average(np.sum(y * np.conj(a), axis=0), phase_window)
        phi = phi + np.angle(z)
        h = _least_squares_matrix(r * np.exp(-1j * phi), a)
    return phi, h

def least_squares_equalizer(samples, sent, taps: int):
    """Equalizador linear de minimos quadrados, com `taps` coeficientes espacados de meio simbolo.

    samples: (n_pol, 2*n) amostras a 2 por simbolo, com a amostra 2*k no instante de decisao do
    simbolo k. sent: (n_pol, n). Cada polarizacao de saida usa as `taps` amostras em torno do
    instante de decisao, das duas polarizacoes de entrada (n_pol*taps coeficientes complexos,
    resolvidos pelas equacoes normais). As bordas sao tratadas como periodicas.
    Devolve os simbolos equalizados (n_pol, n).
    """
    u = np.atleast_2d(np.asarray(samples, dtype=complex))
    a = np.atleast_2d(np.asarray(sent, dtype=complex))
    n_pol, n = a.shape
    if taps < 1 or taps % 2 == 0:
        raise ValueError(f"taps = {taps}: use um numero impar de coeficientes")
    offsets = np.arange(taps) - taps // 2
    idx = (2 * np.arange(n)[:, None] + offsets[None, :]) % u.shape[1]
    regressors = np.concatenate([u[q][idx] for q in range(n_pol)], axis=1)      # (n, n_pol*taps)
    gram = regressors.conj().T @ regressors
    weights = np.linalg.solve(gram, regressors.conj().T @ a.T)                  # (n_pol*taps, n_pol)
    return (regressors @ weights).T

def data_aided_receiver(signal, symbols, samples_per_symbol: int, guard: int, *, taps: int = 7,
                        phase_window: int = 21, block_length: int = 64):
    """Simbolos equalizados (n_pol, n_medidos) de um sinal (n_pol, N) ja com a dispersao, o desvio
    de frequencia e o filtro casado tratados. Devolve (simbolos, atraso [amostras]).

    Cadeia: instante dos simbolos e desvio de frequencia (`timing_and_frequency`), fase e
    polarizacao (`carrier_and_polarization`),
    equalizador (`least_squares_equalizer`), de novo a fase fina, agora com os simbolos equalizados,
    e por fim o ganho complexo de cada polarizacao (tira o vies do equalizador de minimos
    quadrados). Os `guard` simbolos de cada borda da janela ficam fora.
    """
    x = np.atleast_2d(np.asarray(signal, dtype=complex))
    a_all = np.atleast_2d(np.asarray(symbols, dtype=complex))
    sps = int(samples_per_symbol)
    if sps % 2:
        raise ValueError(f"{sps} amostras por simbolo: o equalizador pede um numero par")
    first, last = guard, a_all.shape[1] - guard
    if last - first < 4 * max(block_length, taps):
        raise ValueError(f"guarda de {guard} simbolos por borda deixa {max(last - first, 0)} "
                         f"simbolos: poucos para o receptor")
    a = a_all[:, first:last]
    delay, rotation = timing_and_frequency(x, a_all, sps, first, last)
    x = x * np.exp(-1j * rotation * np.arange(x.shape[1]) / sps)      # tira o desvio de frequencia
    k = np.arange(first, last)
    on_time = x[:, (k * sps + delay) % x.shape[1]]
    phi, h = carrier_and_polarization(on_time, a, block_length, phase_window)
    # 2 amostras por simbolo, com a fase de cada simbolo aplicada as suas duas amostras
    idx2 = ((np.arange(2 * k.size) * (sps // 2)) + first * sps + delay) % x.shape[1]
    u = np.linalg.solve(h, x[:, idx2] * np.exp(-1j * np.repeat(phi, 2)))
    y = least_squares_equalizer(u, a, taps)
    y = y * np.exp(-1j * np.angle(_moving_average(np.sum(y * np.conj(a), axis=0), phase_window)))
    gains = np.sum(y * np.conj(a), axis=1) / np.sum(np.abs(a) ** 2, axis=1)
    return y / gains[:, None], delay

# --------------------------------------------------------------------------- DSP cego ----------
# Nada daqui para baixo usa os simbolos transmitidos (docs/sistema_coerente.md, secao 12).
#
# Fontes. Agrawal, secao 8.7.3: a matriz de Jones inversa do enlace pode ser construida "a partir
# do proprio sinal recebido, usando um algoritmo conhecido como algoritmo de modulo constante". O
# VPI (documento de parametros, secao 8) usa, nesta ordem: compensacao de dispersao, MIMO
# adaptativo (AlgTDE = MMA, NTaps = 7, 2 amostras por simbolo), separacao de polarizacao no espaco
# de Stokes, relogio, desvio de frequencia e busca cega de fase (NbTestAngles = 64). O professor
# (draft.py, `QAM_receiver`) faz a busca cega de fase em blocos de 16 simbolos com 21 angulos.
# Os algoritmos em si nao estao no Agrawal nem no documento do VPI. Referencias usuais:
#   modulo constante: D. Godard, IEEE Trans. Commun. 28, 1867 (1980);
#   dirigido por raio para 16QAM: I. Fatadin, D. Ives, S. Savory, J. Lightw. Technol. 27, 3042 (2009);
#   espaco de Stokes: B. Szafraniec, B. Nebendahl, T. Marshall, Opt. Express 18, 17928 (2010);
#   busca cega de fase: T. Pfau, S. Hoffmann, R. Noe, J. Lightw. Technol. 27, 989 (2009).

def stokes_polarization_demux(samples) -> np.ndarray:
    """Separa as duas polarizacoes de um sinal (2, N) pelo espaco de Stokes, sem os simbolos.

    Num sinal com multiplexacao em polarizacao, os vetores de Stokes das amostras ficam espalhados
    em torno de um plano que passa pela origem, e a normal desse plano aponta para os dois estados
    de polarizacao transmitidos. A normal e o autovetor de menor autovalor da matriz de covariancia
    dos vetores de Stokes. A matriz de Jones aplicada leva esse estado para a polarizacao X do
    receptor. Nao corrige a fase entre X e Y nem a PMD: isso fica com o equalizador.
    """
    x = np.asarray(samples, dtype=complex)
    if x.shape[0] != 2:
        raise ValueError("a separacao pelo espaco de Stokes pede um sinal de duas polarizacoes")
    ex, ey = x
    s = np.array([np.abs(ex) ** 2 - np.abs(ey) ** 2, 2.0 * np.real(ex * np.conj(ey)),
                  -2.0 * np.imag(ex * np.conj(ey))])
    _, vectors = np.linalg.eigh(s @ s.T / s.shape[1])
    n1, n2, n3 = vectors[:, 0]                    # normal do plano (menor autovalor)
    # vetor de Jones cujo vetor de Stokes e a normal: (cos(chi), sin(chi)*exp(i*delta))
    chi = 0.5 * np.arccos(np.clip(n1, -1.0, 1.0))
    delta = np.arctan2(-n3, n2)
    a, b = np.cos(chi), np.sin(chi) * np.exp(1j * delta)
    jones = np.array([[np.conj(a), np.conj(b)], [-b, a]])      # unitaria: leva (a, b) para (1, 0)
    return jones @ x

def modulus_radii(points) -> tuple[float, np.ndarray]:
    """(R2 do modulo constante, raios ao quadrado distintos) de uma constelacao de potencia media
    unitaria. R2 = E|a|^4/E|a|^2 (Godard). Para o 16QAM os raios ao quadrado sao 0,2, 1,0 e 1,8,
    os 2, 10 e 18 de `RadiiSq` do VPI divididos pela potencia media 10."""
    c = np.asarray(points, dtype=complex)
    c = c / np.sqrt(np.mean(np.abs(c) ** 2))
    p = np.abs(c) ** 2
    return float(np.mean(p ** 2) / np.mean(p)), np.unique(np.round(p, 9))

def _regressors(samples: np.ndarray, taps: int) -> np.ndarray:
    """(n, n_pol*taps): para cada simbolo k, as `taps` amostras em torno da amostra 2*k, de cada
    polarizacao. Bordas periodicas."""
    n = samples.shape[1] // 2
    idx = (2 * np.arange(n)[:, None] + (np.arange(taps) - taps // 2)[None, :]) % samples.shape[1]
    return np.concatenate([row[idx] for row in samples], axis=1)

def _adapt(u: np.ndarray, w: np.ndarray, step: float, radii2: np.ndarray, passes: int) -> np.ndarray:
    """Gradiente estocastico, simbolo a simbolo: w <- w + step*y*(R2 - |y|^2)*conj(u), com R2 o
    raio ao quadrado de `radii2` mais proximo de |y|^2 (um raio so: modulo constante)."""
    single = radii2.size == 1
    for _ in range(passes):
        for k in range(u.shape[0]):
            uk = u[k]
            y = w @ uk
            p = y.real ** 2 + y.imag ** 2
            r = radii2[0] if single else radii2[np.argmin(np.abs(radii2[:, None] - p[None, :]), axis=0)]
            w += (step * (r - p) * y)[:, None] * np.conj(uk)[None, :]
    return w

def adaptive_equalizer(samples, points, taps: int = 7, step: float = 1e-3, cma_passes: int = 4,
                       rde_passes: int = 0) -> np.ndarray:
    """Equalizador adaptativo cego, n_pol x n_pol, com `taps` coeficientes a 2 amostras por simbolo.

    samples: (n_pol, 2*n), a 2 amostras por simbolo, em qualquer fase de amostragem (o equalizador
    fracionario absorve o atraso). points: a constelacao. Devolve (n_pol, n) simbolos com potencia
    media unitaria, ainda com a fase da portadora e o desvio de frequencia, que o criterio do
    modulo nao enxerga.

    1. Cada polarizacao e normalizada para potencia unitaria. Os coeficientes comecam com 1 no
       centro (como `hxx` e `hyy` do VPI).
    2. `cma_passes` passagens pelo bloco com o modulo constante. Se as duas saidas convergirem para
       a mesma polarizacao (matriz dos coeficientes quase singular), a saida Y e refeita ortogonal
       a X e o modulo constante roda de novo.
    3. `rde_passes` passagens com o criterio dirigido por raio (o `MMA` com `RadiiSq` do VPI), que
       so muda algo em constelacoes com mais de um raio. Opcional: com pouca SNR os raios sao
       decididos errado, e quem refina os coeficientes e o ajuste por decisao de
       `blind_receiver`.
    4. Os coeficientes finais sao aplicados ao bloco inteiro: o canal e fixo na janela.
    `step` vale para sinais de potencia unitaria. O 5e-4 do VPI nao e comparavel: a normalizacao
    de la nao esta no documento.
    """
    x = np.atleast_2d(np.asarray(samples, dtype=complex))
    n_pol = x.shape[0]
    if taps < 1 or taps % 2 == 0:
        raise ValueError(f"taps = {taps}: use um numero impar de coeficientes")
    x = x / np.sqrt(np.mean(np.abs(x) ** 2, axis=1, keepdims=True))
    u = _regressors(x, taps)
    r2, radii2 = modulus_radii(points)
    w = np.zeros((n_pol, n_pol * taps), dtype=complex)
    for p in range(n_pol):
        w[p, p * taps + taps // 2] = 1.0
    cma = np.array([r2])
    w = _adapt(u, w, step, cma, cma_passes)
    if n_pol == 2:
        # matriz de Jones equivalente: soma dos coeficientes de cada entrada
        jones = w.reshape(2, 2, taps).sum(axis=2)
        singular = abs(np.linalg.det(jones)) < 0.3 * np.linalg.norm(jones[0]) * np.linalg.norm(jones[1])
        if singular:
            wx = w[0].reshape(2, taps)
            w[1] = np.concatenate([-np.conj(wx[1, ::-1]), np.conj(wx[0, ::-1])])
            w = _adapt(u, w, step, cma, cma_passes)
    if radii2.size > 1 and rde_passes > 0:
        w = _adapt(u, w, step, radii2, rde_passes)
    y = (u @ w.T).T
    return y / np.sqrt(np.mean(np.abs(y) ** 2, axis=1, keepdims=True))

def blind_phase_search(symbols, modulation: str, test_angles: int = 64, window: int = 21) -> np.ndarray:
    """Fase da portadora phi[k] de uma sequencia de simbolos (n,), pela busca cega de fase. Com
    varias linhas (n_pol, n), uma fase so para todas: as distancias das linhas sao somadas.

    Para cada um dos `test_angles` angulos entre 0 e 2*pi/S (S = simetria de rotacao da
    constelacao), gira os simbolos, decide pelo ponto mais proximo e soma a distancia ao quadrado
    numa janela de `window` simbolos. Em cada simbolo vale o angulo de menor soma. A fase sai
    desdobrada com passo 2*pi/S, e por isso tem uma ambiguidade constante de um multiplo de
    2*pi/S, que so os dados transmitidos (ou um codigo diferencial) resolvem.
    Os simbolos devem estar na escala da constelacao de kernels/mapping.py.
    """
    from .mapping import ROTATIONAL_SYMMETRY, nearest_points
    y = np.atleast_2d(np.asarray(symbols, dtype=complex))
    sector = 2.0 * np.pi / ROTATIONAL_SYMMETRY[modulation]
    angles = np.arange(test_angles) * sector / test_angles
    kernel = np.ones(min(int(window), y.shape[1]))
    cost = np.empty((test_angles, y.shape[1]))
    for b, angle in enumerate(angles):
        r = y * np.exp(-1j * angle)
        distance = np.sum(np.abs(r - nearest_points(r, modulation)) ** 2, axis=0)
        cost[b] = np.convolve(distance, kernel, mode="same")
    phase = angles[np.argmin(cost, axis=0)]
    return np.unwrap(phase, period=sector)

def decision_scale(symbols, modulation: str, iterations: int = 2) -> np.ndarray:
    """Leva uma sequencia de simbolos (n,) para a escala da constelacao sem os dados: divide pelo
    ganho de minimos quadrados contra as proprias decisoes, `iterations` vezes."""
    from .mapping import nearest_points
    z = np.asarray(symbols, dtype=complex)
    for _ in range(iterations):
        d = nearest_points(z, modulation)
        z = z / abs(np.sum(z * np.conj(d)) / np.sum(np.abs(d) ** 2))
    return z

def carrier_recovery(symbols, modulation: str, test_angles: int = 64, window: int = 21):
    """Tira o desvio de frequencia e a fase da portadora de simbolos (n_pol, n), sem os dados.

    O desvio sai da raia da quarta potencia dos simbolos ja equalizados, somada nas polarizacoes
    (`estimate_frequency_offset`, valida ate 1/8 da taxa de simbolos). Depois a busca cega de fase.

    Com duas polarizacoes a busca e conjunta, como no VPI (`CombinedPol = Yes`): a fase dos lasers
    e a mesma nas duas, entao uma busca por polarizacao da a diferenca constante entre elas (a
    menos de 2*pi/S, o que basta), e uma segunda busca usa as duas juntas, com o dobro de
    amostras por janela. Medido em 2026-11-15 (16QAM, 12,9 dB de SNR, janela de 81 simbolos): a
    busca por polarizacao saltou 15 vezes numa polarizacao e nenhuma na outra.

    Devolve (simbolos na escala da constelacao, fase total tirada de cada simbolo (n_pol, n)
    [rad], desvio de frequencia [ciclos por simbolo]).
    """
    from .mapping import ROTATIONAL_SYMMETRY, constellation
    y = np.atleast_2d(np.asarray(symbols, dtype=complex))
    rms = np.sqrt(np.mean(np.abs(constellation(modulation)) ** 2))
    offset = estimate_frequency_offset(y, 1.0, 1.0 / 8.0)
    ramp = 2.0 * np.pi * offset * np.arange(y.shape[1])
    rows = [decision_scale(row * np.exp(-1j * ramp) / np.sqrt(np.mean(np.abs(row) ** 2)) * rms,
                           modulation) for row in y]
    phis = [blind_phase_search(row, modulation, test_angles, window) for row in rows]
    if len(rows) == 2:
        symmetry = ROTATIONAL_SYMMETRY[modulation]
        delta = np.angle(np.mean(np.exp(1j * symmetry * (phis[1] - phis[0])))) / symmetry
        joint = blind_phase_search(np.array([rows[0], rows[1] * np.exp(-1j * delta)]), modulation,
                                   test_angles, window)
        phis = [joint, joint + delta]
    out = [decision_scale(row * np.exp(-1j * phi), modulation) for row, phi in zip(rows, phis)]
    return np.array(out), np.array([ramp + phi for phi in phis]), float(offset)

def blind_receiver(samples, modulation: str, *, taps: int = 7, step: float = 1e-3,
                   cma_passes: int = 4, rde_passes: int = 0, dd_iterations: int = 4,
                   test_angles: int = 64, window: int = 21, stokes: bool = True):
    """Simbolos (n_pol, n) de um sinal a 2 amostras por simbolo (n_pol, 2*n), sem os dados
    transmitidos. Devolve (simbolos, desvio de frequencia [ciclos por simbolo]).

    1. Com duas polarizacoes, separacao pelo espaco de Stokes (`stokes_polarization_demux`).
    2. Equalizador adaptativo: modulo constante e depois dirigido por raio (`adaptive_equalizer`).
    3. Desvio de frequencia e busca cega de fase (`carrier_recovery`).
    4. `dd_iterations` vezes, ajuste dirigido por decisao: os coeficientes do equalizador sao
       refeitos por minimos quadrados no bloco, tendo como alvo as proprias decisoes levadas de
       volta a fase recebida, e a portadora e recuperada de novo. E o equivalente em bloco do LMS
       dirigido por decisao (o VPI lista `Mode = Decision-Directed`). Sem este passo o criterio do
       modulo fica 0,7 dB abaixo da solucao de minimo erro (medido em 2026-11-15, 16QAM a 16 dB).
       Com 13 dB de SNR foram precisos 4 ajustes para as duas polarizacoes convergirem.

    O que fica sem resolver, porque so os dados transmitidos resolvem: qual saida e qual
    polarizacao, o atraso em simbolos e um giro de multiplo de 2*pi/S. Abaixo de uma SNR a busca
    cega de fase salta de um giro para outro no meio do bloco: a janela maior reduz os saltos.
    """
    from .mapping import constellation, nearest_points
    x = np.atleast_2d(np.asarray(samples, dtype=complex))
    if x.shape[1] % 2:
        raise ValueError("o receptor cego pede 2 amostras por simbolo (numero par de amostras)")
    if stokes and x.shape[0] == 2:
        x = stokes_polarization_demux(x)
    points = constellation(modulation)
    y = adaptive_equalizer(x, points, taps, step, cma_passes, rde_passes)
    z, phases, offset = carrier_recovery(y, modulation, test_angles, window)
    if dd_iterations > 0:
        u = _regressors(x / np.sqrt(np.mean(np.abs(x) ** 2, axis=1, keepdims=True)), taps)
        gram = u.conj().T @ u
        for _ in range(dd_iterations):
            target = nearest_points(z, modulation) * np.exp(1j * phases)
            weights = np.linalg.solve(gram, u.conj().T @ target.T)
            z, phases, offset = carrier_recovery((u @ weights).T, modulation, test_angles, window)
    return z, offset
