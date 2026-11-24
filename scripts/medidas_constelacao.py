"""Medidas nas constelacoes depois do DSP, com o mesmo codigo para o VPI e para o simcopy
(docs/plano_replicar_vpi.md, etapa A; docs/replicacao_vpi_diario.md).

Entrada: simbolos recebidos x e transmitidos a, das duas polarizacoes. Primeiro o ajuste de
docs/comparacao_vpi.md, x = g*a + h*conj(a) + c + e. Tudo o que segue mede o erro e/g, o "resto":

- radial e tangencial por anel: var = aditivo/2 + inclinacao*|a|^2. A inclinacao tangencial e
  ruido de fase [rad^2], a radial e ruido de amplitude. A media por anel e a deformacao da
  constelacao (compressao e rotacao que dependem do raio).
- memoria linear: quanto do resto e funcao linear dos simbolos vizinhos da mesma polarizacao
  (interferencia entre simbolos) e dos simbolos da outra polarizacao (diafonia).
- autocovariancia do erro de fase, dentro de uma polarizacao e entre as duas. Para um estimador de
  fase por media movel de W simbolos, o erro medido contra o proprio simbolo tem autocovariancia
  negativa ate o atraso (W-1)/2, positiva dai ate W-1 e nula depois: a posicao do salto da a janela.
- saltos de fase: giro de cada bloco de simbolos em relacao aos dados.
- SER medida contra a de ruido gaussiano com a mesma SNR, curtose, erro ao longo do bloco.

    python scripts/medidas_constelacao.py sintetico
    python scripts/medidas_constelacao.py vpi "<pasta>/dados4canais/extraido" --lengths 125 --powers 13 -o saida.json
    python scripts/medidas_constelacao.py simbolos resultados/replica_vpi/simbolos/*.npz -o saida.json
    python scripts/medidas_constelacao.py tabela vpi.json simcopy.json

Os dados do VPI nao ficam no repositorio: a pasta e passada na linha de comando.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import erfc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

MAX_LAG = 128
SLIP_BLOCK = 128

# ------------------------------------------------------------------ ajuste e resto
def widely_linear(x, a):
    """x = g*a + h*conj(a) + c + e. Devolve (g, h, c, e/g)."""
    m = np.column_stack([a, np.conj(a), np.ones_like(a)])
    (g, h, c), *_ = np.linalg.lstsq(m, x, rcond=None)
    return g, h, c, (x - m @ np.array([g, h, c])) / g

# ------------------------------------------------------------------ radial e tangencial
def radial_tangential(eps, a) -> dict:
    """Decompoe o resto em parte radial e tangencial por anel e ajusta var = i + s*|a|^2 (pesos:
    numero de simbolos do anel). Tudo em unidades de 1/SNR (potencia do erro sobre a do sinal)."""
    p = np.abs(a) ** 2
    ps = float(np.mean(p))
    u = a / np.abs(a)
    r, t = np.real(eps * np.conj(u)), np.imag(eps * np.conj(u))
    rings = np.unique(np.round(p, 6))
    count, var_r, var_t, mean_r, mean_t = [], [], [], [], []
    for ring in rings:
        sel = np.abs(p - ring) < 1e-6 * max(1.0, ring)
        count.append(int(sel.sum()))
        mean_r.append(float(r[sel].mean()))
        mean_t.append(float(t[sel].mean()))
        var_r.append(float(r[sel].var()))
        var_t.append(float(t[sel].var()))
    w = np.array(count, dtype=float)
    out = {"rings": (rings / ps).tolist(), "count": count,
           "var_r": (np.array(var_r) / ps).tolist(), "var_t": (np.array(var_t) / ps).tolist(),
           # deformacao: media do erro em cada anel, como fracao do raio do anel
           "compression": (np.array(mean_r) / np.sqrt(rings)).tolist(),
           "rotation_rad": (np.array(mean_t) / np.sqrt(rings)).tolist(),
           "warp": float(np.sum(w * (np.array(mean_r) ** 2 + np.array(mean_t) ** 2)) / w.sum() / ps),
           "total": float(np.mean(np.abs(eps) ** 2) / ps)}
    if rings.size >= 2:
        design = np.column_stack([np.ones_like(rings), rings])
        sw = np.sqrt(w)[:, None]
        (i_r, s_r), *_ = np.linalg.lstsq(design * sw, np.array(var_r) * sw[:, 0], rcond=None)
        (i_t, s_t), *_ = np.linalg.lstsq(design * sw, np.array(var_t) * sw[:, 0], rcond=None)
        out.update(additive=float((i_r + i_t) / ps), phase=float(s_t), amplitude=float(s_r),
                   additive_r=float(2 * i_r / ps), additive_t=float(2 * i_t / ps))
    return out

# ------------------------------------------------------------------ memoria linear
def linear_memory(x, a, b=None, taps: int = 5, edge: int = 0) -> dict:
    """Parte do resto que e funcao linear dos vizinhos. Regressao em tres passos, sobre os
    simbolos k de `edge` a n - edge (bordas circulares se edge = 0):
      0: a_k, conj(a_k), 1                         -> resto (o de docs/comparacao_vpi.md)
      1: + a_{k+j}, conj(a_{k+j}), j != 0          -> o que sai e interferencia entre simbolos
      2: + b_{k+j}, conj(b_{k+j}), todos os j      -> o que sai e diafonia da outra polarizacao
    Devolve as tres potencias sobre a do sinal e os coeficientes dos vizinhos."""
    n = a.size
    k = np.arange(edge, n - edge) if edge else np.arange(n)
    lags = [j for j in range(-taps, taps + 1) if j != 0]
    col = lambda v, j: v[(k + j) % n]
    base = [a[k], np.conj(a[k]), np.ones(k.size, dtype=complex)]
    own = [f(col(a, j)) for j in lags for f in (lambda v: v, np.conj)]
    cross = [] if b is None else [f(col(b, j)) for j in range(-taps, taps + 1) for f in (lambda v: v, np.conj)]
    y = x[k]
    def fit(columns):
        m = np.column_stack(columns)
        coef, *_ = np.linalg.lstsq(m, y, rcond=None)
        return coef, float(np.mean(np.abs(y - m @ coef) ** 2))
    c0, p0 = fit(base)
    c1, p1 = fit(base + own)
    c2, p2 = fit(base + own + cross) if cross else (c1, p1)
    signal = float(np.mean(np.abs(c0[0] * a[k]) ** 2))
    g = c1[0]
    own_taps = {int(j): [float(abs(c1[3 + 2 * i] / g)), float(abs(c1[4 + 2 * i] / g))]
                for i, j in enumerate(lags)}
    # coeficientes complexos, relativos ao principal: [linear, conjugado] de cada vizinho. O
    # conjugado e dividido por conj(g): assim nao depende da rotacao comum da constelacao.
    lin = {j: c1[3 + 2 * i] / g for i, j in enumerate(lags)}
    conj = {j: c1[4 + 2 * i] / np.conj(g) for i, j in enumerate(lags)}
    complex_taps = {int(j): [[float(lin[j].real), float(lin[j].imag)], [float(conj[j].real), float(conj[j].imag)]]
                    for j in lags}
    # parte simetrica e antissimetrica do nucleo conjugado: um skew entre I e Q da um nucleo
    # antissimetrico (derivada do pulso); uma diferenca de resposta entre I e Q, um simetrico
    sym = float(sum(abs(conj[j] + conj[-j]) ** 2 for j in range(1, taps + 1)) / 2.0)
    anti = float(sum(abs(conj[j] - conj[-j]) ** 2 for j in range(1, taps + 1)) / 2.0)
    lin_sym = float(sum(abs(lin[j] + lin[-j]) ** 2 for j in range(1, taps + 1)) / 2.0)
    lin_anti = float(sum(abs(lin[j] - lin[-j]) ** 2 for j in range(1, taps + 1)) / 2.0)
    return {"rest": p0 / signal, "isi": (p0 - p1) / signal, "crosstalk": (p1 - p2) / signal,
            "unexplained": p2 / signal, "n": int(k.size), "regressors": len(base + own + cross),
            "own_taps": own_taps, "complex_taps": complex_taps,
            "conj_symmetric": sym, "conj_antisymmetric": anti,
            "linear_symmetric": lin_sym, "linear_antisymmetric": lin_anti}

# ------------------------------------------------------------------ erro de fase e de amplitude
def phase_amplitude_series(eps, a):
    """v = Im(e*conj(a)) e w = Re(e*conj(a)), divididos por E|a|^2: para ruido de fase comum phi e
    de amplitude m, v ~ phi*|a|^2/E|a|^2 e w ~ m*|a|^2/E|a|^2, mais ruido aditivo branco."""
    ps = np.mean(np.abs(a) ** 2)
    z = eps * np.conj(a) / ps
    return z.imag, z.real

def covariance(u, v, max_lag: int = MAX_LAG, circular: bool = True) -> np.ndarray:
    """E[u_k * v_{k+l}] para l = 0..max_lag (u e v ja sem a media)."""
    n = u.size
    if circular:
        return np.array([np.mean(u * np.roll(v, -l)) for l in range(max_lag + 1)])
    return np.array([np.mean(u[:n - l] * v[l:]) for l in range(max_lag + 1)])

def window_from_covariance(c: np.ndarray, noise: float) -> dict:
    """Janela W de um estimador por media movel, lida da autocovariancia do erro de fase c[l]
    (l = 0..): negativa ate h = (W-1)/2 e positiva logo depois, com um degrau entre h e h + 1.
    Devolve o atraso do maior degrau positivo (a partir do atraso 2), se ele passar de 4 vezes o
    desvio padrao `noise` de cada c[l]. O cruzamento de zero nao serve: o ruido de fase nao linear,
    que tem memoria de alguns simbolos, soma uma parte positiva e o desloca."""
    step = np.diff(c)
    l = int(np.argmax(step[2:])) + 2
    if step[l] < 4.0 * noise * np.sqrt(2.0):
        return {"half_window": None, "window": None, "step": float(step[l])}
    return {"half_window": l, "window": 2 * l + 1, "step": float(step[l])}

# ------------------------------------------------------------------ saltos de fase
def slips(x, a, block: int = SLIP_BLOCK, sector: float = np.pi / 2) -> dict:
    """Giro de cada bloco em relacao aos dados (multiplo de `sector`), numero de mudancas de giro,
    e blocos sem correlacao com os dados (|correlacao normalizada| < 0,5)."""
    n = (a.size // block) * block
    c = np.sum((x[:n] * np.conj(a[:n])).reshape(-1, block), axis=1)
    norm = np.sqrt(np.sum(np.abs(x[:n].reshape(-1, block)) ** 2, axis=1) *
                   np.sum(np.abs(a[:n].reshape(-1, block)) ** 2, axis=1))
    rho = np.abs(c) / norm
    turns = np.round(np.angle(c) / sector).astype(int) % int(round(2 * np.pi / sector))
    good = rho >= 0.5
    changes = int(np.count_nonzero(np.diff(turns[good]))) if good.sum() > 1 else 0
    hist = np.bincount(turns[good], minlength=int(round(2 * np.pi / sector)))
    # giro de cada bloco sem arredondar, em relacao ao giro mais comum: nos arquivos em que o VPI
    # se perde a constelacao fica girada perto de 45 graus, que o arredondamento conta errado
    main = np.bincount(turns).argmax() * sector
    angle = np.degrees(np.angle(c * np.exp(-1j * main)))
    rotated = np.abs(angle) > 22.5
    return {"blocks": int(c.size), "lost_blocks": int((~good).sum()), "slips": changes,
            "rotated_fraction": float(rotated.mean()),
            "rotated_episodes": int(np.count_nonzero(np.diff(rotated.astype(int)) == 1) + int(rotated[0])),
            "angle_abs_max_deg": float(np.abs(angle).max()),
            "turn_histogram": hist.tolist(), "mean_correlation": float(rho.mean()),
            "min_correlation": float(rho.min()),
            "first_bad_block": int(np.flatnonzero((turns != np.bincount(turns).argmax()) | ~good)[0])
            if ((turns != np.bincount(turns).argmax()) | ~good).any() else None}

def best_alignment(x, a, alternatives: dict, max_lag: int = 32) -> dict:
    """Procura, entre `alternatives` (nome -> sequencia candidata a ser a transmitida) e atrasos de
    -max_lag a max_lag, a de maior soma de |correlacao por bloco|^2 (nao depende do giro)."""
    n = (a.size // SLIP_BLOCK) * SLIP_BLOCK
    best = None
    for name, cand in alternatives.items():
        for lag in range(-max_lag, max_lag + 1):
            c = np.sum((np.roll(x, -lag)[:n] * np.conj(cand[:n])).reshape(-1, SLIP_BLOCK), axis=1)
            score = float(np.sum(np.abs(c) ** 2))
            if best is None or score > best[0]:
                best = (score, name, lag)
    reference = float(np.sum(np.abs(a[:n]) ** 2) * np.sum(np.abs(x[:n]) ** 2) / (n // SLIP_BLOCK))
    return {"candidate": best[1], "lag": best[2], "score": best[0] / reference}

# ------------------------------------------------------------------ outras
def ser_gaussian_16qam(snr: float) -> float:
    """SER do 16QAM com ruido gaussiano circular, para a SNR linear (Es/N0)."""
    p = 1.5 * 0.5 * erfc(np.sqrt(snr / 10.0))       # erro numa quadratura: 2*(1 - 1/4)*Q(sqrt(snr/5))
    return float(1.0 - (1.0 - p) ** 2)

def decide_square(x, levels) -> np.ndarray:
    lv = np.asarray(levels, dtype=float)
    q = lambda v: lv[np.argmin(np.abs(v[:, None] - lv[None, :]), axis=1)]
    return q(x.real) + 1j * q(x.imag)

def excess_kurtosis(v) -> float:
    v = v - v.mean()
    return float(np.mean(v ** 4) / np.mean(v ** 2) ** 2 - 3.0)

def along_block(eps, segments: int = 32) -> dict:
    n = (eps.size // segments) * segments
    p = np.mean(np.abs(eps[:n]) ** 2 .reshape(segments, -1), axis=1) if False else \
        np.mean((np.abs(eps[:n]) ** 2).reshape(segments, -1), axis=1)
    p = p / p.mean()
    return {"first": float(p[0]), "last": float(p[-1]), "max": float(p.max()), "min": float(p.min()),
            "std": float(p.std())}

# ------------------------------------------------------------------ todas as medidas de um ponto
def measure(received, sent, circular: bool, edge: int = 0) -> dict:
    """received, sent: (2, n). Devolve as medidas de cada polarizacao e as cruzadas."""
    x, a = np.asarray(received, dtype=complex), np.asarray(sent, dtype=complex)
    levels = np.unique(np.round(a.real, 6))
    pols, series = [], []
    for p in range(a.shape[0]):
        g, h, c, eps = widely_linear(x[p], a[p])
        ps = float(np.mean(np.abs(a[p]) ** 2))
        total = float(np.mean(np.abs(x[p] / g - a[p]) ** 2) / ps)
        v, w = phase_amplitude_series(eps, a[p])
        series.append((v - v.mean(), w - w.mean()))
        rt = radial_tangential(eps, a[p])
        mem = linear_memory(x[p], a[p], a[1 - p] if a.shape[0] == 2 else None, edge=edge)
        ser = float(np.mean(np.abs(decide_square(x[p] / g, levels) - a[p]) > 1e-3 * np.abs(levels).max()))
        cv = covariance(series[-1][0], series[-1][0], circular=circular)
        cw = covariance(series[-1][1], series[-1][1], circular=circular)
        noise = float(cv[0] / np.sqrt(a.shape[1]))
        pols.append({
            "image": float(abs(h / g)), "offset": float(abs(c / g) / np.sqrt(ps)),
            "inv_snr_total": total, "inv_snr_rest": rt["total"], "ser": ser,
            "ser_gaussian": ser_gaussian_16qam(1.0 / total) if levels.size == 4 else None,
            "radial_tangential": rt, "memory": mem,
            "phase_cov": cv.tolist(), "amplitude_cov": cw.tolist(),
            "window": window_from_covariance(cv, noise),
            "kurtosis_re": excess_kurtosis(eps.real), "kurtosis_im": excess_kurtosis(eps.imag),
            "along_block": along_block(eps), "slips": slips(x[p] / g, a[p])})
    out = {"n": int(a.shape[1]), "pol": pols}
    if a.shape[0] == 2:
        cxy = covariance(series[0][0], series[1][0], circular=circular)
        cyx = covariance(series[1][0], series[0][0], circular=circular)
        common = 0.5 * (cxy + cyx)
        sigma = np.sqrt(pols[0]["phase_cov"][0] * pols[1]["phase_cov"][0])
        out["cross"] = {"phase_cov": common.tolist(), "common_phase_var": float(common[0]),
                        "phase_correlation": float(common[0] / sigma),
                        "amplitude_common_var": float(np.mean(series[0][1] * series[1][1])),
                        "window": window_from_covariance(common, float(sigma / np.sqrt(a.shape[1])))}
    return out

def summary(m: dict) -> dict:
    """Media das duas polarizacoes, so com os numeros que vao para a tabela."""
    avg = lambda f: float(np.mean([f(p) for p in m["pol"]]))
    s = {"inv_snr_total": avg(lambda p: p["inv_snr_total"]), "inv_snr_rest": avg(lambda p: p["inv_snr_rest"]),
         "image": avg(lambda p: p["image"]), "offset": avg(lambda p: p["offset"]), "ser": avg(lambda p: p["ser"]),
         "additive": avg(lambda p: p["radial_tangential"]["additive"]),
         "phase": avg(lambda p: p["radial_tangential"]["phase"]),
         "amplitude": avg(lambda p: p["radial_tangential"]["amplitude"]),
         "warp": avg(lambda p: p["radial_tangential"]["warp"]),
         "isi": avg(lambda p: p["memory"]["isi"]), "crosstalk": avg(lambda p: p["memory"]["crosstalk"]),
         "unexplained": avg(lambda p: p["memory"]["unexplained"]),
         "kurtosis": avg(lambda p: 0.5 * (p["kurtosis_re"] + p["kurtosis_im"])),
         "slips": int(sum(p["slips"]["slips"] for p in m["pol"])),
         "lost_blocks": int(sum(p["slips"]["lost_blocks"] for p in m["pol"])),
         "rotated_fraction": avg(lambda p: p["slips"]["rotated_fraction"]),
         "rotated_episodes": int(max(p["slips"]["rotated_episodes"] for p in m["pol"])),
         "min_correlation": float(min(p["slips"]["min_correlation"] for p in m["pol"])),
         "rotation_outer_minus_inner": avg(lambda p: p["radial_tangential"]["rotation_rad"][-1]
                                           - p["radial_tangential"]["rotation_rad"][0]),
         "conj_tap_1": avg(lambda p: 0.5 * (p["memory"]["own_taps"][1][1] + p["memory"]["own_taps"][-1][1])),
         "conj_symmetric": avg(lambda p: p["memory"]["conj_symmetric"]),
         "conj_antisymmetric": avg(lambda p: p["memory"]["conj_antisymmetric"]),
         "linear_symmetric": avg(lambda p: p["memory"]["linear_symmetric"]),
         "linear_antisymmetric": avg(lambda p: p["memory"]["linear_antisymmetric"]),
         "lin_tap_2": avg(lambda p: 0.5 * (p["memory"]["own_taps"][2][0] + p["memory"]["own_taps"][-2][0])),
         "snr_x_minus_y_db": float(10 * np.log10(m["pol"][1]["inv_snr_total"] / m["pol"][0]["inv_snr_total"]))
         if len(m["pol"]) == 2 else None,
         "block_std": avg(lambda p: p["along_block"]["std"]),
         "window_pol": [p["window"]["window"] for p in m["pol"]]}
    ser_g = [p["ser_gaussian"] for p in m["pol"]]
    if all(v for v in ser_g):
        s["ser_over_gaussian"] = float(np.mean([p["ser"] / p["ser_gaussian"] for p in m["pol"]]))
    if "cross" in m:
        s.update(common_phase_var=m["cross"]["common_phase_var"],
                 phase_correlation=m["cross"]["phase_correlation"], window_joint=m["cross"]["window"]["window"])
    return s

# ------------------------------------------------------------------ dados do VPI
def load_vpi_sets(root: Path, cache: Path | None):
    """Devolve (referencia {pol: simbolos decididos}, funcao carrega(L, P) -> (2, n) ou None)."""
    from compara_vpi import VPI_SETS, decide_16qam, load_vpi
    sim = root / "2026_SIMULATION" / "dados_4_canais"
    def cached(path: Path, key: str):
        if cache is None:
            return load_vpi(path)
        f = cache / f"{key}.npy"
        if not f.exists():
            cache.mkdir(parents=True, exist_ok=True)
            np.save(f, load_vpi(path))
        return np.load(f)
    ref = np.array([decide_16qam(cached(sim / "dados_ideais" / "dados_ref_canal2"
                                        / f"dados_4_canais_ref_{pol}_pol_canal2.csv", f"ref_{pol}"))
                    for pol in "XY"])
    def load(length: int, power: int):
        folder, pattern = VPI_SETS[length]
        files = [sim / folder / pattern.format(pol=pol, p=power) for pol in "XY"]
        if not all(f.exists() for f in files):
            return None
        return np.array([cached(f, f"{length}km_{pol}_{power}dBm") for f, pol in zip(files, "XY")])
    return ref, load

# ------------------------------------------------------------------ teste em dados sinteticos
def synthetic(seed: int = 0, n: int = 2 ** 17) -> None:
    """Dados de resposta conhecida, para conferir cada medida."""
    rng = np.random.default_rng(seed)
    lv = np.array([-3.0, -1.0, 1.0, 3.0])
    a = rng.choice(lv, (2, n)) + 1j * rng.choice(lv, (2, n))
    ps = 10.0
    cn = lambda scale, shape: scale * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)
    print("teste em dados sintéticos (16QAM, duas polarizações, n =", n, ")")

    # 1. so ruido aditivo, 1/SNR = 0,02
    x = a + cn(np.sqrt(0.02 * ps), a.shape)
    s = summary(measure(x, a, True))
    print(f"1. só ruído aditivo de 0,0200: aditivo {s['additive']:.4f} | fase {s['phase']:.5f} | amplitude "
          f"{s['amplitude']:.5f} | interferência {s['isi']:.5f} | diafonia {s['crosstalk']:.5f} | SER/gauss "
          f"{s['ser_over_gaussian']:.3f} | curtose {s['kurtosis']:+.3f}")

    # 2. ruido aditivo 0,01 + fase comum branca de 0,004 rad^2 + amplitude 0,002
    phi = np.sqrt(0.004) * rng.standard_normal(n)
    amp = np.sqrt(0.002) * rng.standard_normal((2, n))
    x = a * (1 + amp) * np.exp(1j * phi) + cn(np.sqrt(0.01 * ps), a.shape)
    s = summary(measure(x, a, True))
    print(f"2. aditivo 0,0100 + fase comum 0,0040 + amplitude 0,0020: aditivo {s['additive']:.4f} | fase "
          f"{s['phase']:.4f} | amplitude {s['amplitude']:.4f} | fase comum entre X e Y {s['common_phase_var']:.4f} "
          f"| correlação {s['phase_correlation']:.2f}")

    # 3. interferencia entre simbolos (vizinho +1 com 0,05 e -2 com 0,03j) e diafonia de 0,04
    x = a + 0.05 * np.roll(a, -1, axis=1) + 0.03j * np.roll(a, 2, axis=1) + 0.04 * a[::-1] + cn(np.sqrt(0.01 * ps), a.shape)
    m = measure(x, a, True)
    s = summary(m)
    taps = m["pol"][0]["memory"]["own_taps"]
    print(f"3. interferência 0,05^2 + 0,03^2 = 0,0034 e diafonia 0,04^2 = 0,0016, aditivo 0,0100: interferência "
          f"{s['isi']:.4f} | diafonia {s['crosstalk']:.4f} | sem explicação {s['unexplained']:.4f} | "
          f"coeficientes +1: {taps[1][0]:.3f}, -2: {taps[-2][0]:.3f}, +3: {taps[3][0]:.3f}")

    # 4. estimador de fase por media movel: fase de Wiener lenta, estimativa com janela W, conjunta
    for W in (21, 43):
        wiener = np.cumsum(np.sqrt(2 * np.pi * 1e6 / 59.84e9) * rng.standard_normal(n))
        y = a * np.exp(1j * wiener) + cn(np.sqrt(0.025 * ps), a.shape)
        z = np.sum(y * np.conj(a), axis=0)                         # estimador com os dados, as duas polarizacoes
        zf = np.convolve(z, np.ones(W), mode="same")
        x = y * np.exp(-1j * np.angle(zf))
        m = measure(x, a, False, edge=W)
        s = summary(m)
        c = np.array(m["cross"]["phase_cov"])
        h = (W - 1) // 2
        print(f"4. média móvel conjunta de W = {W}: janela lida {s['window_joint']} (entre X e Y), "
              f"{s['window_pol']} (cada polarização) | covariância entre X e Y nos atrasos 1, {h}, {h + 1}, "
              f"{W - 1}, {W + 5}: " + ", ".join(f"{c[l] * 1e4:+.2f}e-4" for l in (1, h, h + 1, W - 1, W + 5)))

    # 5. saltos: giro de 90 graus em dois trechos
    x = a[0] + cn(np.sqrt(0.02 * ps), n)
    x[30000:50000] *= 1j
    x[90000:] *= -1
    print("5. dois trechos girados (3 mudanças de giro):", {k: v for k, v in slips(x, a[0]).items()
                                                             if k in ("slips", "lost_blocks", "turn_histogram")})

# ------------------------------------------------------------------ linha de comando
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    sub.add_parser("sintetico")
    v = sub.add_parser("vpi")
    v.add_argument("folder", type=Path)
    v.add_argument("--lengths", type=int, nargs="+", default=[125, 130, 135, 140])
    v.add_argument("--powers", type=int, nargs="+", default=list(range(6, 19)))
    v.add_argument("--cache", type=Path)
    v.add_argument("-o", "--out", type=Path, required=True)
    s = sub.add_parser("simbolos")
    s.add_argument("files", type=Path, nargs="+")
    s.add_argument("-o", "--out", type=Path, required=True)
    t = sub.add_parser("tabela")
    t.add_argument("files", type=Path, nargs="+")
    a = ap.parse_args(argv)

    if a.command == "sintetico":
        return synthetic()
    if a.command == "tabela":
        for f in a.files:
            data = json.loads(f.read_text(encoding="utf-8"))
            print(f"== {f.name} ({data['source']})")
            print("L | P | 1/SNR resto | aditivo | fase | amplitude | deform. | interf. | diafonia | fase comum | "
                  "corr. XY | janela | SER/gauss | curtose | saltos | perdidos")
            for q in data["points"]:
                r = q["summary"]
                print(f"{q['length_km']} | {q['power_dbm']:4.1f} | {r['inv_snr_rest']:.5f} | {r['additive']:.5f} | "
                      f"{r['phase']:.5f} | {r['amplitude']:.5f} | {r['warp']:.5f} | {r['isi']:.5f} | "
                      f"{r['crosstalk']:.5f} | {r.get('common_phase_var', float('nan')):.5f} | "
                      f"{r.get('phase_correlation', float('nan')):.2f} | {r.get('window_joint')} | "
                      f"{r.get('ser_over_gaussian', float('nan')):.2f} | {r['kurtosis']:+.2f} | {r['slips']} | "
                      f"{r['lost_blocks']} | girado {r['rotated_fraction']:.3f}")
        return
    points = []
    if a.command == "vpi":
        ref, load = load_vpi_sets(a.folder, a.cache)
        for length in a.lengths:
            for power in a.powers:
                x = load(length, power)
                if x is None:
                    continue
                m = measure(x, ref, True)
                points.append({"length_km": length, "power_dbm": power, "summary": summary(m), "detail": m})
                r = points[-1]["summary"]
                print(f"VPI {length} km {power:2d} dBm: resto {r['inv_snr_rest']:.5f} | aditivo {r['additive']:.5f} | "
                      f"fase {r['phase']:.5f} | interf. {r['isi']:.5f} | diafonia {r['crosstalk']:.5f} | "
                      f"janela {r.get('window_joint')} | girado {r['rotated_fraction']:.3f} em "
                      f"{r['rotated_episodes']} trechos | correlação mínima {r['min_correlation']:.2f}", flush=True)
        result = {"source": "vpi", "points": points}
    else:
        for f in sorted(a.files):
            d = np.load(f)
            m = measure(d["received"], d["sent"], False, edge=8)
            points.append({"length_km": float(d["length_km"]), "power_dbm": float(d["power_dbm"]),
                           "mode": str(d["mode"]), "file": f.name, "summary": summary(m), "detail": m})
            r = points[-1]["summary"]
            print(f"simcopy {f.name}: resto {r['inv_snr_rest']:.5f} | aditivo {r['additive']:.5f} | fase "
                  f"{r['phase']:.5f} | interf. {r['isi']:.5f} | diafonia {r['crosstalk']:.5f} | janela "
                  f"{r.get('window_joint')} | girado {r['rotated_fraction']:.3f}", flush=True)
        result = {"source": "simcopy", "points": points}
    def compact(v):                      # 6 algarismos significativos: o arquivo fica 5 vezes menor
        if isinstance(v, float):
            return float(f"{v:.6g}")
        if isinstance(v, dict):
            return {k: compact(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [compact(x) for x in v]
        return v
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(compact(result), separators=(",", ":")), encoding="utf-8")
    print(f"gravado: {a.out}")

if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
