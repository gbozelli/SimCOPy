"""Experimento da linha A (docs/replicacao_vpi_diario.md): a busca de fase do VPI estima a fase
com mais ruido que a do simcopy?

O excesso de ruido do VPI em potencia baixa (4,7 % a 8 dBm em 125 km, nada de 15 dBm em diante)
tem a forma de uma diferenca so no termo A/P. O ruido de estimacao da fase vem do ruido aditivo e
escala com 1/P, e duas escolhas do bloco `BlindPhaseSearch` do VPI o aumentariam:
  `JointEstimation = No`   cada polarizacao estima a sua fase (o simcopy soma as duas, com o dobro
                           de amostras por janela);
  `Format = ReadFromFile`  as decisoes da busca sao feitas na grade de `256QAM_constellation.txt`,
                           que e mais densa e discrimina menos os angulos.
A escala do 16QAM sobre a grade de 256QAM nao esta no documento. Duas leituras:
  s = 1      os 16 pontos do 16QAM caem sobre os 16 pontos internos do 256QAM;
  s = 4,123  o 16QAM tem a potencia media do 256QAM (170, `ConstMeanPower` do equalizador).

Sobre os sinais guardados na saida do ADC, o receptor cego do simcopy (janela de 43) roda com a
busca de fase trocada. O resto da cadeia nao muda. Fica fora de src/.

Uso, na pasta repo:
  python scripts/experimentos/busca_fase_vpi.py PASTA_SINAIS SAIDA.json [--files 125km_8dBm ...]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from equalizador_gradiente import evaluate  # noqa: E402
from simcopy.kernels import dsp  # noqa: E402
from simcopy.kernels.demodulation import measure_blind  # noqa: E402
from simcopy.kernels.mapping import ROTATIONAL_SYMMETRY, constellation, map_bits, nearest_points  # noqa: E402

MODULATION = "16QAM"
WINDOW = 43
STATE = {"grid_scale": None, "joint": True, "fourth_power_span_deg": None, "coarse_window": None,
         "partition": False, "freq_block": None, "freq_offsets": None, "coarse_period_deg": None}
VARIANTS = {
    "base": dict(grid_scale=None, joint=True),
    "por_polarizacao": dict(grid_scale=None, joint=False),
    "grade256_s1": dict(grid_scale=1.0, joint=True),
    "grade256_s4.123": dict(grid_scale=np.sqrt(17.0), joint=True),
    "por_polarizacao_grade256_s1": dict(grid_scale=1.0, joint=False),
    "por_polarizacao_grade256_s4.123": dict(grid_scale=np.sqrt(17.0), joint=False),
    # dois estagios (`FourthPower = Yes`, `FineBlindEstimation = Yes`): fase grossa pela quarta
    # potencia e busca fina restrita a +-span em torno dela
    "quarta_fina_45": dict(grid_scale=None, joint=True, fourth_power_span_deg=45.0),
    "quarta_fina_22": dict(grid_scale=None, joint=True, fourth_power_span_deg=22.5),
    "quarta_fina_11": dict(grid_scale=None, joint=True, fourth_power_span_deg=11.25),
    "quarta_fina_6": dict(grid_scale=None, joint=True, fourth_power_span_deg=5.625),
    "por_polarizacao_quarta_fina_22": dict(grid_scale=None, joint=False, fourth_power_span_deg=22.5),
    "por_polarizacao_quarta_fina_11": dict(grid_scale=None, joint=False, fourth_power_span_deg=11.25),
    # a janela do estagio grosso separada da do fino (43)
    "quarta_fina_22_grossa21": dict(grid_scale=None, joint=True, fourth_power_span_deg=22.5, coarse_window=21),
    "quarta_fina_22_grossa11": dict(grid_scale=None, joint=True, fourth_power_span_deg=22.5, coarse_window=11),
    "por_polarizacao_quarta_fina_22_grossa85": dict(grid_scale=None, joint=False, fourth_power_span_deg=22.5,
                                                    coarse_window=85),
    "por_polarizacao_quarta_fina_22_grossa129": dict(grid_scale=None, joint=False, fourth_power_span_deg=22.5,
                                                     coarse_window=129),
    "por_polarizacao_quarta_fina_22_grossa257": dict(grid_scale=None, joint=False, fourth_power_span_deg=22.5,
                                                     coarse_window=257),
    "por_polarizacao_quarta_fina_22_grossa55": dict(grid_scale=None, joint=False, fourth_power_span_deg=22.5,
                                                    coarse_window=55),
    "por_polarizacao_quarta_fina_22_grossa65": dict(grid_scale=None, joint=False, fourth_power_span_deg=22.5,
                                                    coarse_window=65),
    # estagio grosso so com os aneis de dentro e de fora do 16QAM (particao em QPSK), fase de y^4/|y|^4
    "particionada_22": dict(grid_scale=None, joint=True, fourth_power_span_deg=22.5, partition=True),
    "por_polarizacao_particionada_22": dict(grid_scale=None, joint=False, fourth_power_span_deg=22.5,
                                            partition=True),
    "por_polarizacao_particionada_22_grossa21": dict(grid_scale=None, joint=False, fourth_power_span_deg=22.5,
                                                     partition=True, coarse_window=21),
    "particionada_22_grossa21": dict(grid_scale=None, joint=True, fourth_power_span_deg=22.5, partition=True,
                                     coarse_window=21),
}
for _b in (1024, 512, 256):           # hipotese 3: frequencia por blocos
    VARIANTS[f"freq_bloco_{_b}"] = dict(freq_block=_b)
    VARIANTS[f"dois_estagios_32_freq_bloco_{_b}"] = dict(joint=False, fourth_power_span_deg=22.5, partition=True,
                                                         coarse_window=32, freq_block=_b)
    VARIANTS[f"dois_estagios_43_freq_bloco_{_b}"] = dict(joint=False, fourth_power_span_deg=22.5, partition=True,
                                                         coarse_window=43, freq_block=_b)
for _w in (43, 32):                     # hipotese 4 (exploratoria): grossa desdobrada com periodo de 45 graus
    for _span in (22.5, 5.625):
        VARIANTS[f"periodo45_grossa{_w}_fina{_span:g}"] = dict(joint=False, fourth_power_span_deg=_span, partition=True,
                                                               coarse_window=_w, coarse_period_deg=45.0)
for _w in (43, 32):                     # so a quarta potencia, sem busca fina, desdobrada com periodo de 45 graus
    VARIANTS[f"quarta45_janela{_w}_aneis"] = dict(joint=False, fourth_power_span_deg=0.0, partition=True,
                                                 coarse_window=_w, coarse_period_deg=45.0)
    VARIANTS[f"quarta45_janela{_w}_todos"] = dict(joint=False, fourth_power_span_deg=0.0, partition=False,
                                                 coarse_window=_w, coarse_period_deg=45.0)
for _w in (43, 32, 27, 25, 21):         # o mesmo, com a fase estimada em conjunto nas duas polarizacoes
    VARIANTS[f"conjunta45_grossa{_w}_fina22.5"] = dict(joint=True, fourth_power_span_deg=22.5, partition=True,
                                                       coarse_window=_w, coarse_period_deg=45.0)
    VARIANTS[f"conjunta45_quarta{_w}_aneis"] = dict(joint=True, fourth_power_span_deg=0.0, partition=True,
                                                    coarse_window=_w, coarse_period_deg=45.0)
for _w in (25, 29, 32, 37):          # janela do estagio grosso entre 21 e 43 (o VPI tem FilterLength = 32)
    VARIANTS[f"por_polarizacao_particionada_22_grossa{_w}"] = dict(
        grid_scale=None, joint=False, fourth_power_span_deg=22.5, partition=True, coarse_window=_w)


def decisions(r: np.ndarray) -> np.ndarray:
    """Decisao da busca de fase: no 16QAM, ou na grade de 256QAM (niveis impares de -15 a 15) com
    o 16QAM (niveis +-1 e +-3) multiplicado pela escala. Devolve na escala de r."""
    s = STATE["grid_scale"]
    if s is None:
        return nearest_points(r, MODULATION)
    k = 3.0 * s                         # o 16QAM do simcopy tem niveis +-1/3 e +-1
    q = lambda v: np.clip(2.0 * np.round((v - 1.0) / 2.0) + 1.0, -15.0, 15.0)
    v = r * k
    return (q(v.real) + 1j * q(v.imag)) / k


def fourth_power_search(y: np.ndarray, test_angles: int, window: int, span: float) -> np.ndarray:
    """Busca em dois estagios. Grossa: arg(-soma de y^4 na janela)/4, desdobrada com periodo de
    90 graus (no 16QAM a media de a^4 e negativa, por isso o sinal de menos). Fina: os test_angles
    angulos entre -span e +span em torno da grossa, com a mesma soma de distancias da busca cega.
    Se o ruido virar o sinal da soma na janela, a grossa erra 45 graus, e a fina so corrige se
    span passar de 45 graus."""
    kernel = np.ones(min(int(window), y.shape[1]))
    coarse_kernel = np.ones(min(int(STATE.get("coarse_window") or window), y.shape[1]))
    y4 = y ** 4
    if STATE.get("partition"):
        # aneis do 16QAM na escala do simcopy: 2/9, 10/9 e 2. Fica com o de dentro e o de fora.
        p = np.abs(y) ** 2
        y4 = np.where((p < 0.667) | (p > 1.556), y4 / np.maximum(p, 1e-12) ** 2, 0.0)
    s4 = np.convolve(np.sum(y4, axis=0), coarse_kernel, mode="same")
    # hipotese 4 (exploratoria): desdobrar com periodo de 45 graus deixa a grossa sem distinguir o estado
    # alinhado do girado de 45 graus, como parece acontecer no VPI (fracao girada perto de 0,5)
    period = np.radians(STATE.get("coarse_period_deg") or 90.0)
    coarse = np.unwrap(np.angle(-s4) / 4.0, period=period)
    if period < np.pi / 2:
        # com periodo de 45 graus a grossa nao sabe em que estado esta. Em cada simbolo ela difere da
        # estimativa de quarta potencia com o sinal certo (periodo de 90 graus, mesma janela) de 0 ou
        # de 45 graus. O ramo e o da maioria num trecho limpo de 2048 simbolos, longe da borda do bloco,
        # onde o sinal ainda nao presta (a guarda do simcopy e de ~500 simbolos). A comparacao e local:
        # o ruido de fase dos lasers anda dezenas de graus em 2048 simbolos.
        ref = slice(1024, 3072) if y.shape[1] >= 4096 else slice(0, y.shape[1])
        if np.mean(np.real(np.exp(4j * (coarse[ref] - np.angle(-s4[ref]) / 4.0)))) < 0:
            coarse = coarse + period
    if span == 0.0:
        return coarse
    offsets = np.linspace(-span, span, test_angles)
    cost = np.empty((test_angles, y.shape[1]))
    for b, off in enumerate(offsets):
        r = y * np.exp(-1j * (coarse + off))
        cost[b] = np.convolve(np.sum(np.abs(r - decisions(r)) ** 2, axis=0), kernel, mode="same")
    return coarse + offsets[np.argmin(cost, axis=0)]


def blind_phase_search(symbols, modulation, test_angles=64, window=21):
    """`dsp.blind_phase_search` com as decisoes de `decisions`, ou a busca em dois estagios."""
    y = np.atleast_2d(np.asarray(symbols, dtype=complex))
    if STATE.get("fourth_power_span_deg") is not None:
        return fourth_power_search(y, test_angles, window, np.radians(STATE["fourth_power_span_deg"]))
    sector = 2.0 * np.pi / ROTATIONAL_SYMMETRY[modulation]
    angles = np.arange(test_angles) * sector / test_angles
    kernel = np.ones(min(int(window), y.shape[1]))
    cost = np.empty((test_angles, y.shape[1]))
    for b, angle in enumerate(angles):
        r = y * np.exp(-1j * angle)
        cost[b] = np.convolve(np.sum(np.abs(r - decisions(r)) ** 2, axis=0), kernel, mode="same")
    return np.unwrap(angles[np.argmin(cost, axis=0)], period=sector)


def carrier_recovery(symbols, modulation, test_angles=64, window=21):
    """`dsp.carrier_recovery`, com a busca conjunta das duas polarizacoes opcional."""
    y = np.atleast_2d(np.asarray(symbols, dtype=complex))
    rms = np.sqrt(np.mean(np.abs(constellation(modulation)) ** 2))
    block = STATE.get("freq_block")
    if block:
        # hipotese 3: o desvio de frequencia estimado por blocos de `block` simbolos (o bloco de
        # frequencia do VPI trabalha por blocos), com a fase continua de um bloco para o outro
        n = y.shape[1]
        ramp, offsets, start = np.empty(n), [], 0.0
        for s0 in range(0, n, block):
            seg = y[:, s0:s0 + block]
            f_b = dsp.estimate_frequency_offset(seg, 1.0, 1.0 / 8.0)
            offsets.append(f_b)
            ramp[s0:s0 + seg.shape[1]] = start + 2.0 * np.pi * f_b * np.arange(seg.shape[1])
            start = ramp[s0 + seg.shape[1] - 1] + 2.0 * np.pi * f_b
        offset = float(np.mean(offsets))
        STATE["freq_offsets"] = offsets
    else:
        offset = dsp.estimate_frequency_offset(y, 1.0, 1.0 / 8.0)
        ramp = 2.0 * np.pi * offset * np.arange(y.shape[1])
    rows = [dsp.decision_scale(row * np.exp(-1j * ramp) / np.sqrt(np.mean(np.abs(row) ** 2)) * rms, modulation)
            for row in y]
    if len(rows) == 2 and STATE["joint"] and STATE.get("coarse_period_deg"):
        # com periodo de 45 graus as buscas por polarizacao podem estar em estados diferentes, e a
        # diferenca entre elas sairia errada. A fase relativa entre as polarizacoes e constante no
        # bloco (os lasers sao os mesmos): sai da quarta potencia cruzada, sem depender de busca.
        delta = np.angle(np.sum(rows[1] ** 4 * np.conj(rows[0] ** 4))) / 4.0
        joint = blind_phase_search(np.array([rows[0], rows[1] * np.exp(-1j * delta)]), modulation,
                                   test_angles, window)
        out = [dsp.decision_scale(row * np.exp(-1j * phi), modulation) for row, phi in zip(rows, (joint, joint + delta))]
        return np.array(out), np.array([ramp + joint, ramp + joint + delta]), float(offset)
    phis = [blind_phase_search(row, modulation, test_angles, window) for row in rows]
    if len(rows) == 2 and STATE["joint"]:
        symmetry = ROTATIONAL_SYMMETRY[modulation]
        delta = np.angle(np.mean(np.exp(1j * symmetry * (phis[1] - phis[0])))) / symmetry
        joint = blind_phase_search(np.array([rows[0], rows[1] * np.exp(-1j * delta)]), modulation,
                                   test_angles, window)
        phis = [joint, joint + delta]
    out = [dsp.decision_scale(row * np.exp(-1j * phi), modulation) for row, phi in zip(rows, phis)]
    return np.array(out), np.array([ramp + phi for phi in phis]), float(offset)


dsp.carrier_recovery = carrier_recovery     # blind_receiver procura o nome no modulo a cada chamada


def rotation(z: np.ndarray, a: np.ndarray, guard: int, block: int = 128) -> dict:
    """Giro de cada bloco de `block` simbolos em relacao aos dados, modulo 90 graus, depois do
    alinhamento: fracao girada mais de 22,5 graus e o giro medio desses blocos [graus]."""
    from simcopy.kernels.demodulation import align_to_transmitted
    y, _, _ = align_to_transmitted(z, a, MODULATION)
    y, at = y[:, guard:a.shape[1] - guard], a[:, guard:a.shape[1] - guard]
    n = (at.shape[1] // block) * block
    c = np.sum((y[:, :n] * np.conj(at[:, :n])).reshape(at.shape[0], -1, block), axis=2).sum(axis=0)
    ang = np.degrees(np.angle(c ** 4)) / 4.0
    rot = np.abs(ang) > 22.5
    return {"rotated_fraction": float(np.mean(rot)),
            "rotated_mean_abs_deg": float(np.mean(np.abs(ang[rot]))) if rot.any() else None}


def entries(z: np.ndarray, a: np.ndarray, guard: int, window: int = 64) -> dict:
    """Trocas de estado de fase (passos de 45 graus) por polarizacao, com o mesmo detector usado nos
    dados do VPI (`fronteiras_vpi.phase_states`), a fracao do tempo em 45 graus e em saltos de 90 ou
    180 graus. Tambem a SNR com a fase corrigida pelos simbolos conhecidos numa janela deslizante de
    `window`: nao depende de a busca cega ter girado ou saltado, entao serve de eixo comum para o VPI
    e para todas as variantes."""
    from simcopy.kernels.demodulation import align_to_transmitted
    from fronteiras_vpi import phase_states
    y, _, _ = align_to_transmitted(z, a, MODULATION)
    y, at = y[:, guard:a.shape[1] - guard], a[:, guard:a.shape[1] - guard]
    per_pol, noise = [], []
    for q in range(at.shape[0]):
        st = phase_states(y[q], at[q], window)
        per_pol.append({"changes": len(st["changes"]), "frac_45": st["frac_45"], "frac_slipped": st["frac_slipped"]})
        c = np.convolve(y[q] * np.conj(at[q]), np.ones(window), mode="same")
        yy = y[q] * np.exp(-1j * np.angle(c))
        g = np.sum(yy * np.conj(at[q])) / np.sum(np.abs(at[q]) ** 2)
        noise.append(np.mean(np.abs(yy / g - at[q]) ** 2) / np.mean(np.abs(at[q]) ** 2))
    mean = lambda k: float(np.mean([q[k] for q in per_pol]))
    return {"changes": mean("changes"), "frac_45": mean("frac_45"), "frac_slipped": mean("frac_slipped"),
            "states_per_pol": per_pol, "n_symbols": int(at.shape[1]),
            "snr_known_phase_db": float(-10 * np.log10(np.mean(noise)))}


def base_equalized(x: np.ndarray) -> np.ndarray:
    """Saida do equalizador do receptor cego de hoje (`dsp.blind_receiver` com a busca de hoje), antes
    da ultima recuperacao de portadora. Com `--fixed-equalizer`, cada variante so troca a recuperacao
    de portadora aplicada sobre esta saida: o equalizador nao depende da busca de fase, como no VPI
    (no simcopy o ajuste por decisao usa a saida da busca, e uma busca que falha estragaria o
    equalizador e realimentaria a falha)."""
    STATE.update({"grid_scale": None, "joint": True, "fourth_power_span_deg": None, "coarse_window": None,
                  "partition": False, "freq_block": None, "freq_offsets": None, "coarse_period_deg": None})
    xs = np.atleast_2d(np.asarray(x, dtype=complex))
    if xs.shape[0] == 2:
        xs = dsp.stokes_polarization_demux(xs)
    y = dsp.adaptive_equalizer(xs, constellation(MODULATION), 7, 1e-3, 4, 0)
    z, phases, _ = dsp.carrier_recovery(y, MODULATION, 64, WINDOW)
    u = dsp._regressors(xs / np.sqrt(np.mean(np.abs(xs) ** 2, axis=1, keepdims=True)), 7)
    gram = u.conj().T @ u
    for _ in range(4):
        target = nearest_points(z, MODULATION) * np.exp(1j * phases)
        y = (u @ np.linalg.solve(gram, u.conj().T @ target.T)).T
        z, phases, _ = dsp.carrier_recovery(y, MODULATION, 64, WINDOW)
    return y


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("folder", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--files", nargs="*", default=None)
    ap.add_argument("--variants", nargs="*", default=None)
    ap.add_argument("--fixed-equalizer", action="store_true",
                    help="equaliza uma vez com o receptor de hoje e aplica so a recuperacao de portadora de cada variante")
    ap.add_argument("--scales", type=float, nargs="*", default=[],
                    help="escalas a mais da grade de 256QAM (variantes conjunta e por polarizacao)")
    args = ap.parse_args()
    for s in args.scales:
        VARIANTS[f"grade256_s{s:g}"] = dict(grid_scale=s, joint=True)
        VARIANTS[f"por_polarizacao_grade256_s{s:g}"] = dict(grid_scale=s, joint=False)
    if args.variants is None:
        args.variants = list(VARIANTS)
    sys.stdout.reconfigure(encoding="utf-8")
    results = json.loads(args.out.read_text(encoding="utf-8")) if args.out.exists() else []
    done = {(r["file"], r["variant"]) for r in results}
    files = sorted(args.folder.glob("*.npz"))
    if args.files:
        files = [f for f in files if any(f.name.startswith(p + "_") for p in args.files)]
    for f in files:
        d = np.load(f)
        x = d["adc"].astype(complex)
        a = np.array([map_bits(row, str(d["modulation"])) for row in d["bits"]])
        guard = int(d["guard"])
        y_eq = None
        if args.fixed_equalizer and any((f.name, n) not in done for n in args.variants):
            y_eq = base_equalized(x)
        for name in args.variants:
            if (f.name, name) in done:
                continue
            STATE.update({"grid_scale": None, "joint": True, "fourth_power_span_deg": None, "coarse_window": None,
                          "partition": False, "freq_block": None, "freq_offsets": None, "coarse_period_deg": None,
                          **VARIANTS[name]})
            t0 = time.perf_counter()
            if y_eq is not None:
                z = dsp.carrier_recovery(y_eq, MODULATION, 64, WINDOW)[0]
            else:
                z, _ = dsp.blind_receiver(x, MODULATION, window=WINDOW)
            m, _ = measure_blind(z, d["bits"], MODULATION, guard)
            r = {"file": f.name, "variant": name, "length_km": float(d["length_km"]),
                 "power_dbm": float(d["power_dbm"]), "seed": int(d["seed"]), **evaluate(z, a, guard),
                 "equalizer": "fixo, o de hoje" if y_eq is not None else "o de hoje, com a busca da variante",
                 "ber": float(m.ber_counted), "ser": float(m.ser), "phase_slips": int(m.phase_slips),
                 "freq_std_hz": (float(np.std(STATE["freq_offsets"]) * float(d["symbol_rate"]))
                                 if STATE.get("freq_offsets") and "symbol_rate" in d.files else
                                 (float(np.std(STATE["freq_offsets"]) * 59.84e9) if STATE.get("freq_offsets") else None)),
                 "rotation": rotation(z, a, guard), **entries(z, a, guard), "seconds": time.perf_counter() - t0}
            results.append(r)
            args.out.write_text(json.dumps(results), encoding="utf-8")
            g = r["rotation"]
            mean = f"{g['rotated_mean_abs_deg']:.1f} graus" if g["rotated_mean_abs_deg"] is not None else "-"
            print(f"{f.stem:22s} {name:32s}: BER {r['ber']:.2e} | SER {r['ser']:.2e} | SNR {r['snr_db']:6.2f} | "
                  f"girado {g['rotated_fraction']:.3f} (em media {mean}) | trocas {r['changes']:.1f}, "
                  f"{100 * r['frac_45']:.0f} % em 45, {100 * r['frac_slipped']:.0f} % saltado | "
                  f"SNR com fase conhecida {r['snr_known_phase_db']:.2f} | {r['seconds']:.0f} s", flush=True)


if __name__ == "__main__":
    main()
