"""Compara o simcopy com as constelacoes do VPI no cenario de 4 canais DP-16QAM
(examples/400zr_4ch_75ghz.yaml; docs/comparacao_vpi.md).

Roda o canal sob teste (index 2, o Ch_X2/Ch_Y2 do VPI) para cada potencia e comprimento, e mede nos
simbolos equalizados as mesmas grandezas que mede nas constelacoes exportadas do VPI:

    x = g*a + h*conj(a) + c + ruido

- |h/g|: imagem de quadratura (desbalanco de fase do modulador IQ);
- |c/g|/rms(a): deslocamento constante (vazamento de portadora pela extincao finita dos MZM);
- SNR do resto, sem a imagem e sem o deslocamento;
- SNR total (contra g*a) e taxa de erro de simbolo.

`PowerdBm` do VPI (decisao D-22 do pacote, em aberto): --power-is total (a potencia do booster e a
soma dos 4 canais) ou per_channel. No simcopy `booster.output_power` e por canal.

    python scripts/compara_vpi.py simular --lengths 125 --powers 8 10 12 13 14 15 16 17 \
        --n-symbols 16384 -o resultados/vpi/simcopy_total.json
    python scripts/compara_vpi.py vpi "<pasta>/dados4canais/extraido" -o resultados/vpi/vpi.json
    python scripts/compara_vpi.py figura resultados/vpi/vpi.json resultados/vpi/simcopy_total.json \
        -o resultados/vpi/comparacao.png

Os dados do VPI nao ficam no repositorio: a pasta e passada na linha de comando.
"""
from __future__ import annotations
import argparse
import dataclasses as dc
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

LEVELS = np.array([-3.0, -1.0, 1.0, 3.0])

def decide_16qam(x: np.ndarray) -> np.ndarray:
    """Ponto mais proximo da grade +-1, +-3."""
    q = lambda v: LEVELS[np.clip(np.round((v + 3.0) / 2.0), 0, 3).astype(int)]
    return q(x.real) + 1j * q(x.imag)

def decompose(x: np.ndarray, a: np.ndarray) -> dict:
    """Ajuste x = g*a + h*conj(a) + c de minimos quadrados, e as grandezas da comparacao.
    `a` sao os simbolos transmitidos, em qualquer escala."""
    m = np.column_stack([a, np.conj(a), np.ones_like(a)])
    (g, h, c), *_ = np.linalg.lstsq(m, x, rcond=None)
    signal = np.sum(np.abs(g * a) ** 2)
    g1 = np.sum(x * np.conj(a)) / np.sum(np.abs(a) ** 2)          # ganho so, como a medida do simcopy
    y = x / g1
    scale = np.sqrt(np.mean(np.abs(a) ** 2) / 10.0)               # leva `a` para a grade +-1, +-3
    return {
        "image": float(abs(h / g)),
        "offset": float(abs(c / g) / np.sqrt(np.mean(np.abs(a) ** 2))),
        "snr_residual_db": float(10 * np.log10(signal / np.sum(np.abs(x - m @ np.array([g, h, c])) ** 2))),
        "snr_db": float(10 * np.log10(np.sum(np.abs(a) ** 2) / np.sum(np.abs(y - a) ** 2))),
        "ser": float(np.mean(decide_16qam(y / scale) != np.round(a / scale))),
    }

def mean_of(items: list[dict]) -> dict:
    """Media das polarizacoes: SNR pela media do ruido (1/SNR), o resto pela media simples."""
    out = {}
    for k in items[0]:
        v = np.array([i[k] for i in items])
        out[k] = float(-10 * np.log10(np.mean(10 ** (-v / 10)))) if k.endswith("_db") else float(v.mean())
    return out

# ------------------------------------------------------------------------------- VPI ----------
def load_vpi(path: Path) -> np.ndarray:
    return np.loadtxt(path, delimiter=",", comments=("//", "#")) @ np.array([1, 1j])

VPI_SETS = {     # comprimento [km]: (subpasta, padrao do nome do arquivo)
    125: ("dados_distorcidos/dados_125km/dados_125km_canal_2", "dados_4_canais_canal_2_125km_{pol}_pol_{p}dBm.csv"),
    130: ("dados_distorcidos/dados_130km/dados_130km_canal_2", "dados_4_canais_canal_2_130km_{pol}_pol_{p}dBm.csv"),
    135: ("dados_distorcidos/dados_135km/dados_135km_canal_2", "dados_4_canais_canal_2_135km_{pol}_pol_{p}dBm.csv"),
    140: ("dados_distorcidos/dados_140m/dados_140km_canal_2", "dados_4_canais_canal_2_140km_{pol}_pol_{p}dBm.csv"),
}

def analyse_vpi(folder: Path, lengths, powers) -> dict:
    """Decomposicao das constelacoes do VPI. Os simbolos transmitidos saem da constelacao de
    referencia (dados_ideais), que tem SNR de 22 dB contra a grade: a decisao nela nao erra."""
    base = folder / "2026_SIMULATION" / "dados_4_canais"
    sent = {pol: decide_16qam(load_vpi(base / "dados_ideais" / "dados_ref_canal2"
                                       / f"dados_4_canais_ref_{pol}_pol_canal2.csv")) for pol in "XY"}
    points = []
    for length in lengths:
        sub, pattern = VPI_SETS[length]
        for p in powers:
            files = {pol: base / sub / pattern.format(pol=pol, p=p) for pol in "XY"}
            if not all(f.exists() for f in files.values()):
                continue
            per_pol = [decompose(load_vpi(files[pol]), sent[pol]) for pol in "XY"]
            points.append({"length_km": length, "power_dbm": p, **mean_of(per_pol),
                           "per_pol": per_pol})
            print(f"VPI {length} km {p:2d} dBm: SNR = {points[-1]['snr_db']:6.2f} dB | resto = "
                  f"{points[-1]['snr_residual_db']:6.2f} dB | imagem = {points[-1]['image']:.4f} | "
                  f"deslocamento = {points[-1]['offset']:.4f} | SER = {points[-1]['ser']:.3e}", flush=True)
    return {"source": "vpi", "points": points}

# --------------------------------------------------------------------------- simcopy ----------
def _with_dsp(system, options: dict | None):
    """O sistema com os campos `options` trocados em rx.dsp de todos os canais."""
    if not options:
        return system
    return dc.replace(system, channels=tuple(
        dc.replace(c, rx=dc.replace(c.rx, dsp=dc.replace(c.rx.dsp, **options)))
        for c in system.channels))

def simulate(config: Path, lengths, powers, n_symbols: int, sps: int, power_is: str, channel: int,
             pmd: bool, seed: int, dsp_options: dict | None = None,
             checkpoint: Path | None = None, second: list[tuple[dict, Path]] | None = None,
             save_symbols: Path | None = None, save_signals: Path | None = None) -> dict:
    """Varre potencias e comprimentos. Com `checkpoint`, grava o resultado parcial nesse
    arquivo depois de cada ponto: uma varredura interrompida nao perde o que ja rodou.

    `second` = lista de (campos de rx.dsp, arquivo): outras cadeias de DSP medidas em cada ponto
    sobre o mesmo sinal propagado (os estagios ate o receptor sao reaproveitados), cada uma
    gravada no seu arquivo. Serve para medir varios DSP sem propagar de novo.

    `save_symbols`: pasta em que os simbolos equalizados e os transmitidos de cada ponto sao
    gravados (.npz), para as medidas de scripts/medidas_constelacao.py.

    `save_signals`: pasta em que o sinal na saida do ADC do canal sob teste (2 amostras por
    simbolo, antes do DSP cego) e os bits transmitidos sao gravados (.npz), uma vez por ponto.
    So existe nas cadeias com mode = blind. Serve para testar outro DSP sem propagar de novo."""
    from simcopy.config.io import load
    from simcopy.config.simulation import SimulationConfig
    from simcopy.core.units import dbm_to_w
    from simcopy.kernels.mapping import map_bits
    from simcopy.pipelines.runner import run

    base, simulation = load(config)
    system = _with_dsp(base, dsp_options)
    n = system.n_channels
    rate = system.channels[0].symbol_rate
    simulation = SimulationConfig.from_symbols(
        rate, sps, n_symbols, root_seed=seed, ssfm_max_step=simulation.ssfm_max_step,
        ssfm_max_phase_change_deg=simulation.ssfm_max_phase_change_deg)
    second = list(second or [])
    points, points2 = [], [[] for _ in second]
    header = {"source": "simcopy", "config": str(config), "n_symbols": n_symbols,
              "samples_per_symbol": sps, "power_is": power_is, "channel": channel, "pmd": pmd,
              "root_seed": seed, "dsp_options": dsp_options or {}}
    header2 = [{**header, "dsp_options": {**(dsp_options or {}), **options}} for options, _ in second]

    def measure(r, length, p, t0) -> dict:
        m = r.metrics[channel]
        c = next(x for x in r.system.channels if x.index == channel)
        bits = r.output(f"ch{channel}.tx.bits").data
        k = c.modulation.bits_per_symbol
        sent = np.array([map_bits(row[m.first_symbol * k:m.last_symbol * k], c.modulation.name)
                         for row in bits])
        per_pol = [decompose(r.symbols[channel][q], sent[q]) for q in range(sent.shape[0])]
        q = {"length_km": length, "power_dbm": p, **mean_of(per_pol), "per_pol": per_pol,
             "ber": m.ber_counted, "bit_errors": m.errors, "mi": m.mi, "gmi": m.gmi,
             "phase_slips": m.phase_slips, "n_bits": m.n_bits,
             "received_power_dbm": float(10 * np.log10(r.received_power(channel)) + 30),
             "seconds": time.perf_counter() - t0}
        mode = c.rx.dsp.mode
        if save_signals is not None and mode == "blind":
            name = save_signals / f"{length}km_{p:g}dBm_semente{seed}.npz"
            if not name.exists():
                save_signals.mkdir(parents=True, exist_ok=True)
                rc = next(x for x in r.resolved.channels if x.index == channel)
                np.savez_compressed(name, adc=r.output(f"ch{channel}.rx.adc").samples.astype(np.complex64),
                                    bits=np.asarray(bits, dtype=np.uint8), guard=rc.edge_guard_symbols,
                                    modulation=c.modulation.name, symbol_rate=c.symbol_rate,
                                    length_km=length, power_dbm=p, seed=seed)
        if save_symbols is not None:
            save_symbols.mkdir(parents=True, exist_ok=True)
            name = f"{length}km_{p:g}dBm_{mode}_janela{c.rx.dsp.phase_window}_semente{seed}.npz"
            np.savez_compressed(save_symbols / name, received=r.symbols[channel], sent=sent,
                                length_km=length, power_dbm=p, mode=mode, seed=seed,
                                phase_window=c.rx.dsp.phase_window)
        print(f"simcopy {length} km {p:5.1f} dBm ({power_is}, {mode}): SNR = {q['snr_db']:6.2f} dB | "
              f"resto = {q['snr_residual_db']:6.2f} dB | imagem = {q['image']:.4f} "
              f"| deslocamento = {q['offset']:.4f} | SER = {q['ser']:.3e} | "
              f"BER = {m.ber_counted:.3e} | GMI = {m.gmi:.3f} | saltos = {m.phase_slips} | "
              f"{q['seconds']:.0f} s", flush=True)
        return q

    def save(path, head, pts):
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps({**head, "points": pts}, indent=1), encoding="utf-8")

    for length in lengths:
        fiber = dc.replace(system.link.spans[0].fiber, length=length * 1e3)
        if not pmd:
            fiber = dc.replace(fiber, birefringence=dc.replace(fiber.birefringence, enabled=False))
        previous = None
        for p in powers:
            per_channel = dbm_to_w(p) / (n if power_is == "total" else 1)
            link = dc.replace(system.link, spans=(dc.replace(system.link.spans[0], fiber=fiber),),
                              booster=dc.replace(system.link.booster, output_power=per_channel))
            t0 = time.perf_counter()
            r = run(dc.replace(system, link=link), simulation, previous=previous)
            previous = r
            points.append(measure(r, length, p, t0))
            save(checkpoint, header, points)
            for i, (options, path) in enumerate(second):
                t0 = time.perf_counter()
                r = run(dc.replace(_with_dsp(system, options), link=link), simulation,
                        previous=previous)
                previous = r
                points2[i].append(measure(r, length, p, t0))
                save(path, header2[i], points2[i])
    return {**header, "points": points}

# ---------------------------------------------------------------------- ajuste e figura -------
def fit_abc(points: list[dict], length: int, key: str = "snr_db", window=(8, 16)) -> dict | None:
    """Ajuste 1/SNR = A/P + B + C*P^2 (P em mW, a potencia da varredura), como
    reference/vpi/ajuste_abc.py do pacote, numa janela de potencias em que o DSP do VPI converge."""
    sel = [q for q in points if q["length_km"] == length and window[0] <= q["power_dbm"] <= window[1]]
    if len(sel) < 4:
        return None
    p = 10 ** (np.array([q["power_dbm"] for q in sel]) / 10)
    inv = 10 ** (-np.array([q[key] for q in sel]) / 10)
    m = np.column_stack([1 / p, np.ones_like(p), p ** 2])
    (a, b, c), *_ = np.linalg.lstsq(m, inv, rcond=None)
    residual = 10 * np.log10((m @ np.array([a, b, c])) / inv)
    return {"A": float(a), "B": float(b), "C": float(c), "floor_db": float(-10 * np.log10(b)) if b > 0 else None,
            "rms_db": float(np.sqrt(np.mean(residual ** 2))), "n": len(sel)}

def _label(s: dict, f: Path) -> str:
    """Legenda de um conjunto: o VPI, ou o modo do DSP do simcopy (do cabecalho do JSON)."""
    if s["source"] == "vpi":
        return "VPI"
    options = s.get("dsp_options") or {}
    if options.get("mode") == "blind":
        return f"simcopy, DSP cego (janela de fase de {options.get('phase_window', 21)})"
    return "simcopy, símbolos conhecidos" + ("" if f.stem == "simcopy" else f" ({f.stem})")

POWER_LABEL = "potência total de lançamento, soma dos 4 canais (dBm)"

def _matplotlib():
    """matplotlib sem janela e o estilo das figuras do simcopy (src/simcopy/plotting.py)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from simcopy import plotting
    return matplotlib, plt, plotting

def figure(files: list[Path], out: Path, lengths=None):
    """Grade de comparacao: uma linha por comprimento, com SNR total, SNR sem imagem e sem
    deslocamento, SER e deslocamento, contra a potencia."""
    mpl, plt, pl = _matplotlib()
    sets = [json.loads(Path(f).read_text(encoding="utf-8")) for f in files]
    all_lengths = sorted({q["length_km"] for s in sets for q in s["points"]})
    lengths = lengths or all_lengths
    panels = (("snr_db", "SNR total (dB)"), ("snr_residual_db", "SNR sem imagem e\nsem deslocamento (dB)"),
              ("ser", "SER"), ("offset", "deslocamento |c|/rms"))
    sim_colors = iter(["blue", "red", "green", "orange", "purple", "gray"])
    styles = []
    for s in sets:
        styles.append(dict(marker="o", ms=5, mfc="none", color=pl.COLORS["VPI"]) if s["source"] == "vpi"
                      else dict(marker="s", ms=3, color=next(sim_colors)))
    with mpl.rc_context(pl.FIGURE_STYLE):
        fig, axes = plt.subplots(len(lengths), 4, figsize=(13, 2.9 * len(lengths)), squeeze=False)
        letters = iter("abcdefghijklmnopqrstuvwxyz" * 3)
        for row, length in zip(axes, lengths):
            for s, f, style in zip(sets, files, styles):
                pts = sorted((q for q in s["points"] if q["length_km"] == length), key=lambda q: q["power_dbm"])
                if not pts:
                    continue
                x = [q["power_dbm"] for q in pts]
                for a, (key, ylabel) in zip(row, panels):
                    y = [q[key] for q in pts]
                    (a.semilogy if key == "ser" else a.plot)(x, y, label=_label(s, Path(f)), **style)
                    a.set(xlabel="potência (dBm)", ylabel=ylabel)
            row[0].set_ylim(8, 20)
            row[1].set_ylim(8, 22)
            row[2].set_ylim(1e-4, 1)
            lo, hi = row[3].get_ylim()
            row[3].set_ylim(lo - 0.3 * (hi - lo), hi)       # faixa livre para o rotulo
            # canto livre de cada coluna: SNR em cima (curvas abaixo de 19 dB), SER e deslocamento
            # embaixo (em cima fica o VPI quando perde a fase)
            for a, corner in zip(row, ("upper left", "upper left", "lower left", "lower left")):
                pl.panel_label(a, f"({next(letters)}) {length} km", corner=corner)
                pl.decimal_comma(a)
        pl.legend_above(fig, ncol=3)
        out.parent.mkdir(parents=True, exist_ok=True)
        pdf = pl.save_figure(fig, out, png=True)
    plt.close(fig)
    print(f"figura: {pdf} (e o PNG ao lado)")

def _snr_by_length(files) -> dict[int, dict[float, float]]:
    out: dict[int, dict[float, float]] = {}
    for f in files:
        for q in json.loads(Path(f).read_text(encoding="utf-8"))["points"]:
            out.setdefault(int(q["length_km"]), {})[float(q["power_dbm"])] = float(q["snr_db"])
    return out

def figure_snr(vpi_files: list[Path], sim_files: list[Path], out: Path, ylim=(12.5, 16.5),
               xlim=(5.5, 18.5)):
    """(a) SNR contra a potencia, um comprimento por cor (viridis ate 85 % da escala, um marcador
    por curva, guia de figuras): VPI em marcadores vazados, simcopy em linha. Os pontos em que o VPI
    perde a fase ficam abaixo da escala e aparecem como triangulos na borda de baixo.
    (b) SNR do simcopy menos a do VPI, onde o VPI converge (CONVERGED), com a faixa de +-0,1 dB."""
    mpl, plt, pl = _matplotlib()
    from matplotlib.lines import Line2D
    vpi, sim = _snr_by_length(vpi_files), _snr_by_length(sim_files)
    lengths = sorted(set(vpi) & set(sim))
    markers = "osD^vp<>"
    with mpl.rc_context(pl.FIGURE_STYLE):
        fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.2))
        cmap = plt.get_cmap("viridis")
        handles = []
        for k, length in enumerate(lengths):
            color = cmap(0.85 * k / max(len(lengths) - 1, 1))
            mk = markers[k % len(markers)]
            ps = sorted(sim[length])
            a.plot(ps, [sim[length][p] for p in ps], "-", marker=mk, ms=3, color=color)
            pv = sorted(vpi[length])
            shown = [p for p in pv if vpi[length][p] >= ylim[0]]
            a.plot(shown, [vpi[length][p] for p in shown], ls="", marker=mk, ms=7, mfc="none", mew=1.1,
                   color=color)
            lost = [p for p in pv if vpi[length][p] < ylim[0]]
            a.plot(lost, [ylim[0] + 0.03 * (ylim[1] - ylim[0])] * len(lost), ls="", marker="v", ms=6,
                   color=color)
            lo, hi = CONVERGED.get(length, (min(pv), max(pv)))
            common = [p for p in ps if p in vpi[length] and lo <= p <= hi]
            b.plot(common, [sim[length][p] - vpi[length][p] for p in common], "-", marker=mk, ms=4,
                   color=color)
            handles.append((Line2D([], [], color=color, marker=mk, ms=5), f"{length} km"))
        a.set(xlabel=POWER_LABEL, ylabel="SNR (dB)", ylim=ylim, xlim=xlim)
        b.axhline(0.0, color="black", lw=0.8)
        for v in (-0.1, 0.1):
            b.axhline(v, color="black", ls="--", lw=0.8)
        b.set(xlabel=POWER_LABEL, ylabel="SNR do simcopy menos a do VPI (dB)", ylim=(-0.3, 0.3))
        pl.panel_label(a, "(a)", corner="upper left")
        pl.panel_label(b, "(b) onde o VPI converge", corner="upper left")
        handles += [(Line2D([], [], ls="", marker="o", mfc="none", ms=7, color="black"), "VPI"),
                    (Line2D([], [], marker="o", ms=3, color="black"), "simcopy"),
                    (Line2D([], [], ls="", marker="v", ms=6, color="black"), "VPI abaixo da escala (fase perdida)")]
        pl.legend_above(fig, [h for h, _ in handles], [t for _, t in handles], ncol=4)
        for ax in (a, b):
            pl.decimal_comma(ax)
        out.parent.mkdir(parents=True, exist_ok=True)
        pdf = pl.save_figure(fig, out, png=True)
    plt.close(fig)
    print(f"figura: {pdf} (e o PNG ao lado)")

def _unit_grid(x: np.ndarray, sent: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Simbolos recebidos divididos pelo ganho complexo, na escala de mapping.constellation
    (niveis +-1/3 e +-1), e os transmitidos na mesma escala."""
    a = sent / np.max(np.abs(sent.real))
    return x / (np.sum(x * np.conj(a)) / np.sum(np.abs(a) ** 2)), a

def figure_constellations(vpi_folder: Path, symbols_folder: Path, out: Path, length: int = 125,
                          powers=(6, 14, 17), mode: str = "blind", window: int = 43, lim: float = 1.5):
    """Constelacoes da polarizacao X, VPI em cima e simcopy embaixo, como mapa de densidade com a
    mesma escala de cores. Do VPI entra uma amostra sorteada do arquivo, com tantos simbolos
    quantos o simcopy tem, para as densidades serem comparaveis. A SNR do rotulo e a do arquivo
    inteiro."""
    mpl, plt, pl = _matplotlib()
    from simcopy.kernels.mapping import constellation
    base = vpi_folder / "2026_SIMULATION" / "dados_4_canais"
    sent_vpi = decide_16qam(load_vpi(base / "dados_ideais" / "dados_ref_canal2" / "dados_4_canais_ref_X_pol_canal2.csv"))
    sub, pattern = VPI_SETS[length]
    rows = []
    for p in powers:
        d = np.load(symbols_folder / f"{length}km_{p}dBm_{mode}_janela{window}_semente0.npz")
        y_sim, a_sim = _unit_grid(d["received"][0], d["sent"][0])
        n = y_sim.size
        x_vpi = load_vpi(base / sub / pattern.format(pol="X", p=p))
        # amostra sorteada do arquivo inteiro: os trechos em que a fase do VPI se perde ficam no
        # meio do bloco, e os primeiros simbolos nao os mostrariam
        pick = np.sort(np.random.default_rng(0).choice(x_vpi.size, n, replace=False))
        y_vpi, a_vpi = _unit_grid(x_vpi[pick], sent_vpi[pick])
        rows.append((p, ("VPI", y_vpi, decompose(x_vpi, sent_vpi)["snr_db"]),
                     ("simcopy", y_sim, decompose(d["received"][0], d["sent"][0])["snr_db"])))
    ideal = constellation("16QAM")
    bins = int(np.clip(np.sqrt(rows[0][1][1].size), 60, 250))
    top = max(np.histogram2d(y.real, y.imag, bins=bins, range=[[-lim, lim], [-lim, lim]])[0].max()
              for _, *pair in rows for _, y, _ in pair)
    with mpl.rc_context(pl.FIGURE_STYLE):
        fig, axes = plt.subplots(2, len(powers), figsize=(3.1 * len(powers) + 0.9, 6.6), squeeze=False)
        letters = iter("abcdefghijklmnop")
        for r in range(2):
            for c, (p, *pair) in enumerate(rows):
                source, y, snr = pair[r]
                ax = axes[r, c]
                mesh = pl.constellation_density(ax, y, ideal, lim, bins, vmax=top)
                ax.set_ylim(-lim - 0.22 * 2 * lim, lim)
                pl.panel_label(ax, f"({next(letters)}) {source}, {pl.number(p, 0)} dBm\n"
                                   f"SNR = {pl.number(snr)} dB")
                if c:
                    ax.set(ylabel="")
                if r == 0:
                    ax.set(xlabel="")
                pl.decimal_comma(ax)
        fig.colorbar(mesh, ax=list(axes.ravel()), label="símbolos por célula", shrink=0.7)
        out.parent.mkdir(parents=True, exist_ok=True)
        pdf = pl.save_figure(fig, out, png=True)
    plt.close(fig)
    print(f"figura: {pdf} (e o PNG ao lado)")

# --------------------------------------------------------------------------- placar ----------
# Janela de potencias [dBm] em que o DSP do VPI converge, por comprimento [km]: e onde os ajustes
# e as comparacoes ponto a ponto sao feitos (docs/comparacao_vpi.md, secoes 4 e 10).
CONVERGED = {120: (7, 16), 125: (8, 16), 130: (9, 17), 135: (10, 17), 140: (11, 17)}
# Alvos de docs/plano_replicar_vpi.md, secao 1: (rotulo, unidade, tolerancia)
TARGETS = [("imagem", "%", 2.0), ("deslocamento", "%", 2.0), ("A", "dB", 0.1), ("B", "dB", 0.1),
           ("resto do piso", "dB", 0.15), ("C", "dB", 0.15), ("SNR ponto a ponto", "dB", 0.1),
           ("potência ótima", "dBm", 0.0), ("SER sobre a gaussiana", "%", 10.0)]

def _ser_gaussian_16qam(snr_db: float) -> float:
    from scipy.special import erfc
    q = 1.5 * 0.5 * erfc(np.sqrt(10 ** (snr_db / 10) / 10.0))
    return float(1.0 - (1.0 - q) ** 2)

def scoreboard(vpi_points: list[dict], sim_points: list[dict]) -> dict:
    """Linhas do placar por comprimento: simcopy menos VPI em cada grandeza (ou a faixa, para as
    que sao ponto a ponto), na janela em que o VPI converge."""
    out = {}
    key = lambda q: (int(q["length_km"]), float(q["power_dbm"]))
    v_all, s_all = {key(q): q for q in vpi_points}, {key(q): q for q in sim_points}
    for length in sorted({k[0] for k in v_all} & {k[0] for k in s_all} & set(CONVERGED)):
        lo, hi = CONVERGED[length]
        common = sorted(k for k in v_all if k in s_all and k[0] == length and lo <= k[1] <= hi)
        if len(common) < 4:
            continue
        v, s = [v_all[k] for k in common], [s_all[k] for k in common]
        rel = lambda name: [100.0 * (b[name] / a[name] - 1.0) for a, b in zip(v, s)]
        dsnr = [b["snr_db"] - a["snr_db"] for a, b in zip(v, s)]
        excess = lambda q: q["ser"] / _ser_gaussian_16qam(q["snr_db"])
        ser = [100.0 * (excess(b) / excess(a) - 1.0) for a, b in zip(v, s)]
        row = {"imagem": (min(rel("image")), max(rel("image"))),
               "deslocamento": (min(rel("offset")), max(rel("offset"))),
               "SNR ponto a ponto": (min(dsnr), max(dsnr)),
               "SER sobre a gaussiana": (min(ser), max(ser)),
               "potência ótima": (v[int(np.argmax([q["snr_db"] for q in v]))]["power_dbm"],
                                  s[int(np.argmax([q["snr_db"] for q in s]))]["power_dbm"]),
               "pontos": len(common)}
        for name, snr_key, term in (("A", "snr_db", "A"), ("B", "snr_db", "B"), ("C", "snr_db", "C"),
                                    ("resto do piso", "snr_residual_db", "B")):
            fv = fit_abc([v_all[k] for k in v_all if k[0] == length], length, snr_key, (lo, hi))
            fs = fit_abc([s_all[k] for k in s_all if k[0] == length], length, snr_key, (lo, hi))
            row[name] = 10.0 * np.log10(fs[term] / fv[term]) if fv and fs and fv[term] > 0 and fs[term] > 0 else None
        out[length] = row
    return out

def print_scoreboard(board: dict, label: str) -> None:
    lengths = sorted(board)
    print(f"### Placar: {label}")
    print("| Grandeza (simcopy menos VPI) | " + " | ".join(f"{L} km" for L in lengths) + " | Alvo | Estado |")
    print("|---|" + "---|" * (len(lengths) + 2))
    for name, unit, tol in TARGETS:
        cells, ok = [], True
        for L in lengths:
            value = board[L][name]
            if value is None:
                cells.append("sem ajuste")
                ok = False
            elif name == "potência ótima":
                cells.append(f"{value[0]:g} e {value[1]:g}")
                ok &= value[0] == value[1]
            elif isinstance(value, tuple):
                cells.append(f"{value[0]:+.2f} a {value[1]:+.2f}" if unit == "dB" else f"{value[0]:+.0f} a {value[1]:+.0f}")
                ok &= max(abs(value[0]), abs(value[1])) <= tol
            else:
                cells.append(f"{value:+.2f}")
                ok &= abs(value) <= tol
        target = "igual" if name == "potência ótima" else f"±{tol:g} {unit}"
        print(f"| {name} [{unit}] | " + " | ".join(cells) + f" | {target} | {'dentro' if ok else 'fora'} |")
    print()

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    s = sub.add_parser("simular")
    s.add_argument("--config", type=Path, default=ROOT / "examples" / "400zr_4ch_75ghz.yaml")
    s.add_argument("--lengths", type=int, nargs="+", default=[125])
    s.add_argument("--powers", type=float, nargs="+", default=[8, 10, 12, 13, 14, 15, 16, 17])
    s.add_argument("--n-symbols", type=int, default=16384)
    s.add_argument("--sps", type=int, default=16)
    s.add_argument("--power-is", choices=["total", "per_channel"], default="total")
    s.add_argument("--channel", type=int, default=2)
    s.add_argument("--no-pmd", action="store_true")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--dsp", nargs="*", default=[], metavar="CAMPO=VALOR",
                   help="campos de rx.dsp, ex.: mode=blind phase_window=81 matched_filter=false")
    s.add_argument("--second-dsp", nargs="*", action="append", default=None, metavar="CAMPO=VALOR",
                   help="outra cadeia de DSP medida no mesmo sinal propagado, ex.: mode=blind. Pode "
                        "ser repetida, com um --second-out para cada uma")
    s.add_argument("--second-out", type=Path, action="append", help="arquivo de cada cadeia a mais")
    s.add_argument("--save-signals", type=Path, metavar="PASTA",
                   help="grava o sinal na saida do ADC de cada ponto (.npz), antes do DSP cego")
    s.add_argument("--save-symbols", type=Path, metavar="PASTA",
                   help="grava os simbolos de cada ponto (.npz) nesta pasta")
    s.add_argument("-o", "--out", type=Path, required=True)
    v = sub.add_parser("vpi")
    v.add_argument("folder", type=Path)
    v.add_argument("--lengths", type=int, nargs="+", default=[125, 130, 135, 140])
    v.add_argument("--powers", type=int, nargs="+", default=list(range(0, 21)))
    v.add_argument("-o", "--out", type=Path, required=True)
    f = sub.add_parser("figura")
    f.add_argument("files", type=Path, nargs="+")
    f.add_argument("--lengths", type=int, nargs="+")
    f.add_argument("-o", "--out", type=Path, required=True)
    fs = sub.add_parser("figura-snr", help="SNR contra a potencia, VPI e simcopy, e a diferenca")
    fs.add_argument("--vpi", type=Path, nargs="+", required=True, help="resumos do VPI (subcomando vpi)")
    fs.add_argument("--simcopy", type=Path, nargs="+", required=True, help="resultados do simcopy")
    fs.add_argument("-o", "--out", type=Path, required=True)
    fc = sub.add_parser("figura-constelacoes", help="constelacoes do VPI e do simcopy em mapa de densidade")
    fc.add_argument("vpi_folder", type=Path, help="pasta extraido dos dados do VPI")
    fc.add_argument("symbols_folder", type=Path, help="pasta de simular --save-symbols")
    fc.add_argument("--length", type=int, default=125)
    fc.add_argument("--powers", type=int, nargs="+", default=[6, 14, 17])
    fc.add_argument("--mode", default="blind")
    fc.add_argument("--window", type=int, default=43)
    fc.add_argument("-o", "--out", type=Path, required=True)
    t = sub.add_parser("ajuste")
    t.add_argument("files", type=Path, nargs="+")
    t.add_argument("--window", type=float, nargs=2, default=[8, 16])
    b = sub.add_parser("placar")
    b.add_argument("--vpi", type=Path, nargs="+", required=True, help="resumos do VPI (subcomando vpi)")
    b.add_argument("--simcopy", type=Path, nargs="+", required=True, help="resultados do simcopy")
    a = ap.parse_args(argv)

    if a.command == "placar":
        vpi_points = [q for f in a.vpi for q in json.loads(f.read_text(encoding="utf-8"))["points"]]
        for f in a.simcopy:
            data = json.loads(f.read_text(encoding="utf-8"))
            print_scoreboard(scoreboard(vpi_points, data["points"]), _label(data, f) + f" ({f.name})")
        return

    if a.command == "figura":
        return figure(a.files, a.out, a.lengths)
    if a.command == "figura-snr":
        return figure_snr(a.vpi, a.simcopy, a.out)
    if a.command == "figura-constelacoes":
        return figure_constellations(a.vpi_folder, a.symbols_folder, a.out, a.length, a.powers,
                                     a.mode, a.window)
    if a.command == "ajuste":
        for file in a.files:
            data = json.loads(file.read_text(encoding="utf-8"))
            for length in sorted({q["length_km"] for q in data["points"]}):
                for key in ("snr_db", "snr_residual_db"):
                    fit = fit_abc(data["points"], length, key, tuple(a.window))
                    if fit:
                        print(f"{file.stem} {length} km [{key}]: A = {fit['A']:.4f} | B = {fit['B']:.5f} "
                              f"(piso {fit['floor_db'] if fit['floor_db'] is None else round(fit['floor_db'], 2)} dB) | "
                              f"C = {fit['C']:.3e} | resíduo rms = {fit['rms_db']:.3f} dB | {fit['n']} pontos")
        return
    if a.command == "vpi":
        result = analyse_vpi(a.folder, a.lengths, a.powers)
    else:
        def parse(items) -> dict:
            options = {}
            for item in items:
                key, value = item.split("=", 1)
                if value.lower() in ("true", "false"):
                    options[key] = value.lower() == "true"
                else:
                    for cast in (int, float, str):
                        try:
                            options[key] = cast(value)
                            break
                        except ValueError:
                            pass
            return options
        second = []
        if a.second_dsp is not None:
            if a.second_out is None or len(a.second_out) != len(a.second_dsp):
                ap.error("cada --second-dsp pede o seu --second-out")
            second = [(parse(items), out) for items, out in zip(a.second_dsp, a.second_out)]
        result = simulate(a.config, a.lengths, a.powers, a.n_symbols, a.sps, a.power_is, a.channel,
                          not a.no_pmd, a.seed, parse(a.dsp), checkpoint=a.out, second=second,
                          save_symbols=a.save_symbols, save_signals=a.save_signals)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=1), encoding="utf-8")
    print(f"gravado: {a.out}")

if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
