"""Graficos de interesse (docs/mapa_arquitetura.md, 5.12; docs/sistema_imdd.md, secao 5).

Uma funcao serve a qualquer componente porque os paineis sao escolhidos pelo dominio do sinal,
nao pelo componente:
- sinal eletrico: forma de onda no tempo e espectro;
- sinal optico: potencia e chirp no tempo e espectro optico;
- bits: degraus no tempo.
Os blocos nao plotam (CLAUDE.md, regra 6): estes graficos leem as saidas que o executor guarda
(Result.output) e as contas de simcopy.analysis. Precisa do matplotlib (pip install -e ".[figures]").

    from simcopy.plotting import plot_stage, plot_chain, plot_eye, plot_sweep
    from simcopy.plotting import plot_constellation      # canal coerente
    plot_stage(r, "link.span0.fiber")     # entrada e saida da fibra, tempo e frequencia

Estilo (guia de figuras do Gabriel, TCC/codigos/guia_figuras.md, secoes 1 e 2): fonte Times,
matematica em stix, ticks para dentro nos quatro lados, grade tracejada, sem titulo dentro do
grafico (o que era titulo virou rotulo de painel, "(a) ...", com fundo branco), legenda sem
moldura acima dos eixos, cores primarias nomeadas, virgula decimal, constelacoes como mapa de
densidade sem grade. Tamanho: o dos notebooks, com as medidas de fonte e linha do relatorio.
Cada funcao aplica o estilo so nas figuras que cria; `apply_style()` aplica nas demais (por
exemplo, as figuras feitas direto nos notebooks). Tamanhos por documento e idioma:
docs/proposta_figuras.md.
"""
from __future__ import annotations
import functools
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
from . import analysis
from .core.signals import Bits, DigitalSignal, Domain, Waveform
from .kernels.demodulation import sample_symbols
from .kernels.mapping import constellation

# Uma cor segue sempre a mesma coisa em todas as figuras (guia, secao 2).
COLORS = {"entrada": "blue", "saída": "red", "I": "blue", "Q": "red", "bits": "black",
          "ideal": "red", "limiar": "gray", "VPI": "black", "simcopy": "blue"}
IN_COLOR, OUT_COLOR = COLORS["entrada"], COLORS["saída"]
DENSITY_CMAP = "viridis"

# rcParams do guia (secao 1), com a fonte 2 pt acima da do relatorio (11, 10 e 9,5), porque as
# figuras dos notebooks sao mais largas que os 6,3 pol do relatorio (pedido de 2026-11-18).
# A familia vai pelo nome, e nao por "serif": o nome generico so e resolvido na hora de desenhar,
# que num notebook acontece depois que a funcao de grafico ja saiu do contexto do estilo. DejaVu
# Sans fica no fim para os simbolos que a Times nao tem (por exemplo a estrela).
FIGURE_STYLE = {
    "font.family": ["Times New Roman", "DejaVu Serif", "DejaVu Sans"],
    "mathtext.fontset": "stix",
    "font.size": 13, "axes.labelsize": 13, "axes.titlesize": 13,
    "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 11.5,
    "axes.linewidth": 1.0,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 5, "ytick.major.size": 5, "xtick.major.width": 1.0, "ytick.major.width": 1.0,
    "xtick.minor.size": 3, "ytick.minor.size": 3, "xtick.minor.width": 0.8, "ytick.minor.width": 0.8,
    "axes.grid": True, "axes.grid.which": "both",
    "grid.linestyle": "--", "grid.linewidth": 0.4, "grid.color": "#b0b0b0",
    "axes.axisbelow": True,
    "axes.prop_cycle": cycler(color=["blue", "red", "green", "orange", "purple", "black", "gray"]),
    "lines.linewidth": 1.6, "lines.markersize": 5,
    "legend.frameon": False,
    "figure.constrained_layout.use": True,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
    "pdf.fonttype": 42,
}


def apply_style(comma: bool = True, inline: str | None = "svg") -> None:
    """Aplica o estilo a todas as figuras seguintes (rcParams globais): e o que os notebooks
    chamam para as figuras feitas direto nas celulas.

    Sem o layout automatico (constrained layout): as celulas chamam tight_layout, e os dois nao
    convivem numa figura com barra de cores. Com comma, os eixos lineares saem com virgula decimal,
    pelo formato numerico do portugues (so a categoria numerica do locale, LC_NUMERIC; o str() e o
    format() do Python nao mudam). inline: formato das figuras mostradas no notebook. O notebook
    nao mostra PDF; o SVG e o formato vetorial que ele mostra. None deixa o padrao (PNG)."""
    mpl.rcParams.update({k: v for k, v in FIGURE_STYLE.items() if k != "figure.constrained_layout.use"})
    if inline:
        try:
            from IPython import get_ipython
            if get_ipython() is not None:
                from matplotlib_inline.backend_inline import set_matplotlib_formats
                set_matplotlib_formats(inline)
        except ImportError:
            pass
    if comma:
        import locale
        for name in ("pt_BR.UTF-8", "pt_BR.utf8", "pt_BR", "Portuguese_Brazil.1252"):
            try:
                locale.setlocale(locale.LC_NUMERIC, name)
            except locale.Error:
                continue
            mpl.rcParams["axes.formatter.use_locale"] = True
            break


def save_figure(fig, path, png: bool = False, dpi: int = 150):
    """Grava a figura em PDF vetorial, com a fonte embutida (guia, secao 2), no caminho dado com a
    extensao trocada para .pdf. png=True grava tambem um PNG ao lado, para os documentos em
    Markdown, que nao mostram PDF. Devolve o caminho do PDF."""
    from pathlib import Path
    pdf = Path(path).with_suffix(".pdf")
    pdf.parent.mkdir(parents=True, exist_ok=True)
    with mpl.rc_context({k: FIGURE_STYLE[k] for k in ("savefig.bbox", "savefig.pad_inches", "pdf.fonttype")}):
        fig.savefig(pdf)
        if png:
            fig.savefig(pdf.with_suffix(".png"), dpi=dpi)
    return pdf


def _tick(v: float) -> str:
    """Rotulo de tick em portugues: virgula decimal, ponto no milhar, sinal de menos."""
    if float(v).is_integer() and abs(v) >= 1000:
        text = f"{int(v):,}".replace(",", ".")
    else:
        text = f"{v:g}".replace(".", ",")
    return text.replace("-", "−")


def decimal_comma(ax, axes: str = "xy") -> None:
    """Virgula decimal, ponto no milhar e sinal de menos nos eixos lineares (guia, secao 2)."""
    fmt = FuncFormatter(lambda v, _: _tick(v))
    for e in axes:
        axis = ax.xaxis if e == "x" else ax.yaxis
        if axis.get_scale() == "linear":
            axis.set_major_formatter(fmt)


def panel_label(ax, text: str, corner: str = "lower left"):
    """Rotulo de painel dentro do eixo, com fundo branco: "(a)" ou "(a) 120 km"."""
    v, h = corner.split()
    x, ha = (0.03, "left") if h == "left" else (0.97, "right")
    y, va = (0.04, "bottom") if v == "lower" else (0.96, "top")
    return ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, zorder=10,
                   bbox=dict(facecolor="white", edgecolor="none", pad=1.5))


def legend_above(fig_or_ax, handles=None, labels=None, ncol: int = 4):
    """Legenda sem moldura acima dos eixos. Numa figura (varios paineis), uma legenda so, sem
    repetir rotulos."""
    if isinstance(fig_or_ax, mpl.figure.Figure):
        if handles is None:
            seen = {}
            for ax in fig_or_ax.axes:
                for hd, lb in zip(*ax.get_legend_handles_labels()):
                    seen.setdefault(lb, hd)
            labels, handles = list(seen), list(seen.values())
        if not handles:
            return None
        return fig_or_ax.legend(handles, labels, loc="outside upper center",
                                ncol=min(ncol, len(handles)))
    kw = {} if handles is None else dict(handles=handles, labels=labels)
    return fig_or_ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=ncol, **kw)


def number(x: float, digits: int = 2) -> str:
    """Numero com virgula decimal."""
    return f"{x:.{digits}f}".replace(".", ",").replace("-", "−")


def scientific(x: float, digits: int = 1) -> str:
    """1,2 x 10^-3, com virgula decimal."""
    if x == 0 or not np.isfinite(x):
        return number(x, 0) if x == 0 else str(x)
    mantissa, exponent = f"{x:.{digits}e}".split("e")
    return f"{mantissa.replace('.', ',')}" + r"$\,\times\,10^{" + str(int(exponent)) + "}$"


def _bottom_room(ax, fraction: float = 0.22) -> None:
    """Espaco embaixo dos dados para o rotulo do painel."""
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo - fraction * (hi - lo), hi)


def _top_room(ax, fraction: float = 0.3) -> None:
    """Espaco em cima dos dados, para o rotulo do painel onde embaixo nao cabe (histogramas)."""
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + fraction * (hi - lo))


def _styled(fn):
    """Cria a figura dentro do estilo e poe virgula decimal em todos os eixos lineares."""
    @functools.wraps(fn)
    def wrapper(*args, **kw):
        with mpl.rc_context(FIGURE_STYLE):
            fig = fn(*args, **kw)
            for ax in fig.axes:
                decimal_comma(ax)
        return fig
    return wrapper


def light_colormap(name: str, start: float = 0.3):
    """Mapa de cores sequencial do matplotlib sem a ponta quase branca, para a celula de um ponto
    so nao sumir no fundo branco."""
    from matplotlib.colors import ListedColormap
    return ListedColormap(mpl.colormaps[name](np.linspace(start, 1.0, 256)), name=f"{name}_claro")


POLARIZATION_CMAPS = ("Blues", "Reds")          # X em azul, Y em vermelho


def polarization_overlay(ax, symbols, ideal, angle_deg: float, lim: float, bins: int | None = None,
                         vmax: float | None = None):
    """As duas polarizacoes no mesmo plano: X em azul e Y em vermelho, cada uma com o seu mapa de
    densidade, e Y girada de angle_deg so para as duas nao cairem uma sobre a outra (o giro e do
    desenho, nao do sinal). Pontos ideais das duas em preto, sem os limiares. symbols: (2, n).
    Devolve os dois mapas, para as barras de cores."""
    meshes = []
    turn = np.exp(1j * np.radians(angle_deg))
    for p, (row, cmap) in enumerate(zip(np.asarray(symbols), POLARIZATION_CMAPS)):
        rot = turn if p else 1.0
        meshes.append(density_map(ax, row * rot, lim, bins, vmax, cmap=light_colormap(cmap)))
        points = np.asarray(ideal) * rot
        ax.plot(points.real, points.imag, "+", color="black", ms=6, mew=1.0)
    ax.set(xlabel="em fase (I)", ylabel="quadratura (Q)")
    return meshes


def density_map(ax, points, lim: float | None = None, bins: int | None = None,
                vmax: float | None = None, cmap=None):
    """Mapa de densidade de pontos complexos (histograma 2D, escala logaritmica, celulas vazias em
    branco), num quadrado de lado 2*lim, sem grade. Serve a qualquer nuvem de pontos no plano I-Q,
    em qualquer unidade. lim: padrao 1,1 vezes o percentil 99,9 de |I| e |Q|. bins: padrao a raiz
    do numero de pontos, entre 60 e 250. cmap: padrao viridis. Devolve o mapa, para a barra de
    cores."""
    y = np.ravel(np.asarray(points))
    if lim is None:
        lim = 1.1 * float(np.percentile(np.abs(np.concatenate([y.real, y.imag])), 99.9))
    if bins is None:
        bins = int(np.clip(np.sqrt(y.size), 60, 250))
    h, xe, ye = np.histogram2d(y.real, y.imag, bins=bins, range=[[-lim, lim], [-lim, lim]])
    h = np.ma.masked_equal(h, 0)
    mesh = ax.pcolormesh(xe, ye, h.T, cmap=cmap or DENSITY_CMAP, rasterized=True,
                         norm=LogNorm(vmin=1, vmax=vmax or max(float(h.max()), 2.0)))
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), aspect="equal")
    ax.grid(False)
    return mesh


def constellation_density(ax, symbols, ideal, lim: float, bins: int | None = None,
                          vmax: float | None = None):
    """Mapa de densidade dos simbolos recebidos (density_map), com os limiares de decisao
    tracejados e os pontos ideais. Devolve o mapa, para a barra de cores."""
    mesh = density_map(ax, symbols, lim, bins, vmax)
    ideal = np.asarray(ideal)
    for levels, line in ((np.unique(np.round(ideal.real, 9)), ax.axvline),
                         (np.unique(np.round(ideal.imag, 9)), ax.axhline)):
        for t in (levels[:-1] + levels[1:]) / 2:
            line(t, color=COLORS["limiar"], ls="--", lw=0.7)
    ax.plot(ideal.real, ideal.imag, "+", color=COLORS["ideal"], ms=6, mew=1.0)
    ax.set(xlabel="em fase (I)", ylabel="quadratura (Q)", xlim=(-lim, lim), ylim=(-lim, lim),
           aspect="equal")
    ax.grid(False)
    return mesh


def _bits_trace(b: Bits) -> tuple[np.ndarray, np.ndarray]:
    # com mais de um bit por simbolo (QPSK), cada bit ocupa uma fracao do simbolo
    sps = b.geometry.samples_per_symbol // b.geometry.bits_per_symbol
    x = np.repeat(b.data[0].astype(float), sps)
    return np.arange(x.size), x

def _signals(pairs):
    """Separa (sinal, rotulo, cor, unidade) por dominio, ignorando None."""
    el, op, bt = [], [], []
    for s, lab, cor, un in pairs:
        if s is None:
            continue
        if isinstance(s, Bits):
            bt.append((s, lab, cor, un))
        elif s.domain is Domain.OPTICAL:
            op.append((s, lab, cor, un))
        else:
            el.append((s, lab, cor, un))
    return el, op, bt

@_styled
def plot_io(entrada, saida, labels=("entrada", "saída"), title: str = "",
            units=(None, None), time_window: float | None = None,
            freq_span: float | None = None, resolution: float = 0.5e9):
    """Entrada e saida de um componente, no tempo e na frequencia, com os paineis escolhidos pelo
    dominio de cada sinal (Waveform eletrico ou optico, Bits ou None).

    time_window: duracao mostrada no tempo [s] (padrao: 32 simbolos, ou a janela se for menor).
    freq_span: faixa de frequencia mostrada, +-freq_span no optico e 0..freq_span no eletrico [Hz].
    resolution: media movel dos espectros [Hz] (analysis.optical_spectrum).
    title: vai no rotulo do primeiro painel (o guia nao poe titulo dentro do grafico).
    """
    pairs = list(zip((entrada, saida), labels, (IN_COLOR, OUT_COLOR), units))
    el, op, bt = _signals(pairs)
    order = []
    for s, *_ in pairs:
        if s is None:
            continue
        kind = "bits" if isinstance(s, Bits) else s.domain.value
        if kind not in order:
            order.append(kind)
    n_time = sum({"bits": 1, "electrical": 1, "optical": 2}[k] for k in order)
    n_freq = sum(1 for k in order if k != "bits")
    fig, axes = plt.subplots(n_time + n_freq, 1, figsize=(9, 2.0 * (n_time + n_freq)))
    axes = list(np.atleast_1d(axes))
    panels = list(axes)
    ref = next(s for s, *_ in pairs if isinstance(s, Waveform))
    t_ns = ref.grid.time * 1e9
    if time_window is None:
        time_window = min(ref.grid.duration, 32 * ref.geometry.samples_per_symbol * ref.grid.dt)
    t_max = time_window * 1e9

    for kind in order:
        if kind == "bits":
            a = axes.pop(0)
            for b, lab, cor, _ in bt:
                idx, x = _bits_trace(b)
                a.step(t_ns[idx % t_ns.size], x, where="post", color=cor, label=lab)
            a.set(xlabel="t (ns)", ylabel="bit", xlim=(0, t_max))
            _bottom_room(a)
        elif kind == "electrical":
            a = axes.pop(0)
            for w, lab, cor, un in el:
                a.plot(t_ns, w.samples[0].real, color=cor, label=f"{lab} ({un or 'u.a.'})")
            a.set(xlabel="t (ns)", ylabel="sinal elétrico", xlim=(0, t_max))
            _bottom_room(a)
        else:
            a = axes.pop(0)
            for w, lab, cor, _ in op:
                a.plot(t_ns, analysis.optical_power(w) * 1e3, color=cor, label=lab)
            a.set(xlabel="t (ns)", ylabel="P (mW)", xlim=(0, t_max))
            _bottom_room(a)
            a = axes.pop(0)
            for j, (w, lab, cor, _) in enumerate(op):
                a.plot(t_ns, analysis.chirp(w) / 1e9, color=cor, ls="-" if j == 0 else "--")
            a.set(xlabel="t (ns)", ylabel=r"chirp $\delta\nu$ (GHz)", xlim=(0, t_max))
            _bottom_room(a)
    for kind in order:
        if kind == "electrical":
            a = axes.pop(0)
            for w, lab, cor, _ in el:
                f, s = analysis.electrical_spectrum(w, resolution)
                a.plot(f / 1e9, 10 * np.log10(s / s.max() + 1e-15), color=cor)
            span = (freq_span or f.max()) / 1e9
            a.set(xlabel="f (GHz)", ylabel="espectro elétrico (dB)", ylim=(-70, 3), xlim=(0, span))
        elif kind == "optical":
            a = axes.pop(0)
            top = max(analysis.optical_spectrum(w, resolution)[1].max() for w, *_ in op)
            for j, (w, lab, cor, _) in enumerate(op):
                f, s = analysis.optical_spectrum(w, resolution)
                a.plot(f / 1e9, 10 * np.log10(s / top + 1e-15), color=cor,
                       ls="-" if j == 0 else "--")
            span = (freq_span or abs(f).max()) / 1e9
            a.set(xlabel="frequência óptica relativa ao centro da banda (GHz), positivo = azul",
                  ylabel="espectro óptico (dB)", ylim=(-70, 3), xlim=(-span, span))
    for k, a in enumerate(panels):
        panel_label(a, f"({'abcdefgh'[k]})" + (f" {title}" if title and k == 0 else ""))
    legend_above(fig, ncol=2)
    return fig

def _unit(block):
    return getattr(block, "output_unit", None)

def plot_stage(result, stage_id: str, **kw):
    """plot_io da entrada e da saida de um estagio do resultado (Result.ids lista os estagios).
    No mux, a entrada mostrada e a do primeiro canal."""
    st = result.stage(stage_id)
    prev = result.stage(st.inputs[0]) if st.inputs else None
    if isinstance(st.output, DigitalSignal) or (prev and isinstance(prev.output, DigitalSignal)):
        raise ValueError(f"{stage_id}: estagio na taxa do DSP (DigitalSignal). Use plot_chain ou "
                         f"plot_constellation.")
    entrada = prev.output if prev is not None else None
    labels = (f"entrada ({prev.id})" if prev else "entrada", f"saída ({st.id})")
    units = (_unit(prev.block) if prev else None, _unit(st.block))
    return plot_io(entrada, st.output, labels, title=f"{st.id}: {type(st.block).__name__}",
                   units=units, **kw)

@_styled
def plot_chain(result, time_window: float | None = 2e-9, channel: int = 0):
    """Visao da cadeia no tempo: um painel por estagio no caminho do canal (bits, sinal eletrico,
    potencia optica, corrente), como a celula do MZM no notebook 07."""
    stages = result.path(channel)
    fig, axes = plt.subplots(len(stages), 1, figsize=(9, 1.7 * len(stages)), sharex=True)
    axes = np.atleast_1d(axes)
    has_iq = False
    for k, (a, st) in enumerate(zip(axes, stages)):
        s = st.output
        if isinstance(s, Bits):
            ref = next(x.output for x in stages if isinstance(x.output, Waveform))
            idx, x = _bits_trace(s)
            a.step(ref.grid.time[idx % ref.grid.n_samples] * 1e9, x, where="post",
                   color=COLORS["bits"])
            a.set(ylabel="bit")
        elif isinstance(s, DigitalSignal):       # taxa do DSP: pontos, na polarizacao X
            t = np.arange(s.samples.shape[1]) / s.sampling_rate * 1e9
            a.plot(t, s.samples[0].real, ".", ms=3, color=COLORS["I"])
            a.plot(t, s.samples[0].imag, ".", ms=3, color=COLORS["Q"])
            has_iq = True
            a.set(ylabel=_unit(st.block) or "u.a.")
        elif s.domain is Domain.OPTICAL:
            a.plot(s.grid.time * 1e9, analysis.optical_power(s) * 1e3, color=OUT_COLOR)
            a.set(ylabel="P (mW)")
        else:
            un = _unit(st.block)
            y, lab = s.samples[0].real, un or "u.a."
            scale = 1e3 if un == "A" else 1.0
            if un == "A":
                y, lab = y * 1e3, "mA"
            a.plot(s.grid.time * 1e9, y, color=COLORS["I"])
            if np.any(s.samples[0].imag):          # sinal I + i*Q do coerente
                a.plot(s.grid.time * 1e9, s.samples[0].imag * scale, color=COLORS["Q"], ls="--")
                has_iq = True
            a.set(ylabel=lab)
        _bottom_room(a, 0.35)
        panel_label(a, f"({'abcdefghijklmnopqrstuvwxyz'[k % 26]}) {st.id} ({type(st.block).__name__})")
    axes[-1].set(xlabel="t (ns)")
    if time_window is not None:
        axes[-1].set_xlim(0, time_window * 1e9)
    if has_iq:
        legend_above(fig, [Line2D([], [], color=COLORS["I"]), Line2D([], [], color=COLORS["Q"], ls="--")],
                     ["em fase (I)", "quadratura (Q)"], ncol=2)
    return fig

@_styled
def plot_eye(result, channel: int = 0, max_traces: int = 400):
    """Olho da corrente no fotodetector, com o instante de decisao, os niveis I1 e I0, o limiar
    otimo e o fator Q medidos (kernels/decision.py)."""
    m = result.metrics[channel]
    i = result.output(f"ch{channel}.rx.photodiode").samples[0].real
    sps = result.resolved.channels[channel].samples_per_symbol
    if hasattr(m, "eye_q"):               # PAM de varios niveis (kernels/decision.py, measure_pam)
        return _plot_eye_pam(result, channel, m, i, sps, max_traces)
    last = min(m.last_symbol, m.first_symbol + max_traces)
    t, tr = analysis.eye(i, sps, m.delay, m.first_symbol, last)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), gridspec_kw={"width_ratios": [3, 1]})
    ax[0].plot(t, tr.T * 1e3, color=COLORS["entrada"], lw=0.4, alpha=0.35,
               rasterized=True)          # 400 tracos: raster dentro do PDF/SVG, eixos e texto vetoriais
    for y, lab, ls in ((m.level_one, "$I_1$", "-"), (m.level_zero, "$I_0$", "-"),
                       (m.threshold, "$I_D$", "--")):
        ax[0].axhline(y * 1e3, color="black", ls=ls, lw=0.8)
        ax[0].text(t[-1], y * 1e3, f" {lab}", va="center")
    ax[0].axvline(0, color="black", lw=0.6, ls=":")
    ax[0].set(xlabel="tempo (símbolos), 0 = instante de decisão", ylabel="corrente (mA)")
    bits = result.output(f"ch{channel}.tx.bits").data[0].astype(bool)
    k = np.arange(m.first_symbol, m.last_symbol)
    s = i[(k * sps + m.delay) % i.size]
    ax[1].hist(s[bits[k]] * 1e3, bins=40, orientation="horizontal", color=COLORS["saída"],
               alpha=0.6, label="bit 1")
    ax[1].hist(s[~bits[k]] * 1e3, bins=40, orientation="horizontal", color=COLORS["entrada"],
               alpha=0.6, label="bit 0")
    ax[1].axhline(m.threshold * 1e3, color="black", ls="--", lw=0.8)
    ax[1].set(xlabel="contagem")
    ax[1].sharey(ax[0])
    _bottom_room(ax[0], 0.3)
    panel_label(ax[0], f"(a) canal {channel}: Q = {number(m.q)}, BER por Q = {scientific(m.ber_q)}")
    panel_label(ax[1], "(b) na decisão")
    legend_above(fig, ncol=2)
    return fig

def _plot_eye_pam(result, channel, m, current, sps, max_traces):
    """Olho de um sinal PAM de varios niveis: os niveis, os limiares e o Q de cada olho."""
    from .kernels.mapping import gray_to_index
    c = next(x for x in result.system.channels if x.index == channel)
    k = c.modulation.bits_per_symbol
    last = min(m.last_symbol, m.first_symbol + max_traces)
    t, tr = analysis.eye(current, sps, m.delay, m.first_symbol, last)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"width_ratios": [3, 1]})
    ax[0].plot(t, tr.T * 1e3, color=COLORS["entrada"], lw=0.4, alpha=0.35,
               rasterized=True)          # 400 tracos: raster dentro do PDF/SVG, eixos e texto vetoriais
    for y in m.levels:
        ax[0].axhline(y * 1e3, color="black", lw=0.8)
    for y in m.thresholds:
        ax[0].axhline(y * 1e3, color="black", ls="--", lw=0.8)
    ax[0].axvline(0, color="black", lw=0.6, ls=":")
    ax[0].set(xlabel="tempo (símbolos), 0 = instante de decisão", ylabel="corrente (mA)")
    bits = result.output(f"ch{channel}.tx.bits").data[0].reshape(-1, k)
    index = gray_to_index(bits)
    sym = np.arange(m.first_symbol, m.last_symbol)
    s = current[(sym * sps + m.delay) % current.size]
    for level in range(2 ** k):
        ax[1].hist(s[index[sym] == level] * 1e3, bins=40, orientation="horizontal", alpha=0.6,
                   label=f"nível {level}")
    for y in m.thresholds:
        ax[1].axhline(y * 1e3, color="black", ls="--", lw=0.8)
    ax[1].set(xlabel="contagem")
    ax[1].sharey(ax[0])
    _bottom_room(ax[0], 0.3)
    eyes = "; ".join(number(q, 1) for q in m.eye_q)
    panel_label(ax[0], f"(a) canal {channel}: Q por olho = {eyes}, BER por Q = {scientific(m.ber_q)}")
    panel_label(ax[1], "(b) na decisão")
    legend_above(fig, ncol=2 ** k)
    return fig

@_styled
def plot_constellation(result, channel: int = 0, max_points: int | None = None):
    """Constelacao de um canal coerente, como mapa de densidade: a corrente complexa do receptor
    amostrada no instante de decisao e dividida pelo ganho complexo medido
    (kernels/demodulation.py), com os pontos ideais e os limiares de decisao. Ao lado, o
    histograma de cada quadratura. max_points: None usa todos os simbolos medidos."""
    m = result.metrics[channel]
    c = next(x for x in result.system.channels if x.index == channel)
    rc = next(x for x in result.resolved.channels if x.index == channel)
    if channel in result.symbols:        # canal com rx.dsp: os simbolos ja equalizados
        return _plot_equalized(result.symbols[channel], constellation(c.modulation.name), m,
                               channel, max_points)
    i = result.output(f"ch{channel}.rx.frontend").samples[0]
    last = m.last_symbol if max_points is None else min(m.last_symbol, m.first_symbol + max_points)
    y = sample_symbols(i, rc.samples_per_symbol, m.delay, m.first_symbol, last)
    y = y / (m.gain * np.exp(1j * m.phase))
    ideal = constellation(c.modulation.name)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"width_ratios": [1, 1.2]})
    lim = 1.15 * max(np.abs(y.real).max(), np.abs(y.imag).max(), 1.0)
    mesh = constellation_density(ax[0], y, ideal, lim)
    fig.colorbar(mesh, ax=ax[0], label="símbolos por célula", shrink=0.85)
    ax[1].hist(y.real, bins=80, color=COLORS["I"], alpha=0.6, label="em fase (I)")
    ax[1].hist(y.imag, bins=80, color=COLORS["Q"], alpha=0.5, label="quadratura (Q)")
    ax[1].axvline(0, color="black", ls="--", lw=0.8)
    ax[1].set(xlabel="valor normalizado, limiar em 0", ylabel="contagem")
    _top_room(ax[1], 0.35)
    _bottom_room(ax[0], 0.12)
    panel_label(ax[0], f"(a) canal {channel}: EVM = {number(100 * m.evm, 1)} %")
    panel_label(ax[1], f"(b) Q = {number(m.q)}, BER por Q = {scientific(m.ber_q)},\n"
                       f"contada = {scientific(m.ber_counted)}", corner="upper left")
    legend_above(fig, ncol=2)
    return fig

def _plot_equalized(symbols, ideal, m, channel: int, max_points: int | None):
    """Uma constelacao por polarizacao, dos simbolos na saida do DSP (Result.symbols), como mapa
    de densidade com a mesma escala de cores nas duas."""
    n_pol = symbols.shape[0]
    fig, axes = plt.subplots(1, n_pol, figsize=(4.6 * n_pol + 0.8, 4.6), squeeze=False)
    # limite dos eixos: 1,25 vezes o maior nivel, ou mais se as nuvens passarem disso (pouca SNR,
    # ou BPSK e QPSK, em que todos os pontos estao na borda)
    shown = symbols if max_points is None else symbols[:, :max_points]
    lim = max(1.25 * max(np.abs(ideal.real).max(), np.abs(ideal.imag).max()),
              1.05 * float(np.percentile(np.abs(np.concatenate([shown.real, shown.imag], axis=None)),
                                         99.9)))
    bins = int(np.clip(np.sqrt(shown.shape[1]), 60, 250))
    top = max(np.histogram2d(y.real, y.imag, bins=bins, range=[[-lim, lim], [-lim, lim]])[0].max()
              for y in shown)
    for p, a in enumerate(axes[0]):
        mesh = constellation_density(a, shown[p], ideal, lim, bins, vmax=top)
        _bottom_room(a, 0.2)
        panel_label(a, f"({'ab'[p]}) {'XY'[p]}: SNR = {number(m.snr_db_pol[p])} dB\n"
                       f"SER = {scientific(m.ser_pol[p])}, BER = {scientific(m.ber_pol[p])}")
        if p:
            a.set(ylabel="")
    fig.colorbar(mesh, ax=list(axes[0]), label="símbolos por célula", shrink=0.85)
    return fig

@_styled
def plot_sweep(sweeps, labels=None, channel: int = 0, scale: float = 1e-3,
               xlabel: str = "comprimento da fibra (km)", ber_floor: float = 1e-40):
    """Fator Q e BER ao longo de uma ou mais varreduras (pipelines.runner.sweep).

    Linha: BER estimada por ½erfc(Q/√2), mostrada a partir de `ber_floor` (com Q = 27 a formula
    da 1e-162, que nenhuma medida alcanca). Circulos vazios: BER contada, so nos pontos com erros
    (onde ha erros suficientes para contar, ela e a referencia da estimativa por Q)."""
    if not isinstance(sweeps, (list, tuple)):
        sweeps = [sweeps]
    labels = labels or [s.parameter for s in sweeps]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    handles = []
    any_counted = False
    for s, lab in zip(sweeps, labels):
        x = s.values * scale
        q = s.metric("q", channel)
        line, = ax[0].plot(x, 20 * np.log10(q), "o-", ms=3)
        handles.append((line, lab))
        ber = s.metric("ber_q", channel)
        shown = ber >= ber_floor
        ax[1].semilogy(x[shown], ber[shown], "o-", ms=3, color=line.get_color())
        err = s.metric("errors", channel) > 0
        if err.any():
            any_counted = True
            ax[1].semilogy(x[err], s.metric("ber_counted", channel)[err], "o", mfc="none",
                           ms=6, color=line.get_color())
    ax[0].set(xlabel=xlabel, ylabel=r"$20\,\log_{10}Q$ (dB)")
    ax[1].set(xlabel=xlabel, ylabel="BER", ylim=(ber_floor, 1.0))
    panel_label(ax[0], "(a) fator Q")
    panel_label(ax[1], r"(b) BER por $\frac{1}{2}\mathrm{erfc}(Q/\sqrt{2})$")
    hs, ls = [h for h, _ in handles], [lb for _, lb in handles]
    if any_counted:
        hs.append(Line2D([], [], marker="o", mfc="none", ls="", color="gray"))
        ls.append("BER contada")
    legend_above(fig, hs, ls, ncol=len(hs))
    return fig
