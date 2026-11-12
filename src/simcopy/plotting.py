"""Graficos de interesse (docs/mapa_arquitetura.md, 5.12; docs/sistema_imdd.md, secao 5).

Uma funcao serve a qualquer componente porque os paineis sao escolhidos pelo dominio do sinal,
nao pelo componente:
- sinal eletrico: forma de onda no tempo e espectro;
- sinal optico: potencia e chirp no tempo e espectro optico;
- bits: degraus no tempo.
Os blocos nao plotam (CLAUDE.md, regra 6): estes graficos leem as saidas que o executor guarda
(Result.output) e as contas de simcopy.analysis. Precisa do matplotlib (pip install -e ".[figures]").

    from simcopy.plotting import plot_stage, plot_chain, plot_eye, plot_sweep
    plot_stage(r, "link.span0.fiber")     # entrada e saida da fibra, tempo e frequencia
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
from . import analysis
from .core.signals import Bits, Domain, Waveform

IN_COLOR, OUT_COLOR = "tab:blue", "tab:orange"

def _bits_trace(b: Bits) -> tuple[np.ndarray, np.ndarray]:
    sps = b.geometry.samples_per_symbol
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

def plot_io(entrada, saida, labels=("entrada", "saída"), title: str = "",
            units=(None, None), time_window: float | None = None,
            freq_span: float | None = None, resolution: float = 0.5e9):
    """Entrada e saida de um componente, no tempo e na frequencia, com os paineis escolhidos pelo
    dominio de cada sinal (Waveform eletrico ou optico, Bits ou None).

    time_window: duracao mostrada no tempo [s] (padrao: 32 simbolos, ou a janela se for menor).
    freq_span: faixa de frequencia mostrada, +-freq_span no optico e 0..freq_span no eletrico [Hz].
    resolution: media movel dos espectros [Hz] (analysis.optical_spectrum).
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
            a.set(xlabel="t [ns]", ylabel="bit", title="tempo: bits", xlim=(0, t_max))
            a.legend(loc="upper right", fontsize=8)
        elif kind == "electrical":
            a = axes.pop(0)
            for w, lab, cor, un in el:
                a.plot(t_ns, w.samples[0].real, color=cor, label=f"{lab} [{un or 'u.a.'}]")
            a.set(xlabel="t [ns]", ylabel="sinal elétrico", title="tempo: sinal elétrico",
                  xlim=(0, t_max))
            a.legend(loc="upper right", fontsize=8)
        else:
            a = axes.pop(0)
            for w, lab, cor, _ in op:
                a.plot(t_ns, analysis.optical_power(w) * 1e3, color=cor, label=lab)
            a.set(xlabel="t [ns]", ylabel="P [mW]", title="tempo: potência óptica",
                  xlim=(0, t_max))
            a.legend(loc="upper right", fontsize=8)
            a = axes.pop(0)
            for j, (w, lab, cor, _) in enumerate(op):
                a.plot(t_ns, analysis.chirp(w) / 1e9, color=cor, label=lab,
                       ls="-" if j == 0 else "--")
            a.set(xlabel="t [ns]", ylabel="δν [GHz]", title="tempo: chirp (−dφ/dt / 2π)",
                  xlim=(0, t_max))
            a.legend(loc="upper right", fontsize=8)
    for kind in order:
        if kind == "electrical":
            a = axes.pop(0)
            for w, lab, cor, _ in el:
                f, s = analysis.electrical_spectrum(w, resolution)
                a.plot(f / 1e9, 10 * np.log10(s / s.max() + 1e-15), color=cor, label=lab)
            span = (freq_span or f.max()) / 1e9
            a.set(xlabel="f [GHz]", ylabel="dB (norm.)", ylim=(-70, 3), xlim=(0, span),
                  title="frequência: espectro elétrico")
            a.legend(loc="upper right", fontsize=8)
        elif kind == "optical":
            a = axes.pop(0)
            top = max(analysis.optical_spectrum(w, resolution)[1].max() for w, *_ in op)
            for j, (w, lab, cor, _) in enumerate(op):
                f, s = analysis.optical_spectrum(w, resolution)
                a.plot(f / 1e9, 10 * np.log10(s / top + 1e-15), color=cor, label=lab,
                       ls="-" if j == 0 else "--")
            span = (freq_span or abs(f).max()) / 1e9
            a.set(xlabel="frequência óptica relativa ao centro da banda [GHz] (positivo = azul)",
                  ylabel="dB (rel. ao maior)", ylim=(-70, 3), xlim=(-span, span),
                  title="frequência: espectro óptico")
            a.legend(loc="upper right", fontsize=8)
    if title:
        fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    return fig

def _unit(block):
    return getattr(block, "output_unit", None)

def plot_stage(result, stage_id: str, **kw):
    """plot_io da entrada e da saida de um estagio do resultado (Result.ids lista os estagios).
    No mux, a entrada mostrada e a do primeiro canal."""
    st = result.stage(stage_id)
    prev = result.stage(st.inputs[0]) if st.inputs else None
    entrada = prev.output if prev is not None else None
    labels = (f"entrada ({prev.id})" if prev else "entrada", f"saída ({st.id})")
    units = (_unit(prev.block) if prev else None, _unit(st.block))
    return plot_io(entrada, st.output, labels, title=f"{st.id}: {type(st.block).__name__}",
                   units=units, **kw)

def plot_chain(result, time_window: float | None = 2e-9, channel: int = 0):
    """Visao da cadeia no tempo: um painel por estagio no caminho do canal (bits, sinal eletrico,
    potencia optica, corrente), como a celula do MZM no notebook 07."""
    stages = result.path(channel)
    fig, axes = plt.subplots(len(stages), 1, figsize=(9, 1.7 * len(stages)), sharex=True)
    axes = np.atleast_1d(axes)
    for a, st in zip(axes, stages):
        s = st.output
        if isinstance(s, Bits):
            ref = next(x.output for x in stages if isinstance(x.output, Waveform))
            idx, x = _bits_trace(s)
            a.step(ref.grid.time[idx % ref.grid.n_samples] * 1e9, x, where="post", color="k")
            a.set(ylabel="bit")
        elif s.domain is Domain.OPTICAL:
            a.plot(s.grid.time * 1e9, analysis.optical_power(s) * 1e3, color=OUT_COLOR)
            a.set(ylabel="P [mW]")
        else:
            un = _unit(st.block)
            y, lab = s.samples[0].real, un or "u.a."
            if un == "A":
                y, lab = y * 1e3, "mA"
            a.plot(s.grid.time * 1e9, y, color=IN_COLOR)
            a.set(ylabel=lab)
        a.set_title(f"{st.id} ({type(st.block).__name__})", fontsize=9)
    axes[-1].set(xlabel="t [ns]")
    if time_window is not None:
        axes[-1].set_xlim(0, time_window * 1e9)
    fig.tight_layout()
    return fig

def plot_eye(result, channel: int = 0, max_traces: int = 400):
    """Olho da corrente no fotodetector, com o instante de decisao, os niveis I1 e I0, o limiar
    otimo e o fator Q medidos (kernels/decision.py)."""
    m = result.metrics[channel]
    i = result.output(f"ch{channel}.rx.photodiode").samples[0].real
    sps = result.resolved.channels[channel].samples_per_symbol
    last = min(m.last_symbol, m.first_symbol + max_traces)
    t, tr = analysis.eye(i, sps, m.delay, m.first_symbol, last)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw={"width_ratios": [3, 1]})
    ax[0].plot(t, tr.T * 1e3, color=IN_COLOR, lw=0.4, alpha=0.35)
    for y, lab, ls in ((m.level_one, "I₁", "-"), (m.level_zero, "I₀", "-"),
                       (m.threshold, "limiar I_D", "--")):
        ax[0].axhline(y * 1e3, color="k", ls=ls, lw=0.8)
        ax[0].text(t[-1], y * 1e3, f" {lab}", va="center", fontsize=8)
    ax[0].axvline(0, color="k", lw=0.6, ls=":")
    ax[0].set(xlabel="tempo [símbolos] (0 = instante de decisão)", ylabel="corrente [mA]",
              title=f"Diagrama de olho, canal {channel}: Q = {m.q:.2f}, BER por Q = {m.ber_q:.2e}")
    bits = result.output(f"ch{channel}.tx.bits").data[0].astype(bool)
    k = np.arange(m.first_symbol, m.last_symbol)
    s = i[(k * sps + m.delay) % i.size]
    ax[1].hist(s[bits[k]] * 1e3, bins=40, orientation="horizontal", color=OUT_COLOR, alpha=0.7,
               label="bits 1")
    ax[1].hist(s[~bits[k]] * 1e3, bins=40, orientation="horizontal", color=IN_COLOR, alpha=0.7,
               label="bits 0")
    ax[1].axhline(m.threshold * 1e3, color="k", ls="--", lw=0.8)
    ax[1].set(xlabel="contagem", title="no instante de decisão")
    ax[1].sharey(ax[0])
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    return fig

def plot_sweep(sweeps, labels=None, channel: int = 0, scale: float = 1e-3,
               xlabel: str = "comprimento da fibra [km]", ber_floor: float = 1e-40):
    """Fator Q e BER ao longo de uma ou mais varreduras (pipelines.runner.sweep).

    Linha: BER estimada por ½erfc(Q/√2), mostrada a partir de `ber_floor` (com Q = 27 a formula
    da 1e-162, que nenhuma medida alcanca). Circulos vazios: BER contada, so nos pontos com erros
    (onde ha erros suficientes para contar, ela e a referencia da estimativa por Q)."""
    if not isinstance(sweeps, (list, tuple)):
        sweeps = [sweeps]
    labels = labels or [s.parameter for s in sweeps]
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    for s, lab in zip(sweeps, labels):
        x = s.values * scale
        q = s.metric("q", channel)
        line, = ax[0].plot(x, 20 * np.log10(q), "o-", ms=3, label=lab)
        ber = s.metric("ber_q", channel)
        shown = ber >= ber_floor
        ax[1].semilogy(x[shown], ber[shown], "o-", ms=3, color=line.get_color(),
                       label=f"{lab}: ½erfc(Q/√2)")
        err = s.metric("errors", channel) > 0
        if err.any():
            ax[1].semilogy(x[err], s.metric("ber_counted", channel)[err], "o", mfc="none",
                           ms=6, color=line.get_color(), label=f"{lab}: contada")
    ax[0].set(xlabel=xlabel, ylabel="20·log10(Q) [dB]", title="fator Q")
    ax[1].set(xlabel=xlabel, ylabel="BER", ylim=(ber_floor, 1.0),
              title=f"BER por Q e contada (abaixo de {ber_floor:.0e}: fora do gráfico)")
    for a in ax:
        a.grid(True, which="both", alpha=0.3)
        a.legend(fontsize=8)
    fig.tight_layout()
    return fig
