"""CLI minima.

    python -m simcopy.cli validate examples/400zr_4ch_75ghz.yaml
    python -m simcopy.cli run examples/imdd_ook_mzm.yaml -o resultados/mzm
    python -m simcopy.cli sweep examples/imdd_ook_dml.yaml --param link.spans[0].fiber.length \
        --range 20e3 40e3 1e3 -o resultados/dml

`run` e `sweep` aceitam YAML ou JSON. Com -o, gravam as metricas em JSON e, se o matplotlib
estiver instalado (pip install -e ".[figures]"), um grafico por estagio, a cadeia e o olho.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import numpy as np
from .config.io import load
from .config.resolve import resolve
from .validation.rules import validate
from .validation.base import Severity, ValidationFailed

def _describe(system, simulation, res):
    print(f"sistema : {system.n_channels} canais, espacamento "
          f"{system.channel_spacing/1e9:.0f} GHz, centro "
          f"{system.center_frequency/1e12:.2f} THz")
    print(f"solver  : {'vetorial (Manakov)' if system.link_is_vectorial else 'escalar'}")
    print(f"grade   : fs = {simulation.sampling_rate/1e9:.2f} GHz | N = {res.n_samples} "
          f"| janela = {res.window*1e6:.3f} us")
    for c, rc in zip(system.channels, res.channels):
        print(f"  canal {c.index}: f0 = {rc.f0/1e9:+7.1f} GHz | {c.modulation.name} "
              f"{c.pol.name} | SpS = {rc.samples_per_symbol} | "
              f"{rc.n_symbols} simbolos | {rc.geometry.n_bits} bits")
    for i, rf in enumerate(res.fibers):
        print(f"  span {i}: beta2 = {rf.beta2:.3e} s^2/m | gamma = {rf.gamma*1e3:.3f} "
              f"1/W/km | L_eff = {rf.effective_length/1e3:.1f} km")

def _save_json(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, default=float),
                    encoding="utf-8")
    print(f"gravado: {path}")

def _save_figures(out: Path, result):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from .plotting import plot_chain, plot_eye, plot_stage
    except ImportError:
        print("matplotlib ausente: graficos nao gravados (pip install -e \".[figures]\").")
        return
    out.mkdir(parents=True, exist_ok=True)
    figs = []
    for ch in result.metrics:
        figs += [(f"cadeia_ch{ch}", plot_chain(result, channel=ch)),
                 (f"olho_ch{ch}", plot_eye(result, ch))]
    figs += [(f"estagio_{k:02d}_{sid}", plot_stage(result, sid))
             for k, sid in enumerate(result.ids) if result.stage(sid).inputs]
    for name, fig in figs:
        fig.savefig(out / f"{name}.png", dpi=110)
        plt.close(fig)
    print(f"graficos: {out} ({len(figs)} arquivos)")

def main(argv=None):
    ap = argparse.ArgumentParser(prog="simcopy")
    ap.add_argument("command", choices=["validate", "resolve", "run", "sweep"])
    ap.add_argument("config", help="cenario em YAML ou JSON")
    ap.add_argument("-o", "--out", help="pasta para metricas (JSON) e graficos")
    ap.add_argument("--param", help="sweep: campo varrido, ex. link.spans[0].fiber.length")
    ap.add_argument("--range", nargs=3, type=float, metavar=("INICIO", "FIM", "PASSO"),
                    help="sweep: valores de INICIO a FIM (inclusive) com PASSO")
    a = ap.parse_args(argv)

    system, simulation = load(a.config)
    res = resolve(system, simulation)
    _describe(system, simulation, res)

    if a.command == "resolve":
        return 0
    if a.command == "validate":
        viol = validate(system, simulation)
        print()
        if not viol:
            print("validacao: nenhuma violacao.")
            return 0
        for v in viol:
            print(v); print()
        return 1 if any(v.severity is Severity.ERROR for v in viol) else 0

    from .pipelines.runner import run, sweep
    print()
    try:
        if a.command == "run":
            r = run(system, simulation)
            for v in r.violations:
                print(v); print()
            for s in r.stages:
                print(f"  {s.id:<24} {type(s.block).__name__:<22} {s.seconds:8.2f} s")
            print(r.summary())
            if a.out:
                out = Path(a.out)
                _save_json(out / "metricas.json", r.to_dict())
                _save_figures(out, r)
            return 0
        if not (a.param and a.range):
            ap.error("sweep precisa de --param e --range")
        start, stop, step = a.range
        values = np.arange(start, stop + 0.5 * step, step)

        def progresso(i, v, r):
            m = r.metrics[min(r.metrics)]
            print(f"  {a.param} = {v:.6g}: Q = {m.q:.3f} | BER por Q = {m.ber_q:.3e}")

        s = sweep(system, simulation, a.param, values, on_point=progresso)
        if a.out:
            out = Path(a.out)
            _save_json(out / "varredura.json", s.to_dict())
            try:
                import matplotlib
                matplotlib.use("Agg")
                from .plotting import plot_sweep
                fig = plot_sweep(s, labels=[Path(a.config).stem], scale=1.0, xlabel=a.param)
                fig.savefig(out / "varredura.png", dpi=110)
                print(f"grafico: {out / 'varredura.png'}")
            except ImportError:
                print("matplotlib ausente: grafico nao gravado.")
        return 0
    except ValidationFailed as e:
        print(e)
        return 1

if __name__ == "__main__":
    sys.exit(main())
