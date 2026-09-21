"""CLI minima: `python -m simcopy.cli validate examples/400zr_4ch_75ghz.yaml`"""
from __future__ import annotations
import argparse, sys
from .config.io import load
from .config.resolve import resolve
from .validation.rules import validate
from .validation.base import Severity

def main(argv=None):
    ap = argparse.ArgumentParser(prog="simcopy")
    ap.add_argument("command", choices=["validate", "resolve"])
    ap.add_argument("config")
    a = ap.parse_args(argv)

    system, simulation = load(a.config)
    res = resolve(system, simulation)

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

    if a.command == "resolve":
        return 0
    viol = validate(system, simulation)
    print()
    if not viol:
        print("validacao: nenhuma violacao.")
        return 0
    for v in viol:
        print(v); print()
    return 1 if any(v.severity is Severity.ERROR for v in viol) else 0

if __name__ == "__main__":
    sys.exit(main())
