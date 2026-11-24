"""Figuras em PDF, uma por comprimento: BER contra potência, VPI contra três buscas de fase do simcopy.

As curvas são as da última figura da hipótese 4 do notebook 09: o VPI (conjuntos de configuração
conhecida: 2026 de 125 a 140 km e o de 120 km da pasta 400gbps_CH2) e o simcopy com a busca de fase
de hoje, com a da hipótese 1 (dois estágios por polarização, janela grossa de 32) e com o modelo da
hipótese 4 (fase conjunta, quarta potência desdobrada em 45°, janela grossa de 25). Os pontos do
simcopy têm 2^14 símbolos e o equalizador de hoje, calculado uma vez por sinal.

Uso, na pasta repo:
  python scripts/experimentos/figuras_ber_hipotese4.py [PASTA_SAIDA]
Grava ber_120km.pdf ... ber_140km.pdf (padrão: reference/vpi/hipoteses/figuras).
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from simcopy.plotting import apply_style, legend_above, panel_label, save_figure  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
CURVES = {  # chave: (rotulo, cor, marcador)
    "base": ("simcopy, busca de hoje", "gray", "s"),
    "por_polarizacao_particionada_22_grossa32": ("simcopy, hipótese 1", "blue", "^"),
    "conjunta45_grossa25_fina22.5": ("simcopy, hipótese 4 (modelo)", "red", "D"),
}


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "reference" / "vpi" / "hipoteses" / "figuras"
    apply_style(inline=None)
    catalog = json.loads((REPO / "reference/vpi/catalogo/catalogo_vpi.json").read_text(encoding="utf-8"))["points"]
    vpi = {}
    for p in catalog:
        if p["set"] == "2026_4ch" or (p["set"] == "400gbps_CH2" and p["length_km"] == 120):
            vpi.setdefault(int(p["length_km"]), []).append((p["power_dbm"], p["ber"]))
    sim = json.loads((REPO / "reference/vpi/hipoteses/h4_estados.json").read_text(encoding="utf-8"))["points"]
    for L in sorted(vpi):
        fig, ax = plt.subplots(figsize=(6.4, 4.6), layout="constrained")
        v = sorted(vpi[L])
        handles = [ax.semilogy([x for x, _ in v], [b for _, b in v], "-", color="black", marker="o", mfc="white",
                               ms=5.5, mew=1.2, lw=1.6, label="VPI", zorder=3)[0]]
        for key, (label, color, marker) in CURVES.items():
            pts = sorted((r["power_dbm"], r["ber"]) for r in sim if r["variant"] == key and int(r["length_km"]) == L)
            handles.append(ax.semilogy([x for x, _ in pts], [b for _, b in pts], "-", color=color, marker=marker,
                                       ms=4.5, lw=1.3, label=label)[0])
        ax.set(xlabel="potência total de lançamento (dBm)", ylabel="BER", ylim=(1e-4, 1), xlim=(-0.5, 20.5),
               xticks=range(0, 21, 2))
        panel_label(ax, f"{L} km", corner="lower right")
        legend_above(ax, handles, [h.get_label() for h in handles], ncol=2)
        print(save_figure(fig, out / f"ber_{L}km"))
        plt.close(fig)


if __name__ == "__main__":
    main()
