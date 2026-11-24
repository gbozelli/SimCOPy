"""BER contra potencia de lancamento, VPI contra simcopy, de 80 a 140 km.

A BER do VPI e contada aqui, nas constelacoes exportadas pelo grupo, com o mesmo mapeamento de Gray
do simcopy: as tabelas de BER do grupo tem a duvida da contagem (pergunta 13 do documento ao grupo)
e o VPI tem `GrayMapping = No` (pergunta 14), entao nao se comparam direto. Em cada polarizacao: os
simbolos enviados saem da constelacao de referencia (SNR de 22 dB, a decisao nela nao erra), o ganho
e tirado por minimos quadrados (como na medida do simcopy), cada simbolo recebido e decidido no
16QAM e os bits das duas decisoes sao comparados.

Conjuntos do VPI (pasta `extraido` dos dados do grupo, fora do repositorio):
  80 km       400gbps_CH2/Dados_80km: OUTRO cenario (sem os 5 graus, outros bits, configuracao
              desconhecida: pergunta 10 do documento ao grupo). Nao se compara com o simcopy.
  120 km      400gbps_CH2/Dados_120km (com os 5 graus)
  125-140 km  2026_SIMULATION/dados_4_canais

Uso, na pasta repo:
  python scripts/ber_vpi_simcopy.py vpi PASTA_EXTRAIDO -o reference/vpi/comparacao/ber_vpi.json
  python scripts/ber_vpi_simcopy.py figura reference/vpi/comparacao/ber_vpi.json SIMCOPY.json [...]
      -o reference/vpi/comparacao/ber_potencia.png [--dados-simcopy reference/vpi/comparacao/ber_simcopy.json]
Os SIMCOPY.json sao saidas de `scripts/compara_vpi.py simular`. Um ponto repetido fica com o do
ultimo arquivo.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from simcopy.kernels.mapping import constellation, demap  # noqa: E402

LEVELS = np.array([-3.0, -1.0, 1.0, 3.0])
TO_SIMCOPY = np.abs(constellation("16QAM").real).max() / 3.0     # grade +-1, +-3 -> a do simcopy
DIFFERENT_SCENARIO = (80,)


def decide(x: np.ndarray) -> np.ndarray:
    q = lambda v: LEVELS[np.clip(np.round((v + 3.0) / 2.0), 0, 3).astype(int)]
    return q(x.real) + 1j * q(x.imag)


def count(x: np.ndarray, ref: np.ndarray) -> tuple[int, int, int]:
    """(bits errados, simbolos errados, simbolos) de uma polarizacao."""
    a = decide(ref)
    y = x / (np.sum(x * np.conj(a)) / np.sum(np.abs(a) ** 2))
    d = decide(y)
    bits = np.count_nonzero(demap(d * TO_SIMCOPY, "16QAM") != demap(a * TO_SIMCOPY, "16QAM"))
    return int(bits), int(np.count_nonzero(d != a)), int(a.size)


def load_columns(path: Path) -> list[np.ndarray]:
    """Formato 400gbps_CH2: cabecalho X_real,X_imag,Y_real,Y_imag."""
    v = np.loadtxt(path, delimiter=",", skiprows=1)
    return [v[:, 0] + 1j * v[:, 1], v[:, 2] + 1j * v[:, 3]]


def load_single(path: Path) -> np.ndarray:
    return np.loadtxt(path, delimiter=",", comments=("//", "#")) @ np.array([1, 1j])


def vpi_points(folder: Path) -> list[dict]:
    points = []

    def add(length, p, pols, refs):
        e = [count(x, r) for x, r in zip(pols, refs)]
        bits, symbols, n = (sum(v[i] for v in e) for i in range(3))
        points.append({"length_km": length, "power_dbm": p, "ber": bits / (4 * n), "bit_errors": bits,
                       "n_bits": 4 * n, "ser": symbols / n})
        print(f"VPI {length} km {p:2d} dBm: BER = {bits / (4 * n):.3e} ({bits} bits)", flush=True)

    for length in (80, 120):
        base = folder / "400gbps_CH2" / f"Dados_{length}km"
        refs = load_columns(base / f"400gbps_IdealConstellationDiagram_ref_{length}km.csv")
        for p in range(0, 21):
            f = base / f"400gbps_RealConstellationDiagram_{p}dbm_{length}km.csv"
            if f.exists():
                add(length, p, load_columns(f), refs)

    base = folder / "2026_SIMULATION" / "dados_4_canais"
    refs = [load_single(base / "dados_ideais" / "dados_ref_canal2" / f"dados_4_canais_ref_{pol}_pol_canal2.csv")
            for pol in "XY"]
    sub = {125: "dados_125km", 130: "dados_130km", 135: "dados_135km", 140: "dados_140m"}
    for length, folder_name in sub.items():
        d = base / "dados_distorcidos" / folder_name / f"{folder_name.replace('140m', '140km')}_canal_2"
        for p in range(0, 21):
            files = [d / f"dados_4_canais_canal_2_{length}km_{pol}_pol_{p}dBm.csv" for pol in "XY"]
            if all(f.exists() for f in files):
                add(length, p, [load_single(f) for f in files], refs)
    return points


def simcopy_points(files: list[Path]) -> list[dict]:
    merged = {}
    for f in files:
        for q in json.loads(f.read_text(encoding="utf-8"))["points"]:
            n_bits = q.get("n_bits")
            merged[(int(q["length_km"]), int(q["power_dbm"]))] = {
                "length_km": int(q["length_km"]), "power_dbm": int(q["power_dbm"]), "ber": q["ber"],
                "bit_errors": q.get("bit_errors"), "n_bits": n_bits, "snr_db": q.get("snr_db"),
                "phase_slips": q.get("phase_slips")}
    return [merged[k] for k in sorted(merged)]


def figure(vpi: list[dict], sim: list[dict], out: Path) -> None:
    """Um painel por comprimento, no estilo das figuras do simcopy (src/simcopy/plotting.py, guia
    de figuras): VPI em circulos vazados, simcopy em linha azul; ponto do simcopy sem erros no piso
    de meia contagem, como triangulo."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from simcopy import plotting as pl

    vpi_color, sim_color = pl.COLORS["VPI"], pl.COLORS["simcopy"]
    lengths = sorted({p["length_km"] for p in vpi} | {p["length_km"] for p in sim})
    with matplotlib.rc_context(pl.FIGURE_STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(10, 6.2), sharex=True, sharey=True)
        for k, (ax, length) in enumerate(zip(axes.ravel(), lengths)):
            v = sorted((p for p in vpi if p["length_km"] == length), key=lambda p: p["power_dbm"])
            s = sorted((p for p in sim if p["length_km"] == length), key=lambda p: p["power_dbm"])
            floor = lambda p: p["ber"] if p["ber"] > 0 else 0.5 / p["n_bits"]
            if s:
                ax.plot([p["power_dbm"] for p in s], [floor(p) for p in s], "-", color=sim_color, zorder=2)
                ax.plot([p["power_dbm"] for p in s if p["ber"] > 0], [p["ber"] for p in s if p["ber"] > 0],
                        "o", color=sim_color, ms=3.5, zorder=3)
                zero = [p for p in s if p["ber"] == 0]
                if zero:
                    ax.plot([p["power_dbm"] for p in zero], [floor(p) for p in zero], "v", color=sim_color,
                            ms=5, zorder=3)
            ax.plot([p["power_dbm"] for p in v], [p["ber"] for p in v], "o", mfc="none", color=vpi_color,
                    mew=1.1, ms=6, zorder=4)
            ax.set_yscale("log")
            ax.set_ylim(1e-6, 1)        # a decada de baixo fica livre para o rotulo do painel
            ax.set_xlim(-0.5, 20.5)
            ax.set_xticks(range(0, 21, 4))
            note = ", outro cenário no VPI" if length in DIFFERENT_SCENARIO else ""
            pl.panel_label(ax, f"({'abcdef'[k]}) {length} km{note}", corner="lower right")
            pl.decimal_comma(ax, "x")
        for ax in axes[-1]:
            ax.set_xlabel("potência total (dBm)")
        for ax in axes[:, 0]:
            ax.set_ylabel("BER")
        handles = [Line2D([], [], color=vpi_color, marker="o", mfc="none", mew=1.1, ms=6, lw=0),
                   Line2D([], [], color=sim_color, marker="o", ms=3.5)]
        labels = ["VPI", "simcopy"]
        if any(p["ber"] == 0 for p in sim):
            handles.append(Line2D([], [], color=sim_color, marker="v", ms=5, lw=0))
            labels.append("simcopy sem erros (no piso de meia contagem)")
        pl.legend_above(fig, handles, labels, ncol=len(handles))
        pdf = pl.save_figure(fig, out, png=True)
    plt.close(fig)
    print(f"figura: {pdf} (e o PNG ao lado)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("vpi")
    v.add_argument("folder", type=Path)
    v.add_argument("-o", "--out", type=Path, required=True)
    f = sub.add_parser("figura")
    f.add_argument("vpi", type=Path)
    f.add_argument("simcopy", type=Path, nargs="+")
    f.add_argument("-o", "--out", type=Path, required=True)
    f.add_argument("--dados-simcopy", type=Path, help="grava os pontos do simcopy usados na figura")
    args = ap.parse_args(argv)
    sys.stdout.reconfigure(encoding="utf-8")
    if args.cmd == "vpi":
        points = vpi_points(args.folder)
        args.out.write_text(json.dumps({"source": "vpi", "ber": "Gray, contada nas constelacoes",
                                        "points": points}, indent=1), encoding="utf-8")
    else:
        vpi = json.loads(args.vpi.read_text(encoding="utf-8"))["points"]
        sim = simcopy_points(args.simcopy)
        if args.dados_simcopy:
            args.dados_simcopy.write_text(json.dumps({"source": "simcopy", "points": sim}, indent=1),
                                          encoding="utf-8")
        figure(vpi, sim, args.out)


if __name__ == "__main__":
    main()
