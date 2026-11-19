"""Composicao do ruido na saida do DSP: tangencial (fase) contra radial (amplitude), VPI x simcopy.

Para cada ponto em que o VPI converge (2026, 125 a 140 km), alinha a saida com os simbolos
conhecidos, tira um ganho complexo por polarizacao e decompoe o erro e = y/g - a na direcao de a
(radial) e perpendicular a ela (tangencial). Com ruido branco circular a razao
var(tangencial)/var(radial) e 1. Ruido de fase que sobra depois da busca de fase a deixa acima de 1,
mais no anel de fora. O simcopy usa os sinais na saida do ADC guardados por
`scripts/compara_vpi.py simular --save-signals` e o receptor cego com janela de 43.

Uso, na pasta repo (os dados do VPI ficam fora do repo):
  python scripts/experimentos/ruido_radial_tangencial.py [PASTA_VPI]
Grava reference/vpi/hipoteses/ruido_radial_tangencial.json.
"""
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
import catalogo_vpi as cat  # noqa: E402
from simcopy.kernels import dsp  # noqa: E402
from simcopy.kernels.demodulation import measure_blind  # noqa: E402
from simcopy.kernels.mapping import map_bits  # noqa: E402

DEFAULT_VPI = Path(r"C:\Users\bozel\Área de Trabalho\Pesquisa\25-26\dados4canais\extraido")
SIG = REPO / "resultados" / "replica_vpi" / "skew1amostra" / "sinais"
OUT = REPO / "reference" / "vpi" / "hipoteses" / "ruido_radial_tangencial.json"


def split(y, a):
    """Razao var(tangencial)/var(radial), no total e por anel, e a SNR."""
    g = np.sum(y * np.conj(a)) / np.sum(np.abs(a) ** 2)
    e = y / g - a
    u = a / np.abs(a)
    rad, tan = np.real(e * np.conj(u)), np.imag(e * np.conj(u))
    r2 = np.round(np.abs(a) ** 2 / np.min(np.abs(a) ** 2), 1)        # aneis: 1, 5, 9
    out = {"ratio": float(np.var(tan) / np.var(rad)),
           "snr_db": float(10 * np.log10(np.mean(np.abs(a) ** 2) / np.mean(np.abs(e) ** 2)))}
    for name, ring in (("inner", 1.0), ("middle", 5.0), ("outer", 9.0)):
        m = r2 == ring
        out[name] = float(np.var(tan[m]) / np.var(rad[m]))
    return out


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    vpi_folder = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_VPI
    catalog = json.loads((REPO / "reference/vpi/catalogo/catalogo_vpi.json").read_text(encoding="utf-8"))["points"]
    good = {(p["length_km"], p["power_dbm"]) for p in catalog
            if p["set"] == "2026_4ch" and p["rotated_fraction"] <= 0.02 and p["ber"] <= 0.1}
    results = []
    for set_name, L, p, files, refs in cat.datasets(vpi_folder):
        if set_name != "2026_4ch" or (L, p) not in good:
            continue
        sig = SIG / f"{L}km_{p}dBm_semente0.npz"
        if not sig.exists():
            continue
        ref = [cat.to_grid(r) for f in refs for r in cat.as_pols(f)]
        vpi = []
        for x in [x for f in files for x in cat.as_pols(f)]:
            best = None
            for r_index, a in enumerate(ref):
                peak, lag, conj = cat.best_alignment(cat.decide(x / np.sqrt(np.mean(np.abs(x) ** 2) / 10.0)), a)
                if best is None or peak > best[0]:
                    best = (peak, lag, conj, r_index)
            _, lag, conj, r_index = best
            n = min(x.size, ref[r_index].size)
            vpi.append(split((np.conj(x) if conj else x)[:n], np.roll(ref[r_index][:n], lag)))
        d = np.load(sig)
        z, _ = dsp.blind_receiver(d["adc"].astype(complex), "16QAM", window=43)
        _, y = measure_blind(z, d["bits"], "16QAM", int(d["guard"]))
        g = int(d["guard"])
        a = np.array([map_bits(row, "16QAM") for row in d["bits"]])
        a = a[:, g:a.shape[1] - g]
        sim = [split(y[i], a[i]) for i in range(y.shape[0])]
        rec = {"length_km": L, "power_dbm": p}
        for who, parts in (("vpi", vpi), ("simcopy", sim)):
            rec[who] = {k: float(np.mean([q[k] for q in parts])) for k in parts[0]}
        results.append(rec)
        v, s = rec["vpi"], rec["simcopy"]
        print(f"{L} km {p:2d} dBm | SNR {v['snr_db']:5.2f} x {s['snr_db']:5.2f} | tangencial/radial "
              f"{v['ratio']:.3f} x {s['ratio']:.3f} | anel de fora {v['outer']:.3f} x {s['outer']:.3f}", flush=True)
    OUT.write_text(json.dumps({"source": "VPI (2026_4ch) x simcopy (2^14 simbolos, semente 0, janela 43)",
                               "points": results}), encoding="utf-8")


if __name__ == "__main__":
    main()
