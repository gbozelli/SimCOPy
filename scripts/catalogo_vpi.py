"""Catalogo de todos os conjuntos de constelacoes do VPI do grupo (pasta `extraido`, fora do
repositorio): para cada ponto, BER e SER contadas, SNR e a fracao do bloco em que a fase se perdeu.

Serve para comparar o simcopy com o VPI em todos os conjuntos, e para ver onde o DSP do VPI perde a
fase. A BER usa o mesmo codigo de Gray do simcopy (scripts/ber_vpi_simcopy.py). Em cada
polarizacao: a referencia e decidida no 16QAM, o atraso e a orientacao (giro de 90 graus,
conjugado, troca de X e Y) saem do pico da correlacao com a referencia, o ganho complexo e tirado
por minimos quadrados, e cada simbolo e decidido.

Conjuntos (os parametros de simulacao so sao conhecidos nos de 2026 e no 400gbps_CH2 de 120 e 140 km):
  2026_4ch          2026_SIMULATION/dados_4_canais: 125 a 140 km, com os 5 graus (o cenario do repositorio)
  2026_4ch_sem5     2026_SIMULATION/dados_4_canais_sem_5_graus: 120 km
  400gbps_CH2       400gbps_CH2: 80 (outro cenario), 120 e 140 km (140 km igual ao de 2026)
  2025_4ch          daods_2025: 80 a 120 km (configuracao desconhecida)
  2025_1ch          daods_2025/dados_120km_canal1: um canal, 120 km
  2021_1ch          dados_1_canal_175km: um canal, 175 km e um comprimento nao informado
  2021_14GBd        dados_mistos_1_canal: 150 km a 14 GBd

Uso, na pasta repo:
  python scripts/catalogo_vpi.py PASTA_EXTRAIDO -o reference/vpi/catalogo/catalogo_vpi.json
"""
import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from simcopy.kernels.mapping import constellation, demap  # noqa: E402

LEVELS = np.array([-3.0, -1.0, 1.0, 3.0])
TO_SIMCOPY = np.abs(constellation("16QAM").real).max() / 3.0


def read_table(path: Path) -> np.ndarray:
    """Le um CSV numerico do VPI (cabecalhos //, #, linhas de texto) como matriz."""
    text = path.read_text(encoding="utf-8", errors="replace").splitlines()
    k = 0
    while k < len(text) and not re.match(r"^\s*[-+0-9.]", text[k]):
        k += 1
    body = [line for line in text[k:] if line.strip()]
    ncol = body[0].count(",") + 1
    values = np.fromstring(",".join(body), sep=",")
    return values.reshape(-1, ncol)


def as_pols(path: Path) -> list[np.ndarray]:
    """Uma polarizacao (2 colunas) ou duas (4 colunas: X real, X imag, Y real, Y imag)."""
    m = read_table(path)
    if m.shape[1] >= 4:
        return [m[:, 0] + 1j * m[:, 1], m[:, 2] + 1j * m[:, 3]]
    return [m[:, 0] + 1j * m[:, 1]]


def decide(x: np.ndarray) -> np.ndarray:
    q = lambda v: LEVELS[np.clip(np.round((v + 3.0) / 2.0), 0, 3).astype(int)]
    return q(x.real) + 1j * q(x.imag)


def to_grid(ref: np.ndarray) -> np.ndarray:
    """Referencia decidida na grade +-1, +-3."""
    return decide(ref * np.sqrt(10.0 / np.mean(np.abs(ref) ** 2)))


def best_alignment(x: np.ndarray, a: np.ndarray) -> tuple[float, int, bool]:
    """(|pico| normalizado, atraso, conjugado) da correlacao circular de x (ou conj x) com a."""
    n = min(x.size, a.size)
    x, a = x[:n], a[:n]
    best = (0.0, 0, False)
    fa = np.fft.fft(a)
    for conj in (False, True):
        xx = np.conj(x) if conj else x
        c = np.fft.ifft(np.fft.fft(xx) * np.conj(fa))
        k = int(np.argmax(np.abs(c)))
        peak = float(np.abs(c[k]) / np.sqrt(np.sum(np.abs(xx) ** 2) * np.sum(np.abs(a) ** 2)))
        if peak > best[0]:
            best = (peak, k, conj)
    return best


def measure(x: np.ndarray, a: np.ndarray, block: int = 128) -> dict:
    """Erros, SNR e fase de uma polarizacao ja alinhada."""
    g = np.sum(x * np.conj(a)) / np.sum(np.abs(a) ** 2)
    y = x / g
    d = decide(y)
    bits = int(np.count_nonzero(demap(d * TO_SIMCOPY, "16QAM") != demap(a * TO_SIMCOPY, "16QAM")))
    nb = (a.size // block) * block
    c = np.sum((y[:nb] * np.conj(a[:nb])).reshape(-1, block), axis=1)
    ang = np.degrees(np.angle(c ** 4)) / 4.0
    return {"bit_errors": bits, "symbol_errors": int(np.count_nonzero(d != a)), "n": int(a.size),
            "snr_db": float(10 * np.log10(np.mean(np.abs(a) ** 2) / np.mean(np.abs(y - a) ** 2))),
            "rotated_fraction": float(np.mean(np.abs(ang) > 22.5))}


def point(pols: list[np.ndarray], refs: list[np.ndarray]) -> dict:
    """Mede um ponto, escolhendo para cada polarizacao recebida a referencia, o atraso e a
    orientacao de maior correlacao."""
    out = []
    for x in pols:
        best = None
        for r_index, a in enumerate(refs):
            peak, lag, conj = best_alignment(decide(x / np.sqrt(np.mean(np.abs(x) ** 2) / 10.0)), a)
            if best is None or peak > best[0]:
                best = (peak, lag, conj, r_index)
        peak, lag, conj, r_index = best
        n = min(x.size, refs[r_index].size)
        xx = (np.conj(x) if conj else x)[:n]
        a = np.roll(refs[r_index][:n], lag)
        # giro de 90 graus: o ganho complexo do ajuste ja o absorve
        m = measure(xx, a)
        out.append({**m, "lag": lag, "conjugated": conj, "reference": r_index,
                    "correlation": peak})
    bits = sum(o["bit_errors"] for o in out)
    syms = sum(o["symbol_errors"] for o in out)
    n = sum(o["n"] for o in out)
    snr = -10 * np.log10(np.mean([10 ** (-o["snr_db"] / 10) for o in out]))
    return {"ber": bits / (4 * n), "ser": syms / n, "snr_db": float(snr), "bit_errors": bits, "n_bits": 4 * n,
            "rotated_fraction": float(np.mean([o["rotated_fraction"] for o in out])), "per_pol": out}


def power_of(name: str) -> int | None:
    m = re.search(r"(-?\d+)dbm", name, re.I)
    return int(m.group(1)) if m else None


def datasets(root: Path):
    """(conjunto, comprimento [km], potencia [dBm], [arquivos], [referencias])"""
    ref26 = root / "2026_SIMULATION/dados_4_canais/dados_ideais/dados_ref_canal2"
    r26 = [ref26 / f"dados_4_canais_ref_{p}_pol_canal2.csv" for p in "XY"]
    sub = {125: "dados_125km/dados_125km_canal_2", 130: "dados_130km/dados_130km_canal_2",
           135: "dados_135km/dados_135km_canal_2", 140: "dados_140m/dados_140km_canal_2"}
    for L, s in sub.items():
        d = root / "2026_SIMULATION/dados_4_canais/dados_distorcidos" / s
        for p in range(0, 21):
            f = [d / f"dados_4_canais_canal_2_{L}km_{q}_pol_{p}dBm.csv" for q in "XY"]
            if all(x.exists() for x in f):
                yield "2026_4ch", L, p, f, r26
    d = root / "2026_SIMULATION/dados_4_canais_sem_5_graus/dados_distorcidos/dados_120km/dados_120km_canal2"
    for p in range(0, 21):
        f = [d / f"dados_4_canais_sem_5graus_canal_2_120km_{q}_pol_{p}dBm.csv" for q in "XY"]
        if all(x.exists() for x in f):
            yield "2026_4ch_sem5", 120, p, f, r26
    for L in (80, 120, 140):
        d = root / "400gbps_CH2" / f"Dados_{L}km"
        ref = [d / f"400gbps_IdealConstellationDiagram_ref_{L}km.csv"]
        for p in range(0, 21):
            f = d / f"400gbps_RealConstellationDiagram_{p}dbm_{L}km.csv"
            if f.exists():
                yield "400gbps_CH2", L, p, [f], ref
    r25 = [root / "daods_2025" / f"dados_{q}_ref.csv" for q in "XY"]
    for L in (80, 90, 100, 110, 120):
        for p in range(-5, 21):
            d = root / "daods_2025" / f"dados_{L}km" / f"dados_{L}km_{p}dBm"
            f = [d / f"dados_{q}_{L}km_{p}dBm.csv" for q in "XY"]
            if all(x.exists() for x in f):
                yield "2025_4ch", L, p, f, r25
    d = root / "daods_2025/dados_120km_canal1"
    r1 = [root / "daods_2025" / f"dados_1_canal_ref_{q}_pol_canal1.csv" for q in "XY"]
    for p in range(-5, 21):
        f = [d / f"dados_1_canal_canal_1_120km_{q}_pol_{p}dBm_por_canal.csv" for q in "XY"]
        if all(x.exists() for x in f):
            yield "2025_1ch", 120, p, f, r1
    d = root / "dados_1_canal_175km"
    for f in sorted(d.glob("DP_RealConstellationDiagram_*dbm_175km.csv")):
        yield "2021_1ch", 175, power_of(f.name), [f], [d / "DP_IdealConstellationDiagram_0dbm_175km.csv"]
    for f in sorted(d.glob("DP_RealConstellationDiagram_*dbm.csv")):
        yield "2021_1ch", None, power_of(f.name), [f], [d / "DP_IdealConstellationDiagram_0dbm.csv"]
    d = root / "dados_mistos_1_canal"
    for f in sorted(d.glob("DP_RealConstellationDiagram_*dbm_150km_14e9Bd.csv")):
        yield "2021_14GBd", 150, power_of(f.name), [f], [d / "DP_IdealConstellationDiagram_0dbm.csv",
                                                         d / "DP_IdealConstellationDiagram_0dbm_175km.csv"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("folder", type=Path)
    ap.add_argument("-o", "--out", type=Path, required=True)
    args = ap.parse_args()
    sys.stdout.reconfigure(encoding="utf-8")
    cache: dict[tuple, list] = {}
    points = json.loads(args.out.read_text(encoding="utf-8"))["points"] if args.out.exists() else []
    done = {(q["set"], q["length_km"], q["power_dbm"]) for q in points}
    for name, L, p, files, refs in datasets(args.folder):
        if (name, L, p) in done:
            continue
        t0 = time.perf_counter()
        key = tuple(str(r) for r in refs)
        if key not in cache:
            cache[key] = [to_grid(x) for r in refs for x in as_pols(r)]
        pols = [x for f in files for x in as_pols(f)]
        q = {"set": name, "length_km": L, "power_dbm": p, **point(pols, cache[key])}
        points.append(q)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps({"source": "vpi", "points": points}, indent=1), encoding="utf-8")
        corr = min(o["correlation"] for o in q["per_pol"])
        print(f"{name:13s} {str(L):>4s} km {p:3d} dBm: BER {q['ber']:.2e} | SER {q['ser']:.2e} | "
              f"SNR {q['snr_db']:6.2f} dB | girado {q['rotated_fraction']:.2f} | correlacao {corr:.2f} | "
              f"{time.perf_counter() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    main()
