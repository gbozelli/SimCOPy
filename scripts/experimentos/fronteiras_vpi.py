"""Os trechos girados de 45 graus do VPI comecam e terminam em fronteiras de bloco?

Teste da hipotese 3 das bordas (notebook 09): o bloco de frequencia do VPI tem `Ns = 1024` e
`MultipleBlockMode = Reset`. Se a frequencia ou a fase fossem processadas em blocos independentes, as
trocas de estado cairiam nas fronteiras dos blocos.

Para cada ponto em que o VPI perde a fase, alinha cada polarizacao com os simbolos conhecidos,
estima o estado de fase local (angulo da media de (y a*)^4 / 4 numa janela deslizante de 64) e acha
as trocas de estado (alinhado <-> girado) que duram pelo menos 256 simbolos dos dois lados. Depois
mede, para cada periodo N candidato, a concentracao das posicoes das trocas modulo N:
R = |media de exp(2 pi i pos / N)|. Posicoes ao acaso dao R perto de 1/sqrt(numero de trocas);
fronteiras de bloco de N simbolos dao R perto de 1.

Com --estados, em vez disso, segue o estado de fase local em passos de 45 graus (`phase_states`) em
todos os pontos dos conjuntos de configuracao conhecida, e grava as trocas de estado, a fracao do
tempo em 45 graus e a fracao em saltos de 90 ou 180 graus. E o mesmo detector usado no simcopy por
`busca_fase_vpi.py`.

Uso, na pasta repo (os dados do VPI ficam fora do repo):
  python scripts/experimentos/fronteiras_vpi.py [PASTA_VPI]             -> h3_fronteiras_vpi.json
  python scripts/experimentos/fronteiras_vpi.py --estados [PASTA_VPI]   -> estados_fase_vpi.json
(os dois em reference/vpi/hipoteses/).
"""
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
import catalogo_vpi as cat  # noqa: E402

DEFAULT_VPI = Path(r"C:\Users\bozel\Área de Trabalho\Pesquisa\25-26\dados4canais\extraido")
OUT = REPO / "reference" / "vpi" / "hipoteses" / "h3_fronteiras_vpi.json"


def transitions(x, a, window=64, hold=256):
    """Posicoes das trocas de estado (alinhado <-> girado de 45 graus), a fracao girada e a SNR [dB]
    medida so nos trechos alinhados (a SNR do bloco inteiro fica corrompida pelo giro)."""
    g = np.sum(x * np.conj(a)) / np.sum(np.abs(a) ** 2)
    c4 = (x / g * np.conj(a)) ** 4
    s = np.convolve(c4, np.ones(window), mode="same")
    rotated = np.abs(np.degrees(np.angle(s)) / 4.0) > 22.5
    # depois de um trecho girado a fase pode voltar noutro quadrante: tira o multiplo de 90 graus local
    quadrant = np.round(np.angle(np.convolve(x / g * np.conj(a), np.ones(window), mode="same")) / (np.pi / 2))
    xq = x * np.exp(-1j * quadrant * np.pi / 2)
    ok = ~rotated
    g_ok = np.sum(xq[ok] * np.conj(a[ok])) / np.sum(np.abs(a[ok]) ** 2)
    snr = float(10 * np.log10(np.mean(np.abs(a[ok]) ** 2) / np.mean(np.abs(xq[ok] / g_ok - a[ok]) ** 2)))
    change = np.flatnonzero(np.diff(rotated.astype(int)) != 0) + 1
    keep = []
    for k in change:
        before, after = rotated[max(0, k - hold):k], rotated[k:k + hold]
        if before.size < hold or after.size < hold:
            continue
        b, f = before.mean(), after.mean()
        if ((b < 0.1 and f > 0.9) or (b > 0.9 and f < 0.1)) and (not keep or k - keep[-1] > hold // 2):
            keep.append(int(k))
    return keep, float(rotated.mean()), snr


def phase_states(x, a, window=64, hold=256):
    """Estado de fase local em passos de 45 graus e as trocas de estado que duram.

    A fase local vem de sum(y a*) numa janela deslizante de `window` simbolos (com os simbolos
    conhecidos, entao nao depende do giro), desdobrada e arredondada para multiplos de 45 graus.
    Trechos de um mesmo estado com menos de `hold` simbolos sao ignorados (uma passagem rapida por
    45 graus num salto de 90 conta como uma troca so). Devolve as posicoes das trocas, a fracao do
    tempo num multiplo impar de 45 graus em relacao ao estado principal (o estado girado do VPI) e
    a fracao num multiplo par diferente do principal (saltos de 90 ou 180 graus)."""
    c = np.convolve(x * np.conj(a), np.ones(window), mode="same")
    s = np.round(np.unwrap(np.angle(c)) / (np.pi / 4)).astype(int)
    # a fase local e medida contra os simbolos conhecidos, entao a paridade e absoluta: multiplo par
    # de 45 graus e o estado alinhado (a menos de um giro de 90), impar e o girado
    even = s[s % 2 == 0]
    main = int(np.bincount(even - even.min()).argmax() + even.min()) if even.size else 0
    rel = s - main
    starts = np.r_[0, np.flatnonzero(np.diff(s) != 0) + 1]
    ends = np.r_[starts[1:], s.size]
    long = (ends - starts) >= hold
    kept_starts, kept_vals = starts[long], s[starts[long]]
    changes = kept_starts[np.flatnonzero(np.diff(kept_vals) != 0) + 1]
    return {"changes": [int(k) for k in changes], "frac_45": float(np.mean(rel % 2 == 1)),
            "frac_slipped": float(np.mean((rel % 2 == 0) & (rel != 0)))}


def states_of_all(vpi):
    """phase_states em todos os pontos dos conjuntos de configuracao conhecida, por polarizacao."""
    out = []
    for set_name, L, p, files, refs in cat.datasets(vpi):
        if not (set_name == "2026_4ch" or (set_name == "400gbps_CH2" and L == 120)):
            continue
        ref = [cat.to_grid(r) for f in refs for r in cat.as_pols(f)]
        for i, x in enumerate([x for f in files for x in cat.as_pols(f)]):
            best = None
            for r_index, a in enumerate(ref):
                peak, lag, conj = cat.best_alignment(cat.decide(x / np.sqrt(np.mean(np.abs(x) ** 2) / 10.0)), a)
                if best is None or peak > best[0]:
                    best = (peak, lag, conj, r_index)
            _, lag, conj, r_index = best
            n = min(x.size, ref[r_index].size)
            st = phase_states((np.conj(x) if conj else x)[:n], np.roll(ref[r_index][:n], lag))
            out.append({"set": set_name, "length_km": L, "power_dbm": p, "pol": i, "n": int(n), **st})
            print(f"{set_name:12s} {L} km {p:2d} dBm pol {i}: {len(st['changes'])} trocas, "
                  f"{100 * st['frac_45']:.1f} % em 45 graus, {100 * st['frac_slipped']:.1f} % saltado", flush=True)
    return out


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    if "--estados" in sys.argv:
        args = [a for a in sys.argv[1:] if a != "--estados"]
        vpi = Path(args[0]) if args else DEFAULT_VPI
        out = REPO / "reference" / "vpi" / "hipoteses" / "estados_fase_vpi.json"
        out.write_text(json.dumps(states_of_all(vpi)), encoding="utf-8")
        return
    vpi = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_VPI
    catalog = json.loads((REPO / "reference/vpi/catalogo/catalogo_vpi.json").read_text(encoding="utf-8"))["points"]
    lost = {(p["set"], p["length_km"], p["power_dbm"]) for p in catalog
            if p["rotated_fraction"] > 0.02 or p["ber"] > 0.1}
    records = []
    for set_name, L, p, files, refs in cat.datasets(vpi):
        if (set_name, L, p) not in lost or set_name not in ("2026_4ch", "400gbps_CH2", "2025_4ch", "2025_1ch"):
            continue
        if set_name == "400gbps_CH2" and L == 140:      # duplicata de 2026_4ch 140 km
            continue
        pols = [x for f in files for x in cat.as_pols(f)]
        ref = [cat.to_grid(r) for f in refs for r in cat.as_pols(f)]
        for i, x in enumerate(pols):
            best = None
            for r_index, a in enumerate(ref):
                peak, lag, conj = cat.best_alignment(cat.decide(x / np.sqrt(np.mean(np.abs(x) ** 2) / 10.0)), a)
                if best is None or peak > best[0]:
                    best = (peak, lag, conj, r_index)
            peak, lag, conj, r_index = best
            n = min(x.size, ref[r_index].size)
            xx = (np.conj(x) if conj else x)[:n]
            a = np.roll(ref[r_index][:n], lag)
            pos, frac, snr = transitions(xx, a)
            records.append({"set": set_name, "length_km": L, "power_dbm": p, "pol": i, "n": int(n),
                            "rotated": frac, "transitions": pos, "lag": int(lag), "snr_aligned_db": snr})
            print(f"{set_name:12s} {L} km {p:2d} dBm pol {i}: girado {frac:.2f}, {len(pos)} trocas, "
                  f"SNR nos trechos alinhados {snr:.2f} dB", flush=True)

    OUT.write_text(json.dumps(records), encoding="utf-8")
    use = [np.array(r["transitions"]) for r in records if len(r["transitions"]) >= 3]
    print(f"\n{len(use)} polarizacoes com 3 trocas ou mais, {sum(t.size for t in use)} trocas")
    print(f"soma de n R^2 (ao acaso: media {len(use)}, desvio {np.sqrt(len(use)):.1f})")
    for N in (128, 256, 512, 1000, 1024, 2048, 4096, 8192):
        z = sum(abs(np.sum(np.exp(2j * np.pi * t / N))) ** 2 / t.size for t in use)
        print(f"  N = {N:5d}: {z:.1f}")


if __name__ == "__main__":
    main()
