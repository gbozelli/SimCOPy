"""Teste de mutacao: quebra o codigo de proposito, um defeito por vez, e verifica
se algum teste percebe. Mutacao NAO DETECTADA = propriedade sem protecao.

Serve para revisar codigo que voce nao escreveu: em vez de ler linha a linha
perguntando "isso esta certo?", pergunte "se isso estivesse errado, eu saberia?".

Uso (na raiz do repositorio):
    python scripts/mutacoes.py            # todas
    python scripts/mutacoes.py 0 4 5      # so algumas
    python scripts/mutacoes.py --listar

Cada mutacao e restaurada ao final, mesmo se o script for interrompido.
Para criar uma nova, acrescente uma linha em MUTACOES: (descricao, arquivo,
trecho original, trecho mutado). O trecho original precisa existir literalmente.
"""
import pathlib, subprocess, sys

RAIZ = pathlib.Path(__file__).resolve().parents[1]

MUTACOES = [
    ("Sinal de beta2 invertido",
     "src/simcopy/core/units.py",
     "return -D * wavelength**2 / (2.0 * np.pi * C_LIGHT)",
     "return D * wavelength**2 / (2.0 * np.pi * C_LIGHT)"),
    ("Sem o fator 8/9 de Manakov",
     "src/simcopy/kernels/ssfm.py",
     "MANAKOV_FACTOR = 8.0 / 9.0", "MANAKOV_FACTOR = 1.0"),
    ("Sinal da nao linearidade invertido",
     "src/simcopy/kernels/ssfm.py",
     "nl = np.exp(1j * p.gamma", "nl = np.exp(-1j * p.gamma"),
    ("Passo nao linear sem comprimento efetivo",
     "src/simcopy/kernels/ssfm.py",
     "dz_nl = dz if p.alpha_np == 0 else", "dz_nl = dz if True else"),
    ("DGD por passo volta a D_PMD*sqrt(passo) (bug ja corrigido)",
     "src/simcopy/kernels/birefringence.py",
     "return pmd_coefficient * step / np.sqrt(max(step, correlation_length))",
     "return pmd_coefficient * np.sqrt(step)"),
    ("Rotacao de SoP com 2 angulos, como no rascunho",
     "src/simcopy/kernels/birefringence.py",
     "    g = rng.standard_normal(4)\n",
     "    th, ph = rng.random() * np.pi, rng.random() * 2 * np.pi\n"
     "    D = np.diag([np.exp(.5j * ph), np.exp(-.5j * ph)])\n"
     "    return np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]]) @ D\n"
     "    g = rng.standard_normal(4)\n"),
    ("SSFM assimetrico (primeira ordem)",
     "src/simcopy/kernels/ssfm.py",
     "        Ax, Ay = _linear_half(Ax, Ay, omega, p.beta2, p.alpha_np, dz, dtau)",
     "        Ax, Ay = _linear_half(Ax, Ay, omega, p.beta2, p.alpha_np, 2 * dz, dtau)\n"
     "        _ASSIM = True"),
    ("lambda fixo em 1550 nm no resolve",
     "src/simcopy/config/resolve.py",
     "    lam = frequency_to_wavelength(center_frequency)", "    lam = 1550e-9"),
    ("IMDD deixa de recusar formatos coerentes",
     "src/simcopy/config/system.py",
     "        if self.modulation.name not in IMDD_FORMATS:", "        if False:"),
    ("Potencia soma so a polarizacao X (bug C-5 do rascunho)",
     "src/simcopy/core/signals.py",
     "return float(np.sum(np.mean(np.abs(self.samples) ** 2, axis=1)))",
     "return float(2 * np.mean(np.abs(self.samples[0]) ** 2))"),
]

def _aplica(i):
    nome, arq, a, b = MUTACOES[i]
    f = RAIZ / arq
    orig = f.read_text()
    if a not in orig:
        return nome, None, "trecho original nao encontrado (o codigo mudou?)"
    mut = orig.replace(a, b, 1)
    if "_ASSIM = True" in b:    # assimetrico: remove o segundo meio passo linear
        mut = mut.replace(
            "        Ax, Ay = _linear_half(Ax, Ay, omega, p.beta2, p.alpha_np, dz)\n",
            "", 1).replace("        _ASSIM = True\n", "")
    try:
        f.write_text(mut)
        r = subprocess.run([sys.executable, "-m", "pytest", "-q", "--no-header",
                            "-p", "no:cacheprovider"], cwd=RAIZ,
                           capture_output=True, text=True)
        falhas = [l.split(" ")[1] for l in r.stdout.splitlines() if l.startswith("FAILED")]
        return nome, falhas, None
    finally:
        f.write_text(orig)

def main():
    if "--listar" in sys.argv:
        for i, m in enumerate(MUTACOES):
            print(f"[{i}] {m[0]}")
        return 0
    sel = [int(x) for x in sys.argv[1:]] or range(len(MUTACOES))
    lacunas = 0
    for i in sel:
        nome, falhas, erro = _aplica(i)
        if erro:
            print(f"[{i}] {nome}\n     ?? {erro}"); lacunas += 1; continue
        if falhas:
            print(f"[{i}] {nome}\n     detectada por {len(falhas)} teste(s), ex.: {falhas[0]}")
        else:
            print(f"[{i}] {nome}\n     >>> NAO DETECTADA: nenhum teste protege isto <<<")
            lacunas += 1
    print(f"\n{len(list(sel)) - lacunas}/{len(list(sel))} mutacoes detectadas.")
    return 1 if lacunas else 0

if __name__ == "__main__":
    sys.exit(main())
