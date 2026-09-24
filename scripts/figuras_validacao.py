"""Figuras de validacao: solucao ANALITICA contra a NOSSA IMPLEMENTACAO.

Estrutura comum a todas as figuras: a forma fechada da literatura e a saida do
nosso SSFM sao plotadas juntas, variando um parametro por vez, com um painel
dedicado ao erro entre as duas. As figuras 01 a 04 isolam um fenomeno cada; a 05
empilha os fenomenos.

As assercoes numericas estao em tests/. A explicacao painel a painel esta em
docs/figuras.md.

Uso:  python scripts/figuras_validacao.py [pasta_saida]
"""
import sys, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
RAIZ = pathlib.Path(__file__).resolve().parents[1]
sys.path[:0] = [str(RAIZ / "src"), str(RAIZ)]
from simcopy.kernels.ssfm import propagate, FiberKernelParams
from simcopy.core.units import (beta2_from_dispersion, db_per_m_to_np_per_m,
                                effective_length, gamma_from_n2)
from tests.physics._helpers import fwhm, jones_dgd

OUT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "figuras")
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"figure.dpi": 120, "font.size": 9, "axes.grid": True,
                     "grid.alpha": .3, "legend.fontsize": 7.5})
C3 = ["tab:blue", "tab:orange", "tab:green"]
C4 = C3 + ["tab:red"]
ERRO = lambda a, s: np.max(np.abs(s - a)) / np.max(np.abs(a))
LAM = 1550e-9

def salvar(fig, nome, formula):
    fig.suptitle(formula, fontsize=9.5, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT / nome)
    plt.close(fig)
    print("  ", nome)

# ============================================================ 01 ATENUACAO ====
def fig01():
    fig, ax = plt.subplots(2, 2, figsize=(9.5, 6.4))
    zs = np.linspace(2e3, 150e3, 12)
    zz = np.linspace(0, 150e3, 300)
    N, FS = 256, 1e12
    OM = 2 * np.pi * np.fft.fftfreq(N, 1 / FS)
    for i, adb in enumerate([0.15e-3, 0.20e-3, 0.25e-3]):
        a = db_per_m_to_np_per_m(adb)
        sim = np.array([np.mean(np.abs(propagate(np.ones((1, N), complex), OM,
              FiberKernelParams(z, 0, 0, a)))**2) for z in zs])
        ax[0, 0].plot(zz / 1e3, 10 * np.log10(np.exp(-a * zz)), "-", color=C3[i],
                      label=rf"analítico, $\alpha$ = {adb*1e3:.2f} dB/km")
        ax[0, 0].plot(zs / 1e3, 10 * np.log10(sim), "o", ms=4, mfc="none", color=C3[i],
                      label="simulado" if i == 0 else None)
        ax[0, 1].semilogy(zs / 1e3, np.abs(sim - np.exp(-a * zs)) / np.exp(-a * zs),
                          "o-", ms=3, color=C3[i], label=rf"{adb*1e3:.2f} dB/km")
    ax[0, 0].set(xlabel="z [km]", ylabel="P(z)/P(0) [dB]",
                 title=r"(a) decaimento para três $\alpha$")
    ax[0, 0].legend()
    ax[0, 1].axhline(2.2e-16, color="gray", ls=":", lw=.8)
    ax[0, 1].set(xlabel="z [km]", ylabel="erro relativo", title="(b) erro")
    ax[0, 1].legend(title=r"$\alpha$")

    # (c) a perda nao pode depender da grade nem da forma do pulso
    a = db_per_m_to_np_per_m(0.2e-3)
    casos = [("CW, N=256, fs=1 THz", 256, 1e12, None),
             ("CW, N=4096, fs=4 THz", 4096, 4e12, None),
             ("gaussiana 10 ps", 4096, 4e12, 10e-12),
             ("gaussiana 50 ps", 4096, 4e12, 50e-12)]
    for i, (nome, n, fs, t0) in enumerate(casos):
        t = (np.arange(n) - n // 2) / fs
        A = (np.ones(n, complex) if t0 is None
             else np.exp(-t**2 / (2 * t0**2)).astype(complex))[None, :]
        om = 2 * np.pi * np.fft.fftfreq(n, 1 / fs)
        e = [abs(np.sum(np.abs(propagate(A, om, FiberKernelParams(z, 0, 0, a)))**2)
                 / np.sum(np.abs(A)**2) / np.exp(-a * z) - 1) for z in zs]
        ax[1, 0].semilogy(zs / 1e3, np.maximum(e, 1e-17), "o-", ms=3, color=C4[i], label=nome)
    ax[1, 0].set(xlabel="z [km]", ylabel="erro relativo",
                 title=r"(c) invariância: a perda só depende de $\alpha z$")
    ax[1, 0].legend()

    # (d) a atenuacao e exata por passo, logo nao depende do passo
    passos = np.array([50., 200., 1e3, 5e3, 2e4])
    e = [abs(np.mean(np.abs(propagate(np.ones((1, 64), complex),
         2 * np.pi * np.fft.fftfreq(64, 1e-12), FiberKernelParams(100e3, 0, 0, a),
         max_step=h)[0])**2) / np.exp(-a * 100e3) - 1) for h in passos]
    ax[1, 1].loglog(passos / 1e3, np.maximum(e, 1e-17), "o-", ms=4, color=C3[0])
    ax[1, 1].set(xlabel="passo do SSFM [km]", ylabel="erro relativo", ylim=(1e-17, 1e-12),
                 title="(d) independência do passo (operador exato)")
    salvar(fig, "01_atenuacao.png",
           r"01. Atenuação — analítico: $P(z) = P(0)\,e^{-\alpha z}$")

# ============================================================ 02 DISPERSAO ====
N2, FS2 = 1 << 13, 1e12
T2 = (np.arange(N2) - N2 // 2) / FS2
OM2 = 2 * np.pi * np.fft.fftfreq(N2, 1 / FS2)

def gauss_ana(T, T0, z, b2, C=0.0):
    """Agrawal, NLFO 3.2: forma fechada da gaussiana com chirp em meio dispersivo."""
    q = T0**2 - 1j * b2 * z * (1 + 1j * C)
    return np.sqrt(T0**2 / q) * np.exp(-(1 + 1j * C) * T**2 / (2 * q))

def fig02():
    fig, ax = plt.subplots(2, 2, figsize=(9.5, 6.4))
    T0 = 20e-12
    B2 = beta2_from_dispersion(16e-6, LAM)
    LD = T0**2 / abs(B2)
    sim = lambda z, C, b2=B2, t0=T0: propagate(
        gauss_ana(T2, t0, 0, b2, C)[None, :], OM2,
        FiberKernelParams(z, b2, 0.0, 0.0), max_step=z / 200)[0]
    # (a) perfis em 0,5 L_D para tres chirps
    for i, C in enumerate([-2.0, 0.0, 2.0]):
        za = 0.5 * LD
        a, s = gauss_ana(T2, T0, za, B2, C), sim(za, C)
        ax[0, 0].plot(T2 * 1e12, np.abs(a)**2, "-", color=C3[i], label=f"analítico, C = {C:+.0f}")
        ax[0, 0].plot(T2[::22] * 1e12, np.abs(s[::22])**2, "o", ms=3, mfc="none",
                      color=C3[i], label="simulado" if i == 0 else None)
    ax[0, 0].set(xlim=(-140, 140), xlabel="T [ps]", ylabel=r"$|A|^2$ [u.a.]",
                 title=rf"(a) perfil em z = 0,5 $L_D$ ({0.5*LD/1e3:.1f} km), $T_0$ = 20 ps")
    ax[0, 0].legend()
    # (b) variando D com T0 fixo
    for i, D in enumerate([8e-6, 16e-6, 22e-6]):
        b2 = beta2_from_dispersion(D, LAM)
        ld = T0**2 / abs(b2)
        zz = np.linspace(0.02, 2, 60) * (T0**2 / abs(beta2_from_dispersion(8e-6, LAM)))
        w0 = fwhm(T2, np.abs(gauss_ana(T2, T0, 0, b2))**2)
        ax[0, 1].plot(zz / 1e3, np.sqrt(1 + (zz / ld)**2), "-", color=C3[i],
                      label=rf"analítico, D = {D*1e6:.0f} ps/(nm·km)")
        zp = np.linspace(0.15, 2, 7) * (T0**2 / abs(beta2_from_dispersion(8e-6, LAM)))
        ax[0, 1].plot(zp / 1e3, [fwhm(T2, np.abs(sim(z, 0.0, b2))**2) / w0 for z in zp],
                      "o", ms=4, mfc="none", color=C3[i], label="simulado" if i == 0 else None)
    ax[0, 1].set(xlabel="z [km]", ylabel="FWHM(z)/FWHM(0)",
                 title=r"(b) variando D: $L_D = T_0^2/|\beta_2|$ encurta")
    ax[0, 1].legend()
    # (c) colapso: o alargamento so depende de z/L_D
    combos = [(10e-12, 8e-6), (20e-12, 16e-6), (40e-12, 22e-6), (5e-12, 16e-6)]
    xx = np.linspace(0, 2.5, 200)
    ax[1, 0].plot(xx, np.sqrt(1 + xx**2), "k-", lw=2, label=r"analítico, $\sqrt{1+(z/L_D)^2}$")
    for i, (t0, D) in enumerate(combos):
        b2 = beta2_from_dispersion(D, LAM)
        ld = t0**2 / abs(b2)
        w0 = fwhm(T2, np.abs(gauss_ana(T2, t0, 0, b2))**2)
        zp = np.linspace(0.15, 2.5, 7) * ld
        ax[1, 0].plot(zp / ld, [fwhm(T2, np.abs(sim(z, 0.0, b2, t0))**2) / w0 for z in zp],
                      "o", ms=4, mfc="none", color=C4[i],
                      label=rf"$T_0$={t0*1e12:.0f} ps, D={D*1e6:.0f} ($L_D$={ld/1e3:.1f} km)")
    ax[1, 0].set(xlabel=r"$z / L_D$", ylabel="FWHM(z)/FWHM(0)",
                 title="(c) colapso: quatro pares $(T_0, D)$ sobre uma curva")
    ax[1, 0].legend()
    # (d) erro
    for i, (t0, D) in enumerate(combos):
        b2 = beta2_from_dispersion(D, LAM)
        ld = t0**2 / abs(b2)
        zp = np.linspace(0.15, 2.5, 7) * ld
        ax[1, 1].semilogy(zp / ld, [ERRO(gauss_ana(T2, t0, z, b2), sim(z, 0.0, b2, t0))
                                    for z in zp], "o-", ms=3, color=C4[i],
                          label=rf"$T_0$={t0*1e12:.0f} ps, D={D*1e6:.0f}")
    ax[1, 1].set(xlabel=r"$z / L_D$", ylabel="erro relativo no campo", title="(d) erro")
    ax[1, 1].legend()
    salvar(fig, "02_dispersao.png",
           r"02. Dispersão cromática ($\gamma=0$, $\alpha=0$) — analítico: forma fechada da "
           r"gaussiana com chirp, $T_1/T_0=\sqrt{(1+C\beta_2 z/T_0^2)^2+(\beta_2 z/T_0^2)^2}$")

# ===================================================== 03 NAO LINEARIDADE =====
def fig03():
    fig, ax = plt.subplots(2, 3, figsize=(13.5, 6.4))
    T0, L = 20e-12, 50e3
    ALPHA = db_per_m_to_np_per_m(0.2e-3)
    G0 = gamma_from_n2(2.6e-20, 80e-12, LAM)
    LEFF = effective_length(L, 0.2e-3)
    A0f = lambda P0: np.sqrt(P0) * np.exp(-T2**2 / (2 * T0**2))
    def ana(P0, z, g=G0, al=ALPHA):
        A0 = A0f(P0)
        le = z if al == 0 else (1 - np.exp(-al * z)) / al
        return A0 * np.exp(-al * z / 2) * np.exp(1j * g * np.abs(A0)**2 * le)
    def sim(P0, z, g=G0, al=ALPHA):
        return propagate(A0f(P0)[None, :], OM2,
                         FiberKernelParams(z, 0.0, g, al, manakov=False), max_step=1e3)[0]
    def phi(P0, out):
        A0 = A0f(P0)
        m = np.abs(A0) > 1e-6 * np.abs(A0).max()
        p = np.zeros(N2)
        p[m] = np.unwrap(np.angle(out[m] * np.conj(A0[m])))
        return p
    fases = [0.5 * np.pi, 1.5 * np.pi, 3.5 * np.pi]
    for i, pm in enumerate(fases):
        P0 = pm / (G0 * LEFF)
        a, s = ana(P0, L), sim(P0, L)
        ax[0, 0].plot(T2 * 1e12, phi(P0, a) / np.pi, "-", color=C3[i],
                      label=rf"analítico, $\phi_{{max}}$ = {pm/np.pi:.1f}$\pi$")
        ax[0, 0].plot(T2[::22] * 1e12, phi(P0, s)[::22] / np.pi, "o", ms=3, mfc="none",
                      color=C3[i], label="simulado" if i == 0 else None)
        f = np.fft.fftshift(np.fft.fftfreq(N2, 1 / FS2))
        Sa = np.abs(np.fft.fftshift(np.fft.fft(a)))**2
        Ss = np.abs(np.fft.fftshift(np.fft.fft(s)))**2
        ax[0, 1].plot(f / 1e9, 10 * np.log10(Sa / Sa.max() + 1e-12), "-", color=C3[i])
        ax[0, 1].plot(f[::40] / 1e9, 10 * np.log10(Ss[::40] / Sa.max() + 1e-12), "o",
                      ms=2.5, mfc="none", color=C3[i])
        zz = np.linspace(0.05, 1, 8) * L
        ax[0, 2].semilogy(zz / 1e3, [ERRO(ana(P0, z), sim(P0, z)) for z in zz],
                          "o-", ms=3, color=C3[i], label=rf"$\phi_{{max}}$={pm/np.pi:.1f}$\pi$")
    ax[0, 0].set(xlim=(-70, 70), ylim=(-0.2, 4.4), xlabel="T [ps]",
                 ylabel=r"$\phi_{NL}$ [$\pi$ rad]", title="(a) fase não linear")
    ax[0, 0].legend()
    ax[0, 1].set(xlim=(-120, 120), ylim=(-45, 5), xlabel="frequência [GHz]",
                 ylabel="espectro [dB]", title=r"(b) espectro: nº de picos $\approx \phi_{max}/\pi + 1/2$")
    ax[0, 2].set(xlabel="z [km]", ylabel="erro relativo no campo", title="(c) erro")
    ax[0, 2].legend()
    # (d) variando alpha: saturacao em L_eff = 1/alpha
    P0 = 1.5 * np.pi / (G0 * LEFF)
    Ls = np.linspace(5e3, 200e3, 9)
    for i, adb in enumerate([0.0, 0.2e-3, 0.35e-3]):
        al = db_per_m_to_np_per_m(adb)
        LL = np.linspace(1e3, 200e3, 200)
        le = LL if al == 0 else (1 - np.exp(-al * LL)) / al
        ax[1, 0].plot(LL / 1e3, G0 * P0 * le / np.pi, "-", color=C3[i],
                      label=rf"analítico, $\alpha$={adb*1e3:.2f} dB/km")
        ax[1, 0].plot(Ls / 1e3, [phi(P0, sim(P0, z, al=al)).max() / np.pi for z in Ls],
                      "o", ms=4, mfc="none", color=C3[i], label="simulado" if i == 0 else None)
        if al:
            ax[1, 0].axhline(G0 * P0 / al / np.pi, color=C3[i], ls=":", lw=.8)
    ax[1, 0].set(xlabel="L [km]", ylabel=r"$\phi_{max}$ [$\pi$ rad]",
                 title=r"(d) variando $\alpha$: saturação em $L_{eff}\to1/\alpha$")
    ax[1, 0].legend()
    # (e) variando A_eff: phi_max proporcional a gamma
    Ps = np.linspace(0.02, 0.35, 8)
    for i, aeff in enumerate([60e-12, 80e-12, 110e-12]):
        g = gamma_from_n2(2.6e-20, aeff, LAM)
        PP = np.linspace(0, 0.36, 100)
        ax[1, 1].plot(PP * 1e3, g * PP * LEFF / np.pi, "-", color=C3[i],
                      label=rf"analítico, $A_{{eff}}$={aeff*1e12:.0f} µm² ($\gamma$={g*1e3:.2f})")
        ax[1, 1].plot(Ps * 1e3, [phi(P, sim(P, L, g=g)).max() / np.pi for P in Ps],
                      "o", ms=4, mfc="none", color=C3[i], label="simulado" if i == 0 else None)
    ax[1, 1].set(xlabel=r"$P_0$ [mW]", ylabel=r"$\phi_{max}$ [$\pi$ rad]",
                 title=r"(e) variando $A_{eff}$: $\phi_{max} = \gamma P_0 L_{eff}$")
    ax[1, 1].legend()
    # (f) erro dos painéis (d) e (e)
    for i, adb in enumerate([0.0, 0.2e-3, 0.35e-3]):
        al = db_per_m_to_np_per_m(adb)
        ax[1, 2].semilogy(Ls / 1e3, [ERRO(ana(P0, z, al=al), sim(P0, z, al=al)) for z in Ls],
                          "o-", ms=3, color=C3[i], label=rf"$\alpha$={adb*1e3:.2f} dB/km")
    ax[1, 2].set(xlabel="L [km]", ylabel="erro relativo no campo", title="(f) erro do painel (d)")
    ax[1, 2].legend()
    salvar(fig, "03_nao_linearidade.png",
           r"03. Automodulação de fase ($\beta_2=0$, propagação escalar) — analítico: "
           r"$A(z,T)=A(0,T)e^{-\alpha z/2}\exp[i\gamma|A(0,T)|^2 L_{eff}]$, "
           r"$L_{eff}=(1-e^{-\alpha z})/\alpha$")

# ================================================================= 04 PMD =====
DF, NS = 1.25e9, 64
OMP = 2 * np.pi * np.fft.fftfreq(NS, 1 / (NS * DF))

def sigma_ana(z, dp, lc):
    """Agrawal, FOCS cap. 2: sigma_T^2 = 2(dbeta1*lc)^2[exp(-z/lc)+z/lc-1]."""
    db1 = dp / np.sqrt(2 * lc)
    return np.sqrt(2 * (db1 * lc)**2 * (np.exp(-z / lc) + z / lc - 1))

def sigma_sim(L, dp, lc, n=60):
    # O passo so precisa resolver a secao quando z e comparavel a l_c; para
    # z >> l_c o modelo de passo grosso agrega secoes e o passo pode ser maior.
    st = min(max(L / 20, 0.5), 1e3)
    if L < 20 * lc:
        st = min(st, lc / 4)
    tau = jones_dgd(lambda A, r: propagate(A, OMP, FiberKernelParams(L, 0, 0, 0, dp, lc),
                    rng=np.random.default_rng(40_000 + r), max_step=st), n, DF, NS)
    return np.sqrt(np.mean(tau**2)), tau

def fig04():
    fig, ax = plt.subplots(2, 2, figsize=(9.5, 6.4))
    zs = np.array([50., 200., 1e3, 1e4, 8e4])
    zz = np.logspace(1, 5, 200)
    # (a) variando D_PMD, l_c fixo
    for i, dpk in enumerate([0.05, 0.1, 0.2]):
        dp = dpk * 1e-12 / np.sqrt(1e3)
        ax[0, 0].loglog(zz, sigma_ana(zz, dp, 50.) * 1e12, "-", color=C3[i],
                        label=rf"analítico, $D_p$={dpk} ps/$\sqrt{{km}}$")
        s = np.array([sigma_sim(z, dp, 50.)[0] for z in zs])
        ax[0, 0].loglog(zs, s * 1e12, "o", ms=4, mfc="none", color=C3[i],
                        label="simulado" if i == 0 else None)
        ax[1, 0].semilogx(zs, s / sigma_ana(zs, dp, 50.), "o-", ms=3, color=C3[i],
                          label=rf"$D_p$={dpk}")
    ax[0, 0].axvline(50., color="gray", lw=.8)
    ax[0, 0].set(xlabel="z [m]", ylabel=r"$\sigma_T$ [ps]",
                 title=r"(a) variando $D_p$ ($l_c$ = 50 m)")
    ax[0, 0].legend(loc="lower right")
    # (b) variando l_c, D_PMD fixo
    dp = 0.1e-12 / np.sqrt(1e3)
    for i, lc in enumerate([25., 50., 100.]):
        ax[0, 1].loglog(zz, sigma_ana(zz, dp, lc) * 1e12, "-", color=C3[i],
                        label=rf"analítico, $l_c$={lc:.0f} m")
        s = np.array([sigma_sim(z, dp, lc)[0] for z in zs])
        ax[0, 1].loglog(zs, s * 1e12, "o", ms=4, mfc="none", color=C3[i],
                        label="simulado" if i == 0 else None)
        ax[1, 0].semilogx(zs, s / sigma_ana(zs, dp, lc), "s--", ms=3, color=C3[i], alpha=.5)
    ax[0, 1].loglog(zz, dp * np.sqrt(zz) * 1e12, "k:", lw=1, label=r"$D_p\sqrt{z}$")
    ax[0, 1].set(xlabel="z [m]", ylabel=r"$\sigma_T$ [ps]",
                 title=r"(b) variando $l_c$ ($D_p$ = 0,1): só muda a transição")
    ax[0, 1].legend(loc="lower right")
    ax[1, 0].axhline(1, color="k", lw=.8)
    ax[1, 0].axhline(np.sqrt(2), color="gray", ls=":", lw=.8)
    ax[1, 0].text(12, 1.46, r"$\sqrt{2}$", fontsize=8, color="gray")
    ax[1, 0].set(xlabel="z [m]", ylim=(0.6, 1.75), ylabel="simulado / analítico",
                 title=r"(c) razão (círculos: varia $D_p$; quadrados: varia $l_c$)")
    ax[1, 0].legend(ncol=3)
    # (d) distribuicao
    _, tau = sigma_sim(2e4, dp, 50., n=800)
    s3 = np.sqrt(np.mean(tau**2) / 3)
    x = np.linspace(0, tau.max(), 200)
    ax[1, 1].hist(tau * 1e12, bins=40, density=True, alpha=.45, color=C3[1],
                  label="simulado, 20 km")
    ax[1, 1].plot(x * 1e12, np.sqrt(2 / np.pi) * x**2 / s3**3
                  * np.exp(-x**2 / (2 * s3**2)) * 1e-12, "k-", label="Maxwelliana")
    ax[1, 1].set(xlabel="DGD [ps]", ylabel="densidade",
                 title=rf"(d) distribuição ($\langle\tau\rangle/\tau_{{rms}}$ = "
                       rf"{np.mean(tau)/np.sqrt(np.mean(tau**2)):.3f}; teórico 0,921)")
    ax[1, 1].legend()
    salvar(fig, "04_pmd.png",
           r"04. PMD — analítico: $\sigma_T^2 = 2(\Delta\beta_1 l_c)^2"
           r"[e^{-z/l_c} + z/l_c - 1]$, com $D_p = \Delta\beta_1\sqrt{2l_c}$")

# ========================================================== 05 EMPILHANDO =====
def fig05():
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7))
    T0 = 20e-12
    B2 = beta2_from_dispersion(16e-6, LAM)
    LD = T0**2 / abs(B2)
    ALPHA = db_per_m_to_np_per_m(0.2e-3)
    G0 = gamma_from_n2(2.6e-20, 80e-12, LAM)
    LEFF = effective_length(50e3, 0.2e-3)
    A0f = lambda P0: np.sqrt(P0) * np.exp(-T2**2 / (2 * T0**2))
    zz = np.linspace(0.1, 2, 8) * LD
    etapas = [(r"1. só $\alpha$", zz, [ERRO(np.exp(-ALPHA * z / 2) * A0f(1.0),
        propagate(A0f(1.0)[None, :], OM2, FiberKernelParams(z, 0, 0, ALPHA))[0]) for z in zz])]
    etapas.append((r"2. $\alpha + \beta_2$", zz, [ERRO(
        gauss_ana(T2, T0, z, B2) * np.exp(-ALPHA * z / 2),
        propagate(gauss_ana(T2, T0, 0, B2)[None, :], OM2,
                  FiberKernelParams(z, B2, 0.0, ALPHA), max_step=z / 200)[0]) for z in zz]))
    P0 = 1.5 * np.pi / (G0 * LEFF)
    def spm_ana(z):
        A0 = A0f(P0)
        le = (1 - np.exp(-ALPHA * z)) / ALPHA
        return A0 * np.exp(-ALPHA * z / 2) * np.exp(1j * G0 * np.abs(A0)**2 * le)
    etapas.append((r"3. $\alpha + \gamma$", zz, [ERRO(spm_ana(z),
        propagate(A0f(P0)[None, :], OM2,
                  FiberKernelParams(z, 0.0, G0, ALPHA, manakov=False))[0])
        for z in zz]))
    T0s = 10e-12
    n5 = 1 << 11
    FSs = n5 / (80 * T0s)
    Ts = (np.arange(n5) - n5 // 2) / FSs
    OMs = 2 * np.pi * np.fft.fftfreq(n5, 1 / FSs)
    LDs = T0s**2 / abs(B2)
    P0s = abs(B2) / (G0 * T0s**2)
    Z0 = np.pi / 2 * LDs
    A0s = (np.sqrt(P0s) / np.cosh(Ts / T0s)).astype(complex)[None, :]
    zsol = np.linspace(0.25, 4, 8) * Z0
    etapas.append((r"4. $\beta_2 + \gamma$ (sóliton)", zsol / Z0 * (np.pi / 2),
        [ERRO(A0s[0] * np.exp(1j * z / (2 * LDs)),
              propagate(A0s, OMs, FiberKernelParams(z, B2, G0, 0.0, manakov=False))[0])
         for z in zsol]))
    for i, (nome, x, y) in enumerate(etapas):
        ax[0, 0].semilogy(np.asarray(x) / (LD if i < 3 else LDs) * (LD if i < 3 else LDs)
                          / (LD if i < 3 else LDs), y, "o-", ms=3, color=C4[i], label=nome)
    ax[0, 0].set(xlabel=r"$z/L_D$ (o $L_D$ de cada etapa)", ylabel="erro relativo no campo",
                 title="(a) erro contra a forma fechada, fenômeno a fenômeno")
    ax[0, 0].legend()
    ax[0, 1].plot(Ts * 1e12, np.abs(A0s[0])**2 * 1e3, "k-", label="analítico (invariante)")
    for i, n in enumerate([1, 2, 4]):
        o = propagate(A0s, OMs, FiberKernelParams(n * Z0, B2, G0, 0.0, manakov=False))[0]
        ax[0, 1].plot(Ts[::12] * 1e12, np.abs(o[::12])**2 * 1e3, "o", ms=2.5, mfc="none",
                      color=C3[i], label=rf"simulado, {n} $z_0$")
    ax[0, 1].set(xlim=(-45, 45), xlabel="T [ps]", ylabel=r"$|A|^2$ [mW]",
                 title=rf"(b) $\beta_2+\gamma$: sóliton fundamental ($z_0$={Z0/1e3:.1f} km)")
    ax[0, 1].legend()
    ref = propagate(A0s, OMs, FiberKernelParams(2 * Z0, B2, G0, 0.0, manakov=False),
                    max_step=LDs / 2000, max_phase_deg=1e9, min_step=LDs / 2000)[0]
    hs = np.array([LDs / 4, LDs / 8, LDs / 16, LDs / 32, LDs / 64])
    er = [ERRO(ref, propagate(A0s, OMs, FiberKernelParams(2 * Z0, B2, G0, 0.0, manakov=False),
          max_step=h, max_phase_deg=1e9, min_step=h)[0]) for h in hs]
    ax[1, 0].loglog(hs / LDs, er, "o-", ms=4, color=C3[0], label="SSFM simétrico")
    ax[1, 0].loglog(hs / LDs, er[0] * (hs / hs[0])**2, "k--", label=r"$\propto h^2$")
    ax[1, 0].loglog(hs / LDs, er[0] * (hs / hs[0]), "k:", label=r"$\propto h$")
    ax[1, 0].set(xlabel=r"passo $h/L_D$", ylabel="erro vs referência de passo fino",
                 title="(c) sem forma fechada: ordem de convergência")
    ax[1, 0].legend()
    zs4 = np.linspace(10e3, 150e3, 10)
    A2 = np.vstack([A0s[0], A0s[0] * 0.7])
    e0 = np.sum(np.abs(A2)**2)
    en = [np.sum(np.abs(propagate(A2, OMs, FiberKernelParams(z, B2, G0, 0.0,
          0.1e-12 / np.sqrt(1e3), 50.), rng=np.random.default_rng(5)))**2) for z in zs4]
    ax[1, 1].semilogy(zs4 / 1e3, np.abs(np.array(en) - e0) / e0, "o-", ms=4, color=C3[2])
    ax[1, 1].axhline(2.2e-16, color="gray", ls=":", lw=.8)
    ax[1, 1].set(xlabel="z [km]", ylabel="desvio relativo de energia",
                 title=r"(d) $\beta_2+\gamma+$PMD, DP, $\alpha=0$: energia conservada")
    salvar(fig, "05_empilhando.png",
           "05. Empilhando fenômenos — onde existe forma fechada, compara-se com ela; "
           "onde não existe, usam-se ordem de convergência e invariantes")

print("gerando:")
_sel = [a for a in sys.argv[2:]] or ["1", "2", "3", "4", "5"]
for _k in _sel:
    {"1": fig01, "2": fig02, "3": fig03, "4": fig04, "5": fig05}[_k]()
print("figuras em", OUT.resolve())
