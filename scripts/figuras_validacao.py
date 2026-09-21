"""Gera as figuras de validacao fisica. Os testes em tests/physics fazem as
assercoes numericas; este script produz as imagens para inspecao visual e
documentacao. Uso: python scripts/figuras_validacao.py [pasta_saida]"""
import sys, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from simcopy.kernels.ssfm import propagate, FiberKernelParams
from simcopy.core.units import beta2_from_dispersion, db_per_m_to_np_per_m, effective_length
from tests.physics._helpers import gaussian, fwhm, jones_dgd

OUT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "figuras")
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"figure.dpi": 120, "font.size": 10, "axes.grid": True, "grid.alpha": .3})

# ---------------------------------------------------------------- 1. atenuacao
a = db_per_m_to_np_per_m(0.2e-3)
om = 2*np.pi*np.fft.fftfreq(64, 1e-12)
zs = np.linspace(0, 150e3, 16)
P = [np.mean(np.abs(propagate(np.ones((1, 64), complex), om,
      FiberKernelParams(z, 0, 0, a)))**2) if z else 1.0 for z in zs]
fig, ax = plt.subplots(figsize=(5.5, 3.6))
zz = np.linspace(0, 150e3, 300)
ax.plot(zz/1e3, 10*np.log10(np.exp(-a*zz)), "k-", label=r"$e^{-\alpha z}$ (analítico)")
ax.plot(zs/1e3, 10*np.log10(P), "o", mfc="none", label="SSFM")
ax.set(xlabel="z [km]", ylabel="P(z)/P(0) [dB]", title="1. Atenuação (0,2 dB/km)")
ax.legend(); fig.tight_layout(); fig.savefig(OUT/"01_atenuacao.png"); plt.close(fig)

# ---------------------------------------------------------------- 2. dispersao
N, FS = 1 << 14, 1e12
T = (np.arange(N) - N//2)/FS; OM = 2*np.pi*np.fft.fftfreq(N, 1/FS)
T0 = 20e-12; B2 = beta2_from_dispersion(16e-6, 1550e-9); LD = T0**2/abs(B2)
w0 = fwhm(T, np.abs(gaussian(T, T0))**2)
fig, ax = plt.subplots(figsize=(5.5, 3.6))
xs = np.linspace(0, 2.0, 200)
for C, cor in ((0, "k"), (+2, "tab:blue"), (-2, "tab:red")):
    teo = np.sqrt((1 + C*B2*xs*LD/T0**2)**2 + (B2*xs*LD/T0**2)**2)
    ax.plot(xs, teo, "-", color=cor, label=f"C = {C:+d} (Agrawal)")
    pts = np.linspace(0.1, 2.0, 9)
    sim = [fwhm(T, np.abs(propagate(gaussian(T, T0, C=C)[None, :], OM,
           FiberKernelParams(p*LD, B2, 0, 0), max_step=p*LD/100)[0])**2)/w0 for p in pts]
    ax.plot(pts, sim, "o", color=cor, mfc="none")
ax.axhline(1, color="gray", lw=.8)
ax.set(xlabel=r"$z / L_D$", ylabel=r"FWHM$(z)$ / FWHM$(0)$",
       title="2. Dispersão, pulso gaussiano\n" r"($\beta_2<0$: com C>0 o pulso comprime antes de alargar)")
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(OUT/"02_dispersao.png"); plt.close(fig)

# ---------------------------------------------------------------- 3. SPM
N3, FS3 = 4096, 1e12
T3 = (np.arange(N3) - N3//2)/FS3; OM3 = 2*np.pi*np.fft.fftfreq(N3, 1/FS3)
g, L = 1.3e-3, 50e3; Leff = effective_length(L, 0.2e-3)
fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
P0 = 2.0/(g*Leff)
A = gaussian(T3, T0, P0)[None, :]
o = propagate(A, OM3, FiberKernelParams(L, 0, g, a, manakov=False))[0]
I0 = np.abs(A[0])**2; m = I0 > 1e-6*I0.max()
phi = np.full(N3, np.nan); phi[m] = np.unwrap(np.angle(o[m]*np.conj(A[0][m])))
axs[0].plot(T3*1e12, I0/I0.max(), "k-", label=r"$|A(0,T)|^2$ normalizado")
axs[0].plot(T3*1e12, phi/np.nanmax(phi), "r--", label=r"$\phi_{NL}(L,T)$ normalizado")
axs[0].set(xlim=(-80, 80), xlabel="T [ps]",
           title=rf"3a. Fase NL com a forma do pulso: $\phi_{{max}}$ = {np.nanmax(phi):.3f} rad"
                 "\n" rf"(esperado $\gamma P_0 L_{{eff}}$ = {g*P0*Leff:.3f} rad)")
axs[0].legend(fontsize=8)
f = np.fft.fftshift(np.fft.fftfreq(N3, 1/FS3))
for k, pm in enumerate((0.5*np.pi, 1.5*np.pi, 3.5*np.pi)):
    Pk = pm/(g*Leff)
    Ak = gaussian(T3, T0, Pk)[None, :]
    ok = propagate(Ak, OM3, FiberKernelParams(L, 0, g, a, manakov=False))[0]
    S = np.abs(np.fft.fftshift(np.fft.fft(ok)))**2
    axs[1].plot(f/1e9, S/S.max() + 1.2*k, label=rf"$\phi_{{max}}$ = {pm/np.pi:.1f}$\pi$")
axs[1].set(xlim=(-110, 110), xlabel="Frequência [GHz]", yticks=[],
           title="3b. Alargamento espectral por SPM\n(nº de picos ≈ φmax/π + ½: 1, 2 e 4)")
axs[1].legend(fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig(OUT/"03_nao_linearidade.png"); plt.close(fig)

# ---------------------------------------------------------------- 4. PMD
DF, NS = 1.25e9, 64
OMP = 2*np.pi*np.fft.fftfreq(NS, 1/(NS*DF))
DP, LC = 0.1e-12/np.sqrt(1e3), 50.0
def dgd(L, n, st):
    return jones_dgd(lambda A, r: propagate(A, OMP, FiberKernelParams(L, 0, 0, 0, DP, LC),
                     rng=np.random.default_rng(20000+r), max_step=st), n, DF, NS)
zs = np.array([10, 25, 50, 100, 250, 500, 1e3, 2.5e3, 1e4, 4e4])
rms = [np.sqrt(np.mean(dgd(z, 150, min(z/10, 1e3))**2)) for z in zs]
fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
zz = np.logspace(0.8, 4.7, 300)
db1 = DP/np.sqrt(2*LC)
agr = np.sqrt(2*(db1*LC)**2*(np.exp(-zz/LC) + zz/LC - 1))
axs[0].loglog(zz, agr*1e12, "k-", label="Agrawal (correlação exponencial)")
axs[0].loglog(zz, DP*np.sqrt(zz)*1e12, "k:", label=r"$D_p\sqrt{z}$ (assíntota)")
axs[0].loglog(zs, np.array(rms)*1e12, "o", mfc="none", color="tab:blue",
              label="SSFM, 150 realizações")
axs[0].axvline(LC, color="gray", lw=.8); axs[0].text(LC*1.15, 0.3, r"$l_c$ = 50 m", fontsize=8)
axs[0].set(xlabel="z [m]", ylabel="DGD rms [ps]",
           title="4a. DGD rms vs distância\n(modelos só divergem perto de z ~ l_c)")
axs[0].legend(fontsize=8)
tau = dgd(10e3, 1500, 1e3)
s = np.sqrt(np.mean(tau**2)/3)
x = np.linspace(0, tau.max(), 200)
axs[1].hist(tau*1e12, bins=40, density=True, alpha=.5, label="SSFM, 10 km, 1500 real.")
axs[1].plot(x*1e12, np.sqrt(2/np.pi)*x**2/s**3*np.exp(-x**2/(2*s**2))*1e-12, "k-",
            label="Maxwelliana")
axs[1].set(xlabel="DGD [ps]", ylabel="densidade",
           title="4b. Distribuição do DGD a 10 km\n"
                 rf"($\langle\tau\rangle/\tau_{{rms}}$ = {np.mean(tau)/np.sqrt(np.mean(tau**2)):.3f}; Maxwelliana: 0,921)")
axs[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT/"04_pmd.png"); plt.close(fig)
print("figuras em", OUT.resolve())
