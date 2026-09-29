"""Funcoes de apoio dos testes de fisica. Nao sao testes."""
import numpy as np

def gaussian(t, T0, P0=1.0, C=0.0):
    """Pulso gaussiano, possivelmente com chirp (Agrawal, NLFO, sec. 3.2):
    A(0,T) = sqrt(P0) * exp[-(1 + iC) T^2 / (2 T0^2)]"""
    return np.sqrt(P0) * np.exp(-(1 + 1j * C) * t**2 / (2 * T0**2))

def fwhm(t, intensity):
    """Largura a meia altura, com interpolacao linear nas duas bordas."""
    half = intensity.max() / 2
    above = np.where(intensity >= half)[0]
    i0, i1 = above[0], above[-1]
    tl = np.interp(half, [intensity[i0 - 1], intensity[i0]], [t[i0 - 1], t[i0]])
    tr = np.interp(half, [intensity[i1 + 1], intensity[i1]], [t[i1 + 1], t[i1]])
    return tr - tl

def jones_dgd(propagate_fn, n_real, df=1.25e9, n=64):
    """DGD de primeira ordem pelo metodo de Jones: J(w) e J(w+dw) a partir de
    duas entradas ortogonais; tau = |arg(rho1/rho2)| / dw (BIFROST, sec. IV)."""
    out = []
    for r in range(n_real):
        J = np.empty((n, 2, 2), dtype=complex)
        for k in range(2):
            A = np.zeros((2, n), dtype=complex)
            A[k, 0] = 1.0                     # impulso: espectro plano
            o = propagate_fn(A, r)
            J[:, :, k] = np.fft.fft(o, axis=1).T
        ev = np.linalg.eigvals(np.linalg.solve(J[0], J[1]))
        out.append(abs(np.angle(ev[0] / ev[1])) / (2 * np.pi * df))
    return np.array(out)
