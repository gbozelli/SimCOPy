"""Constantes fisicas e conversoes. Dono unico dos valores que no codigo antigo
apareciam espalhados (lambda0 em 17 lugares, c=3e8 em 4)."""
from __future__ import annotations
import numpy as np
from scipy import constants

C_LIGHT: float = constants.c
H_PLANCK: float = constants.h
K_BOLTZMANN: float = constants.k
Q_ELECTRON: float = constants.e

def db_to_lin(x): return 10.0 ** (x / 10.0)
def lin_to_db(x): return 10.0 * np.log10(x)
def dbm_to_w(x):  return 10.0 ** ((x - 30.0) / 10.0)
def w_to_dbm(x):  return 10.0 * np.log10(x) + 30.0
def db_per_m_to_np_per_m(a): return a / (10.0 / np.log(10.0))

def wavelength_to_frequency(lam): return C_LIGHT / lam
def frequency_to_wavelength(f):   return C_LIGHT / f

def beta2_from_dispersion(D, wavelength):
    """D [s/m^2] -> beta2 [s^2/m]. Convencao beta2 = -D*lambda^2/(2*pi*c)."""
    return -D * wavelength**2 / (2.0 * np.pi * C_LIGHT)

def gamma_from_n2(n2, effective_area, wavelength):
    """gamma [1/(W*m)] = 2*pi*n2/(lambda*Aeff)."""
    return 2.0 * np.pi * n2 / (wavelength * effective_area)

def effective_length(length, attenuation_db_per_m):
    a = db_per_m_to_np_per_m(attenuation_db_per_m)
    return length if a == 0.0 else (1.0 - np.exp(-a * length)) / a
