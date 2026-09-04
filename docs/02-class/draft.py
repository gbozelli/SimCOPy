"""
draft_fixed.py - versao consertada e executavel do simulador coerente DP.

FINALIDADE (Fase 0 do plano de migracao): produzir dados de referencia
reproduzives contra os quais a biblioteca modular sera validada.
NAO e a arquitetura final: continua sendo um script de funcoes soltas.

O QUE FOI CONSERTADO (estrutural - o codigo agora roda ponta a ponta):
  F-1  uma unica definicao por funcao (as duplicatas divergentes foram removidas)
  F-2  nenhuma funcao le variavel global (P_laser, t, M, N_inf, SpS, RollOff, ts
       viraram parametros explicitos)
  F-3  verificacao de limites no pico de correlacao -> SyncError com mensagem,
       em vez de IndexError obscuro
  F-5  nenhuma funcao plota; os dados de plotagem vao para RECORD (ver abaixo)
  F-7  np.complex_ -> np.complex128 (roda em NumPy >= 2.0)
  F-8  guarda de limites no filtro optico
  C-13 saida do fotodiodo com shape (N,) sempre, independente das flags de ruido
       (antes era (1,N) com ruido ligado, e thermal_noise=False quebrava o RX DP)
  bug  DPQAM_receiver ignorava os proprios parametros sync_seed_X/Y e passava
       literais 0 e 123 adiante
  bug  fotodiodos da polarizacao X usavam R=0.9 fixo e BW; os de Y usavam R e
       BW/2. Agora usam os mesmos parametros.

O QUE E CONTROLADO POR LEGACY_PHYSICS:
  True  (padrao) reproduz a fisica do notebook original, bugs inclusos.
                 Use para gerar o baseline de regressao.
  False           aplica C-2 (ASE gaussiano circular), C-5 (potencia somando as
                 duas polarizacoes), C-6 (fator 8/9 de Manakov), C-7
                 (birrefringencia Haar-uniforme com comprimento de correlacao e
                 DGD) e C-9 (shot noise instantaneo).

RESPOSTA HONESTA A "isso simula dados perfeitamente?":
  Nao. Com LEGACY_PHYSICS=True os numeros sao auto-consistentes e reproduziveis,
  mas errados em pontos conhecidos (OSNR, potencia de lancamento, nao
  linearidade, despolarizacao). Servem como referencia de regressao, nao para
  publicacao. Com LEGACY_PHYSICS=False os numeros sao melhores mas ainda faltam
  itens do catalogo de correcoes (C-8 split-step simetrico, C-11 filtro optico
  suave, C-12 RRC + filtro casado).

Uso:
    python draft_fixed.py --scenario b2b
    python draft_fixed.py --scenario span --fixed
    python draft_fixed.py --scenario wdm --channels 4
    python draft_fixed.py --scenario sweep-power
    python draft_fixed.py --scenario sweep-length --out resultados/
"""
from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass, asdict, field

import numpy as np
from scipy import signal

try:                                    # F-5: matplotlib e opcional e headless
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:                     # pragma: no cover
    plt = None

# ----------------------------------------------------------------------------
# Constantes fisicas (F-2/C-15: dono unico, em vez de espalhadas pelo codigo)
# ----------------------------------------------------------------------------
C_LIGHT = 299792458.0          # m/s   (o original usava 3e8: erro de 0,07% em beta2)
LAMBDA0 = 1550e-9              # m
H_PLANCK = 6.62607015e-34      # J.s

N_GUARD_SYNC = 10              # simbolos de guarda do gabarito de sincronismo

LEGACY_PHYSICS = True          # sobrescrito pela CLI


class SyncError(RuntimeError):
    """Sincronismo falhou: o pico de correlacao nao permite extrair o payload."""


# ----------------------------------------------------------------------------
# F-5: em vez de plot_flag espalhado, os blocos gravam o que seria plotado.
# Esta e a versao minima do padrao Probe da secao 5.5 do documento.
# ----------------------------------------------------------------------------
RECORD: dict[str, np.ndarray] = {}
RECORD_ENABLED = True


def _rec(key: str, value) -> None:
    if RECORD_ENABLED:
        RECORD[key] = np.asarray(value)


def save_record(path: str) -> None:
    """Exporta tudo que foi gravado. A plotagem e feita depois, fora daqui."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    np.savez_compressed(path, **RECORD)


# ----------------------------------------------------------------------------
# C-7: birrefringencia / PMD
# ----------------------------------------------------------------------------
@dataclass
class PMDConfig:
    """Modelo de passo grosso, desacoplado do passo numerico DeltaL.

    correlation_length: comprimento de correlacao da birrefringencia [m].
    pmd_coefficient:    coeficiente de PMD [s/sqrt(m)]. 0 desliga o DGD e deixa
                        apenas a rotacao unitaria.
    """
    correlation_length: float = 100.0
    pmd_coefficient: float = 0.1e-12 / np.sqrt(1e3)   # 0,1 ps/sqrt(km)
    total_length: float = 100e3


def random_su2() -> np.ndarray:
    """Matriz de Jones uniformemente distribuida em SU(2) (medida de Haar).

    O codigo original usava apenas 2 angulos; SU(2) tem 3. Com 2 angulos e
    entrada linear em X o estado fica preso no plano S3=0 da esfera de Poincare.
    A distribuicao uniforme em sin^2(theta) e o que faz a medida ser de Haar.
    """
    alpha, beta = np.random.uniform(0.0, 2*np.pi, 2)
    theta = np.arcsin(np.sqrt(np.random.random()))
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c*np.exp(1j*alpha),  s*np.exp(1j*beta)],
                     [-s*np.exp(-1j*beta), c*np.exp(-1j*alpha)]])


def _sections_crossed(auxz: int, DeltaL: float, Lcorr: float) -> int:
    z0, z1 = auxz*DeltaL, (auxz + 1)*DeltaL
    return max(0, int(np.floor(z1/Lcorr)) - int(np.floor(z0/Lcorr))) + (1 if auxz == 0 else 0)


def apply_birefringence(Ax, Ay, f, auxz, DeltaL, pmd: PMDConfig | None):
    """LEGACY: rotacao de 2 angulos a cada passo (fisica dependente de DeltaL).
    NOVO:    rotacao Haar-uniforme + DGD, uma vez por comprimento de correlacao.
    """
    if LEGACY_PHYSICS:
        theta = np.random.rand(1)*np.pi
        phi = np.random.rand(1)*2*np.pi
        ax = Ax*np.exp(1j*phi/2)
        ay = Ay*np.exp(-1j*phi/2)
        return (np.cos(theta)*ax + np.sin(theta)*ay,
                -np.sin(theta)*ax + np.cos(theta)*ay)

    pmd = pmd or PMDConfig()
    n_sec = _sections_crossed(auxz, DeltaL, pmd.correlation_length)
    if n_sec == 0:
        return Ax, Ay

    n_total = max(1, int(round(pmd.total_length/pmd.correlation_length)))
    dtau = pmd.pmd_coefficient*np.sqrt(pmd.total_length)*np.sqrt(np.pi/(2*n_total))

    for _ in range(n_sec):
        if dtau > 0:                                   # DGD de 1a ordem
            Ax_f = np.fft.fftshift(np.fft.fft(Ax))*np.exp(-1j*np.pi*f*dtau)
            Ay_f = np.fft.fftshift(np.fft.fft(Ay))*np.exp(+1j*np.pi*f*dtau)
            Ax = np.fft.ifft(np.fft.ifftshift(Ax_f))
            Ay = np.fft.ifft(np.fft.ifftshift(Ay_f))
        U = random_su2()
        Ax, Ay = U[0, 0]*Ax + U[0, 1]*Ay, U[1, 0]*Ax + U[1, 1]*Ay
    return Ax, Ay


def total_power(E) -> float:
    """C-5: potencia media somando as DUAS polarizacoes.
    O original somava a polarizacao X duas vezes."""
    return float(np.mean(np.abs(E[0, :])**2) + np.mean(np.abs(E[1, :])**2))


# ==========================================================================
# Blocos fisicos e de DSP (portados do notebook original)
# ==========================================================================

def QAM_mod(s,M):

    # Preliminary calculations
    N_bits = int(np.log2(M))
    N_sym = int(len(s)/N_bits)
    s_I = np.zeros(N_sym,dtype=int)
    s_Q = np.zeros(N_sym,dtype=int)

    # N_bits = 1, M = 2
    if N_bits == 1:
        s_1 = np.reshape(s,[N_sym,N_bits])
        for ind_sym in range(N_sym):
            if s_1[ind_sym,:]==0:
                s_I[ind_sym], s_Q[ind_sym] = 0,1
            elif s_1[ind_sym,:]==1:
                s_I[ind_sym], s_Q[ind_sym] = 0,-1

    # N_bits = 2, M = 4
    if N_bits == 2:
        s_1 = np.reshape(s,[N_sym,N_bits])
        for ind_sym in range(N_sym):
            if all(s_1[ind_sym,:]==[0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -1,1
            elif all(s_1[ind_sym,:]==[0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,1
            elif all(s_1[ind_sym,:]==[1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-1
            elif all(s_1[ind_sym,:]==[1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-1

    # N_bits = 3, M = 8
    if N_bits == 3:
        s_1 = np.reshape(s,[N_sym,N_bits])
        for ind_sym in range(N_sym):
            if all(s_1[ind_sym,:]==[0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -3,1
            elif all(s_1[ind_sym,:]==[0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,1
            elif all(s_1[ind_sym,:]==[1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,1
            elif all(s_1[ind_sym,:]==[1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 3,1
            elif all(s_1[ind_sym,:]==[0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,-1
            elif all(s_1[ind_sym,:]==[0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-1
            elif all(s_1[ind_sym,:]==[1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-1
            elif all(s_1[ind_sym,:]==[1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,-1

    # N_bits = 4, M = 16
    if N_bits == 4:
        s_1 = np.reshape(s,[N_sym,N_bits])
        for ind_sym in range(N_sym):
            if all(s_1[ind_sym,:]==[0,0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -3,3
            elif all(s_1[ind_sym,:]==[0,1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,3
            elif all(s_1[ind_sym,:]==[1,1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,3
            elif all(s_1[ind_sym,:]==[1,0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 3,3
            elif all(s_1[ind_sym,:]==[0,0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,1
            elif all(s_1[ind_sym,:]==[0,1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -1,1
            elif all(s_1[ind_sym,:]==[1,1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 1,1
            elif all(s_1[ind_sym,:]==[1,0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,1
            elif all(s_1[ind_sym,:]==[0,0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,-1
            elif all(s_1[ind_sym,:]==[0,1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-1
            elif all(s_1[ind_sym,:]==[1,1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-1,
            elif all(s_1[ind_sym,:]==[1,0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,-1
            elif all(s_1[ind_sym,:]==[0,0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -3,-3
            elif all(s_1[ind_sym,:]==[0,1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-3
            elif all(s_1[ind_sym,:]==[1,1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-3
            elif all(s_1[ind_sym,:]==[1,0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 3,-3

    # N_bits = 5, M = 32
    if N_bits == 5:
        s_1 = np.reshape(s,[N_sym,N_bits])
        for ind_sym in range(N_sym):
            if all(s_1[ind_sym,:]==[1,0,1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -3,5
            elif all(s_1[ind_sym,:]==[1,0,1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,5
            elif all(s_1[ind_sym,:]==[1,1,1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,5
            elif all(s_1[ind_sym,:]==[1,1,1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 3,5
            elif all(s_1[ind_sym,:]==[1,0,1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -5,3
            elif all(s_1[ind_sym,:]==[0,0,1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,3
            elif all(s_1[ind_sym,:]==[0,0,1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,3
            elif all(s_1[ind_sym,:]==[0,1,1,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,3
            elif all(s_1[ind_sym,:]==[0,1,1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,3
            elif all(s_1[ind_sym,:]==[1,1,1,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 5,3
            elif all(s_1[ind_sym,:]==[1,0,1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -5,1
            elif all(s_1[ind_sym,:]==[0,0,1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,1
            elif all(s_1[ind_sym,:]==[0,0,1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,1
            elif all(s_1[ind_sym,:]==[0,1,1,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,1
            elif all(s_1[ind_sym,:]==[0,1,1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,1
            elif all(s_1[ind_sym,:]==[1,1,1,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 5,1
            elif all(s_1[ind_sym,:]==[1,0,0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -5,-1
            elif all(s_1[ind_sym,:]==[0,0,0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,-1
            elif all(s_1[ind_sym,:]==[0,0,0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-1
            elif all(s_1[ind_sym,:]==[0,1,0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-1
            elif all(s_1[ind_sym,:]==[0,1,0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,-1,
            elif all(s_1[ind_sym,:]==[1,1,0,0,1]):
                s_I[ind_sym], s_Q[ind_sym] = 5,-1
            elif all(s_1[ind_sym,:]==[1,0,0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -5,-3
            elif all(s_1[ind_sym,:]==[0,0,0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = -3,-3,
            elif all(s_1[ind_sym,:]==[0,0,0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-3
            elif all(s_1[ind_sym,:]==[0,1,0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-3,
            elif all(s_1[ind_sym,:]==[0,1,0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 3,-3
            elif all(s_1[ind_sym,:]==[1,1,0,1,1]):
                s_I[ind_sym], s_Q[ind_sym] = 5,-3
            elif all(s_1[ind_sym,:]==[1,0,0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = -3,-5
            elif all(s_1[ind_sym,:]==[1,0,0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = -1,-5,
            elif all(s_1[ind_sym,:]==[1,1,0,1,0]):
                s_I[ind_sym], s_Q[ind_sym] = 1,-5
            elif all(s_1[ind_sym,:]==[1,1,0,0,0]):
                s_I[ind_sym], s_Q[ind_sym] = 3,-5
    return s_I+1j*s_Q


def QAM_dem(s,M):
    # Preliminary calculations
    N_bits = int(np.log2(M))
    N_sym = len(s)
    N_bits_total = N_bits*N_sym
    sb = np.zeros(N_bits_total,dtype=int)
    sd = np.zeros(N_sym)

    if N_bits == 2:
        # Find the closest constellation point
        const_points = np.array([-1+1j,1+1j,\
                                -1-1j,1-1j])
        sc_ext = np.tile(s,(M,1))

        const_points_ext = np.transpose(np.tile(const_points,(N_sym,1)))

        #print(np.shape(sc_ext))
        #print(np.shape(sc_ext))
        sym_dist = np.abs(sc_ext-const_points_ext)
        min_dist = np.min(sym_dist,axis=0)
        #print(sc_ext)
        #print(const_points_ext)
        #print(sym_dist)

        ind_array = np.zeros(N_sym)
        for ind_sym in range(N_sym):
            ind = np.where(min_dist[ind_sym]==sym_dist[:,ind_sym])
            if ind[0][0] == 0:
                sb_aux = np.array([0,1])
            elif ind[0][0] == 1:
                sb_aux = np.array([0,0])
            elif ind[0][0] == 2:
                sb_aux = np.array([1,1])
            elif ind[0][0] == 3:
                sb_aux = np.array([1,0])
            ind_i = ind_sym*N_bits
            ind_f = (ind_sym+1)*N_bits
            sb[ind_i:ind_f] = sb_aux
            sd[ind_sym]=min_dist[ind_sym]

    if N_bits == 3:
        # Find the closest constellation point
        const_points = np.array([-3+1j,-1+1j,1+1j,3+1j,\
                                -3-1j,-1-1j,1-1j,3-1j])
        sc_ext = np.tile(s,(M,1))

        const_points_ext = np.transpose(np.tile(const_points,(N_sym,1)))

        #print(np.shape(sc_ext))
        #print(np.shape(sc_ext))
        sym_dist = np.abs(sc_ext-const_points_ext)
        min_dist = np.min(sym_dist,axis=0)
        #print(sc_ext)
        #print(const_points_ext)
        #print(sym_dist)

        ind_array = np.zeros(N_sym)
        for ind_sym in range(N_sym):
            ind = np.where(min_dist[ind_sym]==sym_dist[:,ind_sym])
            if ind[0][0] == 0:
                sb_aux = np.array([0,0,0])
            elif ind[0][0] == 1:
                sb_aux = np.array([0,1,0])
            elif ind[0][0] == 2:
                sb_aux = np.array([1,1,0])
            elif ind[0][0] == 3:
                sb_aux = np.array([1,0,0])
            elif ind[0][0] == 4:
                sb_aux = np.array([0,0,1])
            elif ind[0][0] == 5:
                sb_aux = np.array([0,1,1])
            elif ind[0][0] == 6:
                sb_aux = np.array([1,1,1])
            elif ind[0][0] == 7:
                sb_aux = np.array([1,0,1])

            ind_i = ind_sym*N_bits
            ind_f = (ind_sym+1)*N_bits
            sb[ind_i:ind_f] = sb_aux
            sd[ind_sym]=min_dist[ind_sym]

    if N_bits == 4:
        # Find the closest constellation point
        const_points = np.array([-3+3j,-1+3j,1+3j,3+3j,\
                                -3+1j,-1+1j,1+1j,3+1j,\
                                -3-1j,-1-1j,1-1j,3-1j,\
                                -3-3j,-1-3j,1-3j,3-3j])
        sc_ext = np.tile(s,(M,1))
        const_points_ext = np.transpose(np.tile(const_points,(N_sym,1)))
        sym_dist = np.abs(sc_ext-const_points_ext)
        min_dist = np.min(sym_dist,axis=0)
        #print(sc_ext)
        #print(const_points_ext)
        #print(sym_dist)
        #print(min_dist)
        ind_array = np.zeros(N_sym)
        for ind_sym in range(N_sym):
            #print(min_dist[ind_sym])
            #print(sym_dist[:,ind_sym])
            ind = np.where(min_dist[ind_sym]==sym_dist[:,ind_sym])
            #print(ind[0][0])
            #ind_array[ind_sym] = ind[0][0]
            if ind[0][0] == 0:
                sb_aux = np.array([0,0,1,0])
            elif ind[0][0] == 1:
                sb_aux = np.array([0,1,1,0])
            elif ind[0][0] == 2:
                sb_aux = np.array([1,1,1,0])
            elif ind[0][0] == 3:
                sb_aux = np.array([1,0,1,0])
            elif ind[0][0] == 4:
                sb_aux = np.array([0,0,1,1])
            elif ind[0][0] == 5:
                sb_aux = np.array([0,1,1,1])
            elif ind[0][0] == 6:
                sb_aux = np.array([1,1,1,1])
            elif ind[0][0] == 7:
                sb_aux = np.array([1,0,1,1])
            elif ind[0][0] == 8:
                sb_aux = np.array([0,0,0,1])
            elif ind[0][0] == 9:
                sb_aux = np.array([0,1,0,1])
            elif ind[0][0] == 10:
                sb_aux = np.array([1,1,0,1])
            elif ind[0][0] == 11:
                sb_aux = np.array([1,0,0,1])
            elif ind[0][0] == 12:
                sb_aux = np.array([0,0,0,0])
            elif ind[0][0] == 13:
                sb_aux = np.array([0,1,0,0])
            elif ind[0][0] == 14:
                sb_aux = np.array([1,1,0,0])
            elif ind[0][0] == 15:
                sb_aux = np.array([1,0,0,0])

            ind_i = ind_sym*N_bits
            ind_f = (ind_sym+1)*N_bits
            sb[ind_i:ind_f] = sb_aux
            sd[ind_sym]=min_dist[ind_sym]

    if N_bits == 5:
        # Find the closest constellation point
        const_points = np.array([-3+5j,-1+5j,1+5j,3+5j,\
                                -5+3j,-3+3j,-1+3j,1+3j,3+3j,5+3j,\
                                -5+1j,-3+1j,-1+1j,1+1j,3+1j,5+1j,\
                                -5-1j,-3-1j,-1-1j,1-1j,3-1j,5-1j,\
                                -5-3j,-3-3j,-1-3j,1-3j,3-3j,5-3j,\
                                -3-5j,-1-5j,1-5j,3-5j])
        sc_ext = np.tile(s,(M,1))
        const_points_ext = np.transpose(np.tile(const_points,(N_sym,1)))
        sym_dist = np.abs(sc_ext-const_points_ext)
        min_dist = np.min(sym_dist,axis=0)
        #print(sc_ext)
        #print(const_points_ext)
        #print(sym_dist)
        #print(min_dist)
        ind_array = np.zeros(N_sym)
        for ind_sym in range(N_sym):
            #print(min_dist[ind_sym])
            #print(sym_dist[:,ind_sym])
            ind = np.where(min_dist[ind_sym]==sym_dist[:,ind_sym])
            #print(ind[0][0])
            #ind_array[ind_sym] = ind[0][0]
            if ind[0][0] == 0:
                sb_aux = np.array([1,0,1,0,0])
            elif ind[0][0] == 1:
                sb_aux = np.array([1,0,1,1,0])
            elif ind[0][0] == 2:
                sb_aux = np.array([1,1,1,1,0])
            elif ind[0][0] == 3:
                sb_aux = np.array([1,1,1,0,0])
            elif ind[0][0] == 4:
                sb_aux = np.array([1,0,1,1,1])
            elif ind[0][0] == 5:
                sb_aux = np.array([0,0,1,1,1])
            elif ind[0][0] == 6:
                sb_aux = np.array([0,0,1,1,0])
            elif ind[0][0] == 7:
                sb_aux = np.array([0,1,1,1,0])
            elif ind[0][0] == 8:
                sb_aux = np.array([0,1,1,1,1])
            elif ind[0][0] == 9:
                sb_aux = np.array([1,1,1,1,1])
            elif ind[0][0] == 10:
                sb_aux = np.array([1,0,1,0,1])
            elif ind[0][0] == 11:
                sb_aux = np.array([0,0,1,0,1])
            elif ind[0][0] == 12:
                sb_aux = np.array([0,0,1,0,0])
            elif ind[0][0] == 13:
                sb_aux = np.array([0,1,1,0,0])
            elif ind[0][0] == 14:
                sb_aux = np.array([0,1,1,0,1])
            elif ind[0][0] == 15:
                sb_aux = np.array([1,1,1,0,1])
            elif ind[0][0] == 16:
                sb_aux = np.array([1,0,0,0,1])
            elif ind[0][0] == 17:
                sb_aux = np.array([0,0,0,0,1])
            elif ind[0][0] == 18:
                sb_aux = np.array([0,0,0,0,0])
            elif ind[0][0] == 19:
                sb_aux = np.array([0,1,0,0,0])
            elif ind[0][0] == 20:
                sb_aux = np.array([0,1,0,0,1])
            elif ind[0][0] == 21:
                sb_aux = np.array([1,1,0,0,1])
            elif ind[0][0] == 22:
                sb_aux = np.array([1,0,0,1,1])
            elif ind[0][0] == 23:
                sb_aux = np.array([0,0,0,1,1])
            elif ind[0][0] == 24:
                sb_aux = np.array([0,0,0,1,0])
            elif ind[0][0] == 25:
                sb_aux = np.array([0,1,0,1,0])
            elif ind[0][0] == 26:
                sb_aux = np.array([0,1,0,1,1])
            elif ind[0][0] == 27:
                sb_aux = np.array([1,1,0,1,1])
            elif ind[0][0] == 28:
                sb_aux = np.array([1,0,0,0,0])
            elif ind[0][0] == 29:
                sb_aux = np.array([1,0,0,1,0])
            elif ind[0][0] == 30:
                sb_aux = np.array([1,1,0,1,0])
            elif ind[0][0] == 31:
                sb_aux = np.array([1,1,0,0,0])

            ind_i = ind_sym*N_bits
            ind_f = (ind_sym+1)*N_bits
            sb[ind_i:ind_f] = sb_aux
            sd[ind_sym]=min_dist[ind_sym]

    return [sb,sd]


def DAC_Nyquist(s=None,SpS=16,RollOff=0.2,plot_flag=False):

    # Initial calculations
    N_sym = len(s)
    N_samples = N_sym*SpS

    # Conversion to the frequency domain
    S = np.fft.fftshift(np.fft.fft(s))
    S_os = np.zeros(N_samples,dtype=complex)

    # Ideal interpolation (Nyquist limit)
    ind_i = int(N_samples/2-N_sym/2)
    ind_f = ind_i+N_sym
    S_os[ind_i:ind_f] = S
    s_os = SpS*np.fft.ifft(np.fft.ifftshift(S_os))
    S_nyq = np.zeros(N_samples,dtype=complex)
    ind_i = int(N_samples/2-N_sym-N_sym/2)
    ind_f = ind_i+N_sym
    S_nyq[ind_i:ind_f] = S
    ind_i = int(N_samples/2-N_sym/2)
    ind_f = ind_i+N_sym
    S_nyq[ind_i:ind_f] = S
    ind_i = int(N_samples/2+N_sym-N_sym/2)
    ind_f = ind_i+N_sym
    S_nyq[ind_i:ind_f] = S

    ind = np.arange(N_samples)
    N_aux = int(RollOff*N_sym)
    ind_aux = np.arange(int(N_aux))
    s_aux_1 = 0.5*(1-np.cos(np.pi*ind_aux/N_aux))
    s_aux_2 = 0.5*(1+np.cos(np.pi*ind_aux/N_aux))
    window = np.zeros(N_samples)
    ind_i_1 = int(N_samples/2-N_sym/2-N_aux/2)
    ind_f_1 = ind_i_1+int(N_aux)
    window[ind_i_1:ind_f_1] = s_aux_1
    ind_i_2 = int(N_samples/2+N_sym/2-N_aux/2)
    ind_f_2 = ind_i_2+int(N_aux)
    window[ind_i_2:ind_f_2] = s_aux_2
    window[ind_f_1:ind_i_2] = np.ones(ind_i_2-ind_f_1)

    if plot_flag:
        plt.figure()
        plt.plot(np.abs(S_nyq)/max(np.abs(S_nyq)))
        plt.plot(np.abs(window)/max(window))

    S_nyq = window*S_nyq
    s_nyq = SpS*np.fft.ifft(np.fft.ifftshift(S_nyq))

    return s_nyq


def QAM_transmitter(N_sync=128,sync_seed=0,\
                    N_MIMO=32,MIMO_i=0,\
                    N_inf=1024,M=16,SpS=16,RollOff=0.2,ts=1e-12,
                    N_zeros_init=20,N_zeros_final=20,plot_flag=True):

    ## Sequence generation

    # Initial zeros
    s_zeros_init = np.zeros(N_zeros_init)
    # Syncronization sequence
    np.random.seed(sync_seed)
    s_sync_i = np.random.randint(0,2,N_sync,dtype=int)*2-1
    # MIMO sequence
    #print(MIMO_i)
    if MIMO_i == 1:
        s_MIMO_i = np.append(np.zeros(int(N_MIMO/2)),np.ones(int(N_MIMO/2)))
    else:
        s_MIMO_i = np.append(np.ones(int(N_MIMO/2)),np.zeros(int(N_MIMO/2)))
    s_MIMO_i = s_MIMO_i*np.mod(np.arange(N_MIMO)+MIMO_i,2)
    #print(s_MIMO_i)

    # Information sequence
    N_bits = int(np.log2(M))
    N_bits_total = int(N_bits*N_inf)
    s_b = np.random.randint(0,2,N_bits_total,dtype=int)
    s_inf = QAM_mod(s_b,M)
    s_inf = s_inf/np.max(np.real(s_inf))

    # Final zeros
    s_zeros_final = np.zeros(N_zeros_final)

    # Concatenate
    s = np.concatenate([s_zeros_init,s_sync_i,s_inf,s_zeros_final])

    ## Conversion to the continuous time
    s_cont = DAC_Nyquist(s=s,SpS=SpS,RollOff=RollOff)

    ind = np.arange(len(s))*SpS

    t = ts*(np.arange(len(s_cont)))

    if plot_flag:
        plt.figure()
        plt.plot(t*1e9,np.real(s_cont))
        plt.plot(ind*ts*1e9,np.real(s),'.')
        plt.xlabel('Time [ns]')
        plt.ylabel('In-phase component [a.u.]')
        plt.xlim([0,np.max(t)*1e9])

        plt.figure()
        plt.plot(np.real(s_cont),np.imag(s_cont),alpha=0.2)
        plt.plot(np.real(s),np.imag(s),'.')
        plt.xlabel('In-phase')
        plt.ylabel('Quadrature')
        plt.axis('square')

    return [s_cont,t,s_b]


def single_frequency_laser(t = None, P = 0.01, Delta_nu = 1e6, Freq_offset = 0, X_fraction = 1, phi_pol = 0):
    # Calculate the amplitude
    A = np.sqrt(P)
    # Calculate the phase noise
    ts = t[1]-t[0]
    freq_noise = np.sqrt(2*np.pi*Delta_nu*ts)*np.random.randn(len(t))
    phase_noise = np.cumsum(freq_noise)
    # Frequency offset
    phase = 2*np.pi*(Freq_offset)*t+phase_noise
    # Generation of the signal
    s = A*np.exp(1j*phase)
    Ax = np.sqrt(X_fraction)*np.exp(1j*phi_pol/2)*s
    Ay = np.sqrt(1-X_fraction)*np.exp(-1j*phi_pol/2)*s
    A = np.zeros((2,len(t)),dtype=complex)
    #print(np.shape(A))
    A[0,:] = Ax
    A[1,:] = Ay
    return(A)


def MZM(A,Vrf_upper,Vrf_lower,Splitter_IL=1.0,Upper_IL=2.0,Lower_IL=3.0,Combiner_IL=1.1,Vpi = 5,Vbias=2.5):
    Aupper_1 = A/2*10**(-Splitter_IL/20)
    Alower_1 = A/2*10**(-Splitter_IL/20)
    Aupper_2 = Aupper_1*np.exp(1j*(-Vrf_upper+Vbias/2)/Vpi*np.pi)
    Alower_2 = Alower_1*np.exp(1j*(-Vrf_lower-Vbias/2)/Vpi*np.pi)
    Aupper_3 = Aupper_2*10**(-Upper_IL/20)
    Alower_3 = Alower_2*10**(-Lower_IL/20)
    Aupper_4 = Aupper_3*10**(-Combiner_IL/20)
    Alower_4 = Alower_3*10**(-Combiner_IL/20)
    Aout = Aupper_4+Alower_4
    return Aout


def HybridNetwork(A_sig,A_LO):
    A_1 = 1/2*(A_sig+A_LO)
    A_2 = 1/2*(A_sig-A_LO)
    A_3 = 1/2*(A_sig+1j*A_LO)
    A_4 = 1/2*(A_sig-1j*A_LO)
    return [A_1,A_2,A_3,A_4]


def photodiode(s=None,t=None,R=0.9,BW=10e9,shot_noise=True,thermal_noise=True):
    ts = t[1]-t[0]                                               # era t[2]-t[1]
    fs = 1/ts
    BW_sim = 1/ts
    # Extract each polarization
    N = len(t)
    Ax = s[:]

    #print('The length of Ax is',np.shape(Ax))
    #print('The length of t is',np.shape(t))

    # Photodetection
    #I_opt = (np.abs(Ax))**2
    I_opt = np.multiply(abs(Ax),abs(Ax))
    i_pd = R*I_opt

    # Noise addition
    kB = 1.38e-23       # Boltzmann's constant in J/K
    Temp = 300          # Temperature in K
    Rl = 50             # Load resistance in Ohm
    q = 1.6e-19         # Electron charge in C
    if thermal_noise:
        # Thermal noise
        PSD_n_thermal = 4*kB*Temp/Rl
        P_n_thermal = PSD_n_thermal*BW_sim
        n_thermal = np.sqrt(P_n_thermal)*np.random.randn(N)          # C-13: era (1,N)
    else:
        n_thermal = 0

    if shot_noise:
        # Shot noise
        # C-9: shot noise dependente do sinal instantaneo quando LEGACY_PHYSICS=False
        PSD_n_shot = 2*q*(np.mean(np.abs(i_pd)) if LEGACY_PHYSICS else np.abs(i_pd))
        P_n_shot = PSD_n_shot*BW_sim
        n_shot = np.sqrt(P_n_shot)*np.random.randn(N)
    else:
        n_shot = 0
    i_pd = i_pd + n_thermal+n_shot

    # Filtering
    wp = BW
    ws = 1.5*wp      # Begin of the stop band
    gpass = 3        # Maximum loss in the pass band in dB
    gstop = 20       # Minimum loss in the stop band in dB
    analog = False   #

    # Finding the order
    No, Wn = signal.ellipord(wp, ws, gpass, gstop, analog, fs)

    # Create the filter
    b, a = signal.ellip(No, gpass, gstop, Wn, btype='low', analog=False, fs=fs)

    i_pd = signal.lfilter(b,a,i_pd)

    return i_pd


def fiber(A=0,t=0,L=100e3,DeltaL=1e3,D=16e-6,lambda0=1550e-9,alpha_dB=0.2e-3,n2=2.6e-20,Aeff=80e-12,SoP_rotation=None,pmd=None):
    # E: input electrical double-polarization field in V/m
    # t: time vector in [s]
    # L: fiberl length in [m]
    # DeltaL: the longitudinal step value in [m]
    # D: dispersion parameter in [s/m^2]
    # lambda0: operation wavelength
    # alpha_dB: attenuation coefficient in dB/m
    # n2: nonlinear refractive index in [W-1]
    # Aeff: effective modal area in [mˆ2]

    # Constants
    c = C_LIGHT         # C-15: constante unica
    pi = np.pi          # Pi number

    # Preliminary calculations
    Nt = len(t)                    # Number of time samples
    Deltat = t[1]-t[0]             # Time step in s
    fs = 1/Deltat                  # Sampling frequency in Hz
    f = np.arange(Nt)/(Nt/fs)-fs/2 # Frequency vector in Hz
    Nz = int(L/DeltaL)             # Number of segments
    z = np.linspace(0,L-DeltaL,Nz) # Propagation distance vector in m
    alpha_Np = alpha_dB/4.343      # Attenuation of in Np/m
    k = 2*pi/lambda0               # Wave number in m-1
    gamma = n2*k/Aeff              # Nonlinear coefficient in W-1m-1

    # Extraction of the polarizations
    Ax = A[0,:]                    # We extract the x polarization
    Ay = A[1,:]                    # We extract the y polarization

    # Propagate through the fiber
    for auxz in range(len(z)):

        # Scalar linear effects are simulated in the frequency domain

        # Conversion to frequency domain
        Ax_f = np.fft.fftshift(np.fft.fft(Ax))
        Ay_f = np.fft.fftshift(np.fft.fft(Ay))

        # Attenuation
        Ax_f = Ax_f*np.exp(-alpha_Np/2*DeltaL)
        Ay_f = Ay_f*np.exp(-alpha_Np/2*DeltaL)

        # Dispersion
        Ax_f *= np.exp(-1j*(D*lambda0**2/(4*np.pi*c)*DeltaL*((2*np.pi)*f)**2))
        Ay_f *= np.exp(-1j*(D*lambda0**2/(4*np.pi*c)*DeltaL*((2*np.pi)*f)**2))

        # Conversion to time domain
        Ax = np.fft.ifft(np.fft.ifftshift(Ax_f))
        Ay = np.fft.ifft(np.fft.ifftshift(Ay_f))

        # Nonlinear effects are simulated in the time domain
        # C-6: fator 8/9 de Manakov para propagacao com duas polarizacoes
        Ixy = ((abs(Ax))**2+(abs(Ay))**2)*(1.0 if LEGACY_PHYSICS else 8.0/9.0)
        Ax *= np.exp(1j*gamma*Ixy*DeltaL);
        Ay *= np.exp(1j*gamma*Ixy*DeltaL);

        # Polarization rotation can be performed either in time or frequency domain
        # C-7: birrefringencia / PMD  (ver secao 10 do documento de arquitetura)
        if SoP_rotation:
            Ax, Ay = apply_birefringence(Ax, Ay, f, auxz, DeltaL, pmd)

    # Create the output field from the two polarizations
    A = np.zeros((2,len(t)),dtype = np.complex128)
    A[0,:] = Ax
    A[1,:] = Ay
    return A


def EDFA(s_in,t,G_edfa_dB,NF_dB,lambda0=1550e-9):
    h = 6.626e-34     # Planck's constant in J/Hz
    c = C_LIGHT
    f0 = c/lambda0    # Operation frequency in Hz
    ts = t[1]-t[0]    # Sampling time in s
    fs = 1/ts         # Sampling frequency in Hz
    BW_sim = fs       # Simulation bandwidth in Hz

    # Gain
    G_edfa = 10**(G_edfa_dB/10)
    sx = s_in[0]
    sy = s_in[1]
    N = len(sx)

    sx = np.sqrt(G_edfa)*sx
    sy = np.sqrt(G_edfa)*sy

    # Noise
    F_n = 10**(NF_dB/10)
    n_sp = F_n/2
    S_ASE = n_sp*h*f0*(G_edfa-1)
    P_ASE = S_ASE*BW_sim

    # C-2: ASE deve ser gaussiano circular complexo, nao uniforme
    if LEGACY_PHYSICS:
        s_ASEx = np.sqrt(P_ASE)/2*np.random.rand(N)+1j*np.sqrt(P_ASE)/2*np.random.rand(N)
        s_ASEy = np.sqrt(P_ASE)/2*np.random.rand(N)+1j*np.sqrt(P_ASE)/2*np.random.rand(N)
    else:
        s_ASEx = np.sqrt(P_ASE/2)*(np.random.randn(N)+1j*np.random.randn(N))
        s_ASEy = np.sqrt(P_ASE/2)*(np.random.randn(N)+1j*np.random.randn(N))
    s_outx = sx + s_ASEx
    s_outy = sy + s_ASEy

    s = np.zeros((2,len(t)),dtype = np.complex128)
    s[0,:] = s_outx
    s[1,:] = s_outy
    return s


def optical_filter(s,t,f0,BW):
    ts = t[1]-t[0]                       # Sampling period in s
    fs = 1/ts                            # Sampling frequancy in Hz
    Deltaf = 1/(ts*len(t))               # Frequency resolution in Hz
    f = np.arange(len(t))*Deltaf-fs/2    # Frequency vector in Hz
    ind1 = np.where(f<(f0-BW/2))
    ind2 = np.where(f<(f0+BW/2))
    if len(ind1[0]) == 0 or len(ind2[0]) == 0:                   # F-8: guarda de limites
        raise ValueError(f"filtro optico fora da banda de simulacao: f0={f0/1e9:.1f} GHz, "
                         f"BW={BW/1e9:.1f} GHz, banda disponivel=+-{fs/2e9:.1f} GHz")
    #print(ind1[0][-1])
    #print(ind2[0][-1])
    #print(f[ind1[0][-1]]/1e9)
    #print(f[ind2[0][-1]]/1e9)

    S = np.fft.fftshift(np.fft.fft(s))

    H = np.zeros(len(f))
    H[ind1[0][-1]:ind2[0][-1]]=1
    Sf = S*H

    #plt.figure()
    #plt.plot(f,np.abs(S[0,:]))
    #plt.plot(f,np.abs(Sf[0,:]))

    sf = np.fft.ifft(np.fft.ifftshift(Sf))

    return sf


def DPQAM_transmitter(P_laser = 0.01, M = 16, SpS = 16, RollOff = 0.2, ts = 1e-9, sync_seed_X=0,sync_seed_Y=123,\
                           N_MIMO = 128,N_sync = 128,N_inf = 4096,N_zeros_init=10,N_zeros_final=10,Delta_nu=0,Freq_offset=0,\
                           ind_mod = 0.1, Splitter_IL=1.0,Upper_IL=2.0,Lower_IL=2.0,Combiner_IL=1.1,\
                           Vpi = 2.5, Vbias=2.5, plot_flag = False):

    #------------------------------------- Electrical signal --------------------------------------------%
    # Complex signal for X polarization
    s_tx_X, t, s_b_X = QAM_transmitter(M = M, SpS = SpS, RollOff = RollOff,ts = ts,sync_seed=sync_seed_X,\
                           N_sync = N_sync,N_MIMO=N_MIMO,MIMO_i=0,N_inf = N_inf,N_zeros_init=N_zeros_init,N_zeros_final=N_zeros_final,\
                           plot_flag = False)

    # Complex signal for Y polarization
    s_tx_Y, t, s_b_Y = QAM_transmitter(M = M, SpS = SpS, RollOff = RollOff,ts = ts,sync_seed=sync_seed_Y,\
                           N_sync = N_sync,N_MIMO=N_MIMO,MIMO_i=1,N_inf = N_inf,N_zeros_init=N_zeros_init,N_zeros_final=N_zeros_final,\
                           plot_flag = False)

    # Combination of the signals for both polarizations
    # Binary sequences
    s_b = np.zeros((2,len(s_b_X)),dtype=int)
    s_b[0,:] = s_b_X
    s_b[1,:] = s_b_Y
    # Electrical signals
    s_tx = np.zeros((2,len(t)),dtype=complex)
    s_tx[0,:] = s_tx_X
    s_tx[1,:] = s_tx_Y

    #--------------------------------------- Optical signal ---------------------------------------------%
    E_carrier = single_frequency_laser(t = t, P = P_laser,      # F-2: era lida do escopo global Delta_nu = Delta_nu,\
                                       Freq_offset = Freq_offset,X_fraction = 0.5, phi_pol = 0)

    # Polarization splitting
    E_carrier_X = E_carrier[0,:]
    E_carrier_Y = E_carrier[1,:]

    #------------------------------------- Optical modulation --------------------------------------------%
    # Modulation depth
    Kmod = ind_mod

    # Modulation of the X component
    # Modulation of the in-phase component of the X polarization
    E_carrier_X_I = E_carrier_X/np.sqrt(2)
    s_tx_X_I = np.real(s_tx_X)
    Vrf_upper = Kmod*s_tx_X_I
    Vrf_lower = -Kmod*s_tx_X_I
    E_mod_X_I = MZM(E_carrier_X_I,Vrf_upper,Vrf_lower,Splitter_IL=Splitter_IL,Upper_IL=Upper_IL,\
                    Lower_IL=Lower_IL,Combiner_IL=Combiner_IL,Vpi = Vpi,Vbias=Vbias)
    # Modulation of the in-phase component of the X polarization
    E_carrier_X_Q = E_carrier_X/np.sqrt(2)*np.exp(1j*np.pi/2)
    s_tx_X_Q = np.imag(s_tx_X)
    Vrf_upper = Kmod*s_tx_X_Q
    Vrf_lower = -Kmod*s_tx_X_Q
    E_mod_X_Q = MZM(E_carrier_X_Q,Vrf_upper,Vrf_lower,Splitter_IL=Splitter_IL,Upper_IL=Upper_IL,\
                    Lower_IL=Lower_IL,Combiner_IL=Combiner_IL,Vpi = Vpi,Vbias=Vbias)
    # Combine the components
    E_mod_X = 1/np.sqrt(2)*E_mod_X_I+1/np.sqrt(2)*E_mod_X_Q

    # Modulation of the Y component
    # Modulation of the in-phase component of the Y polarization
    E_carrier_Y_I = E_carrier_Y/np.sqrt(2)
    s_tx_Y_I = np.real(s_tx_Y)
    Vrf_upper = Kmod*s_tx_Y_I
    Vrf_lower = -Kmod*s_tx_Y_I
    E_mod_Y_I = MZM(E_carrier_Y_I,Vrf_upper,Vrf_lower,Splitter_IL=Splitter_IL,Upper_IL=Upper_IL,\
                    Lower_IL=Lower_IL,Combiner_IL=Combiner_IL,Vpi = Vpi,Vbias=Vbias)
    # Modulation of the in-phase component of the X polarization
    E_carrier_Y_Q = E_carrier_Y/np.sqrt(2)*np.exp(1j*np.pi/2)
    s_tx_Y_Q = np.imag(s_tx_Y)
    Vrf_upper = Kmod*s_tx_Y_Q
    Vrf_lower = -Kmod*s_tx_Y_Q
    E_mod_Y_Q = MZM(E_carrier_Y_Q,Vrf_upper,Vrf_lower,Splitter_IL=Splitter_IL,Upper_IL=Upper_IL,\
                    Lower_IL=Lower_IL,Combiner_IL=Combiner_IL,Vpi = Vpi,Vbias=Vbias)
    # Combine the components
    E_mod_Y = 1/np.sqrt(2)*E_mod_Y_I+1/np.sqrt(2)*E_mod_Y_Q

    if plot_flag:
        plt.figure()
        plt.subplot(2,1,1)
        plt.plot(np.real(E_mod_X))
        plt.plot(np.imag(E_mod_X))
        plt.subplot(2,1,2)
        plt.plot(np.real(E_mod_Y))
        plt.plot(np.imag(E_mod_Y))

    # Combine polarizations
    E_opt = np.zeros((2,len(t)),dtype=complex)
    #print(np.shape(A))
    E_opt[0,:] = E_mod_X
    E_opt[1,:] = E_mod_Y

    return [t, E_opt, E_carrier, s_tx, s_b]


def QAM_receiver_DP(s,t,s_b,M,N_sync=128,N_MIMO=32,N_inf=4096,SpS=16,RollOff=0.2,sync_seed_X=0,sync_seed_Y=123,ts=1e-12,plot_flag=False,L=10e3,D=16e-6):
    #print('The length is',L)
    #print('The sipersion parameters is',D)
    s_in = s
    s_x = s[0,:]
    s_y = s[1,:]

    # Generation of time vector
    N_samples = len(s_x)
    t = ts*(np.arange(N_samples))
    fs = 1/ts
    f = (np.arange(N_samples)/N_samples)*fs-fs/2

    # Coarse frequency offset compensation
    N_samples = len(s_x)
    N_sym = int(N_samples/SpS)
    S_nyq = np.zeros(N_samples,dtype=complex)
    ind_i = int(N_samples/2-N_sym/2)
    ind_f = int(ind_i+N_sym)
    window = np.zeros(N_samples)
    window[ind_i:ind_f] = np.ones(ind_f-ind_i)

    f_min, f_max, N_f = -500e6,500e6,101
    f_offset_array = np.linspace(f_min,f_max,N_f)
    P_array_x = np.linspace(f_min,f_max,N_f)
    P_array_y = np.linspace(f_min,f_max,N_f)
    for counter, offset in enumerate(f_offset_array):
        s_shift_x = np.exp(1j*2*np.pi*offset*t)*s_x
        S_shift_x = np.abs(np.fft.fftshift(np.fft.fft(s_shift_x)))
        S_f_x = S_shift_x*window
        P_array_x[counter] = sum(np.abs(S_f_x)**2)

        s_shift_y = np.exp(1j*2*np.pi*offset*t)*s_y
        S_shift_y = np.abs(np.fft.fftshift(np.fft.fft(s_shift_y)))
        S_f_y = S_shift_y*window
        P_array_y[counter] = sum(np.abs(S_f_y)**2)
    P_array = (P_array_x+P_array_y)/2

    Npolord = 10
    max_offset_opt_x = max(P_array_x)
    ind_offset_opt_x = np.where(P_array_x==max_offset_opt_x)
    f_offset_opt_x = f_offset_array[ind_offset_opt_x[0][0]]
    p_x = np.polyfit(f_offset_array,P_array_x,Npolord)
    max_offset_opt_y = max(P_array_y)
    ind_offset_opt_y = np.where(P_array_y==max_offset_opt_y)
    f_offset_opt_y = f_offset_array[ind_offset_opt_y[0][0]]
    p_y = np.polyfit(f_offset_array,P_array_y,Npolord)
    max_offset_opt = max(P_array)
    ind_offset_opt = np.where(P_array==max_offset_opt)
    f_offset_opt = f_offset_array[ind_offset_opt[0][0]]
    p = np.polyfit(f_offset_array,P_array,Npolord)

    N_f = 2001
    f_offset_array_ext = np.linspace(f_min,f_max,N_f)
    P_array_ext_x = np.polyval(p_x,f_offset_array_ext)
    P_array_ext_y = np.polyval(p_y,f_offset_array_ext)
    P_array_ext = np.polyval(p,f_offset_array_ext)

    max_offset_opt_x = max(P_array_ext_x)
    ind_offset_opt_x = np.where(P_array_ext_x==max_offset_opt_x)
    f_offset_opt_x = f_offset_array_ext[ind_offset_opt_x[0][0]]

    max_offset_opt_y = max(P_array_ext_y)
    ind_offset_opt_y = np.where(P_array_ext_y==max_offset_opt_y)
    f_offset_opt_y = f_offset_array_ext[ind_offset_opt_y[0][0]]

    max_offset_opt = max(P_array_ext)
    ind_offset_opt = np.where(P_array_ext==max_offset_opt)
    f_offset_opt = f_offset_array_ext[ind_offset_opt[0][0]]

    if plot_flag:
        plt.figure()
        plt.plot(f_offset_array/1e6,P_array_x,'o')
        plt.plot(f_offset_array/1e6,P_array_y,'o')
        plt.plot(f_offset_array/1e6,P_array,'o')
        plt.plot(f_offset_array_ext/1e6,P_array_ext_x)
        plt.plot(f_offset_array_ext/1e6,P_array_ext_y)
        plt.plot(f_offset_array_ext/1e6,P_array_ext)
        plt.plot(f_offset_array_ext[ind_offset_opt_x[0][0]]/1e6,max_offset_opt_x,'^r')
        plt.plot(f_offset_array_ext[ind_offset_opt_y[0][0]]/1e6,max_offset_opt_y,'^r')
        plt.plot(f_offset_array_ext[ind_offset_opt[0][0]]/1e6,max_offset_opt,'^r')
        plt.ylabel('Power [a.u.]')
        plt.xlabel('Frequency offset [MHz]')

    #print('The frequency offset for x is:',f_offset_opt_x/1e6)
    #print('The frequency offset for y is:',f_offset_opt_y/1e6)
    #print('The frequency offset for x-y is:',f_offset_opt/1e6)
    s = np.exp(1j*2*np.pi*f_offset_opt*t)*s
    #s = s_in

    # Filter the signal using a square filter
    N_samples = len(s[0,:])
    N_sym = int(N_samples/SpS)
    S_nyq = np.zeros(N_samples,dtype=complex)
    ind_i = int(N_samples/2-N_sym-N_sym/2)
    ind_f = ind_i+N_sym

    ind = np.arange(N_samples)
    N_aux = int(RollOff*N_sym)
    ind_aux = np.arange(int(N_aux))
    s_aux_1 = 1*np.ones(N_aux)   #0.5*(1-np.cos(np.pi*ind_aux/N_aux))
    s_aux_2 = 1*np.ones(N_aux)   #0.5*(1+np.cos(np.pi*ind_aux/N_aux))
    window = np.zeros(N_samples)
    ind_i_1 = int(N_samples/2-N_sym/2-N_aux/2)
    ind_f_1 = ind_i_1+int(N_aux)
    window[ind_i_1:ind_f_1] = s_aux_1
    ind_i_2 = int(N_samples/2+N_sym/2-N_aux/2)
    ind_f_2 = ind_i_2+int(N_aux)
    window[ind_i_2:ind_f_2] = s_aux_2
    window[ind_f_1:ind_i_2] = np.ones(ind_i_2-ind_f_1)

    sx = s[0,:]
    Sx = np.fft.fftshift(np.fft.fft(sx))
    S_f_x = Sx*window
    s_f_x = SpS*np.fft.ifft(np.fft.ifftshift(S_f_x))
    sy = s[1,:]
    Sy = np.fft.fftshift(np.fft.fft(sy))
    S_f_y = Sy*window
    s_f_y = SpS*np.fft.ifft(np.fft.ifftshift(S_f_y))

    s[0,:] = s_f_x
    s[1,:] = s_f_y
    #s = s_in

    #s = s_in

    if plot_flag:
        epsilon = 1e-12
        plt.figure()
        plt.plot(f/1e9,20*np.log10(np.abs(Sx)/max(np.abs(Sx))+epsilon),label='signal')
        plt.plot(f/1e9,20*np.log10(np.abs(Sy)/max(np.abs(Sy))+epsilon),label='signal')
        plt.plot(f/1e9,20*np.log10(np.abs(S_f_x)/max(np.abs(S_f_x))+epsilon),label='signal')
        plt.plot(f/1e9,20*np.log10(np.abs(S_f_y)/max(np.abs(S_f_y))+epsilon),label='signal')
        plt.plot(f/1e9,20*np.log10(window/max(np.abs(window))+epsilon),label='filter')
        plt.ylabel('Normalized power / Filter resp. [dB]')
        plt.xlabel('Frequency [GHz]')

    ####### Dispersion compensation #########

    # Conversion to frequency domain
    s_f_x_f = np.fft.fftshift(np.fft.fft(s_f_x))
    s_f_y_f = np.fft.fftshift(np.fft.fft(s_f_y))

    # Dispersion
    c = C_LIGHT
    lambda0 = LAMBDA0
    s_f_x_f *= np.exp(1j*(D*lambda0**2/(4*np.pi*c)*L*((2*np.pi)*f)**2))
    s_f_y_f *= np.exp(1j*(D*lambda0**2/(4*np.pi*c)*L*((2*np.pi)*f)**2))

    # Conversion to time domain
    s_f_x = np.fft.ifft(np.fft.ifftshift(s_f_x_f))
    s_f_y = np.fft.ifft(np.fft.ifftshift(s_f_y_f))
    s[0,:] = s_f_x
    s[1,:] = s_f_y

    # Generation of the sync sequence
    # Initial zeros
    s_zeros_init = np.zeros(10,dtype=complex)
    # Sync sequence
    np.random.seed(sync_seed_X)
    s_sync_i_x_d = np.random.randint(0,2,N_sync,dtype=int)*2-1
    np.random.seed(sync_seed_Y)
    s_sync_i_y_d = np.random.randint(0,2,N_sync,dtype=int)*2-1

    # Final zeros
    s_zeros_final = np.zeros(10,dtype=complex)
    s_b_i_x = np.concatenate([s_zeros_init,s_sync_i_x_d,s_zeros_final])
    s_b_i_y = np.concatenate([s_zeros_init,s_sync_i_y_d,s_zeros_final])

    # Digital-to-analogue conversion of the synchronization sequence
    s_sync_i_x = DAC_Nyquist(s=s_b_i_x,SpS=SpS,RollOff=RollOff)
    s_sync_i_y = DAC_Nyquist(s=s_b_i_y,SpS=SpS,RollOff=RollOff)

    # Correlation between signal and synchronization
    sx = s[0,:]
    sy = s[1,:]
    s_xx_i = np.correlate(sx,s_sync_i_x)
    s_xy_i = np.correlate(sx,s_sync_i_y)
    s_yy_i = np.correlate(sy,s_sync_i_y)
    s_yx_i = np.correlate(sy,s_sync_i_x)

    if plot_flag:
        plt.figure()
        plt.plot(abs(s_xx_i))
        plt.plot(abs(s_xy_i))
        plt.plot(abs(s_yy_i))
        plt.plot(abs(s_yx_i))
        plt.ylabel('Cross-correlation')
        plt.xlabel('Sample index')
        plt.xlim([0,len(s_xx_i)])
    s_corr_total = abs(s_xx_i)+abs(s_xy_i)+abs(s_yy_i)+abs(s_yx_i)

    #print(N_sync)
    #kaixo

    # Generating the index of the maximum correlation
    max_i = np.max(s_corr_total)
    ind_i = np.where(max_i==s_corr_total)
    # Finding the inidices of the sync and information
    # F-3: o literal 10 e a guarda do proprio gabarito de sincronismo (N_GUARD_SYNC).
    # A causa real do IndexError era a AUSENCIA de verificacao de limites: com pico
    # de correlacao espurio, ind_inf ultrapassa o fim do buffer.
    peak = int(ind_i[0][0])
    ind_sync = (np.arange(N_sync)+N_GUARD_SYNC)*SpS+peak
    ind_inf  = (np.arange(N_inf)+N_GUARD_SYNC+N_sync)*SpS+peak
    if int(ind_inf[-1]) >= N_samples:
        raise SyncError(
            f"pico de correlacao em {peak} exige {int(ind_inf[-1])+1} amostras, "
            f"mas o sinal tem {N_samples}. Pico provavelmente espurio: excesso de "
            f"ruido, demux mal centrado, ou distancia alem do alcance.")
    #print(N_sync)

    #print(max(abs(s_xx_i[ind_i])))
    #print(max(abs(s_xy_i[ind_i])))
    #print(max(abs(s_yy_i[ind_i])))
    #print(max(abs(s_yx_i[ind_i])))
    #print(max(np.angle(s_xx_i[ind_i])))
    #print(max(np.angle(s_xy_i[ind_i])))
    #print(max(np.angle(s_yy_i[ind_i])))
    #print(max(np.angle(s_yx_i[ind_i])))

    H_hat = np.zeros((2,2),dtype=complex)
    H_hat[0][0] = s_xx_i[ind_i][0]
    H_hat[0][1] = s_xy_i[ind_i][0]
    H_hat[1][0] = s_yx_i[ind_i][0]
    H_hat[1][1] = s_yy_i[ind_i][0]
    H_inv = np.linalg.inv(H_hat)

    #H_hat = np.matrix([[s_xx_i[ind_i][0],s_xy_i[ind_i]][0],[s_yx_i[ind_i][0],s_yy_i[ind_i][0]]])
    #print(H_hat)
    #print(H_inv)

    s_sync = s[:,ind_sync]
    s_sync_x = sx[ind_sync]
    s_sync_y = sy[ind_sync]
    s_inf = s[:,ind_inf]
    s_inf_x = s_inf[0,:]
    s_inf_y = s_inf[1,:]

    s_sync_corr = np.matmul(H_inv,s_sync)
    s_sync_corr_x = s_sync_corr[0,:]
    s_sync_corr_y = s_sync_corr[1,:]

    s_inf_corr = np.matmul(H_inv,s_inf)
    s_inf_corr_x = s_inf_corr[0,:]
    s_inf_corr_y = s_inf_corr[1,:]

    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(np.real(s_sync_corr_x),'.')
    plt.plot(np.imag(s_sync_corr_x),'.')
    plt.plot(np.real(s_inf_corr_x),'.')
    plt.plot(np.imag(s_inf_corr_x),'.')

    plt.subplot(2,1,2)
    plt.plot(np.real(s_sync_corr_y),'.')
    plt.plot(np.imag(s_sync_corr_y),'.')
    plt.plot(np.real(s_inf_corr_y),'.')
    plt.plot(np.imag(s_inf_corr_y),'.')

    if plot_flag:
        # In-phase component time domain
        plt.figure()
        plt.plot(t*1e9,np.real(sx))
        plt.plot(ind_sync*ts*1e9,np.real(s_sync_x),'.')
        plt.plot(ind_inf*ts*1e9,np.real(s_inf_x),'.')
        plt.xlabel('Time [ns]')
        plt.ylabel('In-phase component [a.u.]')
        plt.title('Raw data')
        plt.xlim([0,np.max(t)*1e9])

        # Contellation
        plt.figure()
        plt.plot(np.real(s_inf_corr_x),np.imag(s_inf_corr_x),'.')
        plt.plot(np.real(s_sync_corr_x),np.imag(s_sync_corr_x),'.')
        plt.axis('square')
        plt.xlabel('In-phase [a.u.]')
        plt.ylabel('Quadrature [a.u.]')
        plt.title('Raw data')


    K_amp = np.mean(np.abs(s_sync_corr_x))
    #print(K_amp)

    s_sync_corr_x = s_sync_corr_x/K_amp
    s_sync_corr_y = s_sync_corr_y/K_amp
    s_inf_corr_x = s_inf_corr_x/K_amp
    s_inf_corr_y = s_inf_corr_y/K_amp

    # Contellation
    plt.figure()
    plt.plot(np.real(s_inf_corr_x),np.imag(s_inf_corr_x),'.')
    plt.plot(np.real(s_sync_corr_x),np.imag(s_sync_corr_x),'.')
    plt.axis('square')
    plt.xlabel('In-phase [a.u.]')
    plt.ylabel('Quadrature [a.u.]')
    plt.title('Raw data')

    # Desnormalized the amplitude of the modulated signal
    if M == 16:
        sx_inf_dn = 3*s_inf_corr_x
        sy_inf_dn = 3*s_inf_corr_y
    if M == 32:
        sx_inf_dn = 5*s_inf_corr_x
        sy_inf_dn = 5*s_inf_corr_y
    s_inf_dn = np.append(sx_inf_dn,sy_inf_dn)

    # Blind phase noise and residual frequency offset compensation
    #N_block = 8
    #N_overlap = 4
    N_block = 64
    N_overlap = 32
    #N_block = 32
    #N_overlap = 26

    ind_1 = 0
    ind_2 = ind_1+N_block
    phi = 0
    min_dphi = -np.pi/4
    max_dphi = np.pi/4
    N_dphi = 21
    dphi_array = np.linspace(min_dphi,max_dphi,N_dphi)

    phi_array = np.array([0])
    ind_array = np.array([0])
    counter1 = 0
    while ind_2<=len(sx_inf_dn):
        s_block = sx_inf_dn[ind_1:ind_2]
        d_array_block = np.zeros(N_dphi)
        for counter, dphi in enumerate(dphi_array):
            s_block_rot = s_block*np.exp(1j*(dphi+phi))
            s_hat, dist = QAM_dem(s_block_rot,M)
            d_array_block[counter] = np.var(dist)

        #plt.figure()
        #plt.plot(d_array_block)

        d_min = min(d_array_block)
        ind_min = np.where(d_array_block==d_min)
        dphi_opt = dphi_array[ind_min[0][0]]
        phi = phi+dphi_opt
        counter1 += 1
        phi_array = np.append(phi_array, phi)
        ind_array = np.append(ind_array, int((ind_2+ind_1)/2))
        ind_1 = ind_1+N_block-N_overlap
        ind_2 = ind_1+N_block

    ind_array = np.concatenate(([0],ind_array,[len(sx_inf_dn)]))
    phi_array = np.concatenate(([0],phi_array,[phi_array[-1]]))
    phi_interp = np.interp(np.arange(len(sx_inf_dn)),ind_array,phi_array)

    if plot_flag:
        plt.figure()
        plt.plot(phi_interp)
        #print(len(phi_interp))
    sx_inf_pc = sx_inf_dn*np.exp(1j*(phi_interp))
    sy_inf_pc = sy_inf_dn*np.exp(1j*(phi_interp))

    if plot_flag:
        # Contellation
        Ni = 10
        plt.figure()
        plt.plot(np.real(sx_inf_pc[Ni:]),np.imag(sx_inf_pc[Ni:]),'.')
        plt.plot(np.real(sy_inf_pc[Ni:]),np.imag(sy_inf_pc[Ni:]),'.')
        plt.axis('square')
        plt.xlabel('In-phase [a.u.]')
        plt.ylabel('Quadrature [a.u.]')
        plt.title('After scaling correction')

    # Demapping
    s_hat_x, dist = QAM_dem(sx_inf_pc,M)
    s_hat_y, dist = QAM_dem(sy_inf_pc,M)

    s_b_x = s_b[0,:]
    s_b_y = s_b[1,:]
    #print(s_b_x)
    #print(s_b_y)

    if plot_flag:
        plt.figure()
        plt.stem(s_hat_x!=s_b_x)
        plt.stem(s_hat_y!=s_b_y)

    BER_X = np.mean(s_hat_x!=s_b_x)
    BER_Y = np.mean(s_hat_y!=s_b_y)
    BER = np.array([BER_X,BER_Y])

    return [s_inf, s_sync, s_inf_dn, BER]


def DPQAM_receiver(E_RX, t, s_b, M, N_inf, SpS, RollOff, ts, BW_optical = None, P_laser = 0.01, Delta_nu = 100e3, N_sync= 128, sync_seed_X=0, sync_seed_Y=123,Freq_offset = 1e6, plot_flag = False,\
                  R=0.9, BW=4*20e9, shot_noise = True, thermal_noise = True,L=10e3,D=16e-6):

    #------------------------------------- Optical front end --------------------------------------------%

    _rec('rx_field_in', E_RX[0,:])

    E_RXf = optical_filter(E_RX,t,f0 = Freq_offset, BW = BW if BW_optical is None else BW_optical)

    _rec('rx_field_filtered', E_RXf[0,:])

    # Separate polarizations
    E_RX_X = E_RXf[0,:]
    E_RX_Y = E_RXf[1,:]




    # LO laser
    E_LO = single_frequency_laser(t = t, P = P_laser, Delta_nu = Delta_nu, Freq_offset = Freq_offset, X_fraction = 0.5, phi_pol = 0)



    # Separate the two polarizations of LO
    E_LO_X = E_LO[0,:]
    E_LO_Y = E_LO[1,:]

    # Process X polarization
    # Combine signal and LO
    E_1, E_2, E_3, E_4 = HybridNetwork(E_RX_X,E_LO_X)
    # Detect the outputs of 90-degree hybrid

    ipd_1 = photodiode(s=E_1,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)

    #plt.figure()
    #plt.plot(abs(ipd_1))
    ipd_2 = photodiode(s=E_2,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_3 = photodiode(s=E_3,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_4 = photodiode(s=E_4,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    s_X_I = ipd_1-ipd_2   # Differential amplification
    s_X_Q = ipd_3-ipd_4   # Differential amplification

    #S_aux = np.fft.fftshift(np.fft.fft(s_X_I))
    #plt.figure()
    #plt.plot(np.abs(S_aux))
    #plt.plot(s_X_I)

    # Process Y polarization
    # Combine signal and LO
    E_1, E_2, E_3, E_4 = HybridNetwork(E_RX_Y,E_LO_Y)
    # Detect the outputs of 90-degree hybrid
    ipd_1 = photodiode(s=E_1,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_2 = photodiode(s=E_2,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_3 = photodiode(s=E_3,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    ipd_4 = photodiode(s=E_4,t=t,R=R,BW=BW,shot_noise=shot_noise,thermal_noise=thermal_noise)
    s_Y_I = ipd_1-ipd_2  # Differential amplification
    s_Y_Q = ipd_3-ipd_4  # Differential amplification

    # Cartesian to complex conversion
    s_rx_X = s_X_I+1j*s_X_Q
    s_rx_Y = s_Y_I+1j*s_Y_Q

    #print(np.shape(s_Y_I))

    n = np.shape(s_Y_I)[-1]          # C-13: shape nao depende mais das flags de ruido
    s_rec = np.zeros((2,n),dtype=complex)

    #print(np.shape(s_Y_I))
    #print(np.shape(s_rec))
    #print(np.shape(s_rx_X))

    s_rec[0,:] = s_rx_X
    s_rec[1,:] = s_rx_Y

    #QAM_receiver_DP(s_rec, t = t, M =M, SpS = SpS, N_sync = N_sync, N_inf = N_inf,sync_seed_X=0,sync_seed_Y=123,\
    #                RollOff = RollOff,ts = ts,s_b=s_b[0,:],plot_flag = plot_flag)
    s_inf, s_sync, s_inf_dn, BER = QAM_receiver_DP(s_rec, t = t, M =M, SpS = SpS, N_sync = N_sync, N_inf = N_inf,sync_seed_X=sync_seed_X,sync_seed_Y=sync_seed_Y,\
                    RollOff = RollOff,ts = ts,s_b=s_b,plot_flag = plot_flag,L=L,D=D)
    return BER


def constellation_analysis(s):
    plt.figure()
    plt.plot(np.real(s),np.imag(s),'.')

    H, xedges, yedges = np.histogram2d(np.real(s),np.imag(s),100)
    #plt.figure()
    #plt.contourf(counts)

    fig = plt.figure()
    #ax = fig.add_subplot(131, title='imshow: square bins')
    plt.imshow(np.log10(H), interpolation='nearest', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],cmap="hot")

    fig = plt.figure()
    #ax = fig.add_subplot(131, title='imshow: square bins')
    plt.imshow(H,cmap="hot")

    fig = plt.figure()
    #ax = fig.add_subplot(131, title='imshow: square bins')
    plt.imshow(np.log10(H), interpolation='nearest', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],cmap="Blues")

    fig = plt.figure()
    #ax = fig.add_subplot(131, title='imshow: square bins')
    plt.imshow(H,cmap="Blues")


def spectrum_analysis(s,t):
    ts = t[1]-t[0]                       # Sampling period in s
    fs = 1/ts                            # Sampling frequancy in Hz
    Deltaf = 1/(ts*len(t))               # Frequency resolution in Hz
    f = np.arange(len(t))*Deltaf-fs/2    # Frequency vector in Hz

    plt.figure()
    plt.plot(f/1e9,30+10*np.log10(np.fft.fftshift(np.abs(np.fft.fft(s)))))
    plt.xlabel('Frequency [GHz]')
    plt.ylabel('Power [dBm]')


# ============================================================================
# Cenarios
# ============================================================================
@dataclass
class Scenario:
    # --- grade temporal (unico grau de liberdade real: fs = BaudRate*SpS) ---
    BaudRate: float = 14e9
    SpS: int = 16
    M: int = 16
    RollOff: float = 0.2
    # --- quadro ---
    N_zeros_init: int = 100
    N_sync: int = 256
    N_inf: int = 4*4096
    N_zeros_final: int = 100
    # --- transmissor ---
    P_laser_TX: float = 0.01
    Delta_nu_TX: float = 100e3
    ind_mod: float = 0.1
    # --- receptor ---
    P_laser_RX: float = 0.01
    Delta_nu_RX: float = 100e3
    Freq_offset_LO: float = 55e6
    R: float = 0.9
    BW_elec: float = 4*20e9
    shot_noise: bool = True
    thermal_noise: bool = True
    # --- enlace ---
    L_span: float = 50e3
    N_span: int = 1
    alpha_dB: float = 0.2e-3
    D: float = 16e-6
    n2: float = 2.6e-20
    Aeff: float = 80e-12
    DeltaL: float = 1e3
    SoP_rotation: bool = True
    NF_dB: float = 5.0
    P_LOP: float = 0.002
    # --- WDM ---
    n_channels: int = 1
    channel_spacing: float = 75e9
    BW_optical: float = 75e9
    # --- reprodutibilidade ---
    root_seed: int = 0

    @property
    def ts(self) -> float:   return 1.0/(self.BaudRate*self.SpS)
    @property
    def fs(self) -> float:   return self.BaudRate*self.SpS
    @property
    def n_symbols(self) -> int:
        return self.N_zeros_init + self.N_sync + self.N_inf + self.N_zeros_final
    @property
    def n_samples(self) -> int:  return self.n_symbols*self.SpS
    @property
    def window(self) -> float:   return self.n_samples*self.ts

    def channel_offsets(self) -> list[float]:
        k = np.arange(self.n_channels) - (self.n_channels - 1)/2
        return list(k*self.channel_spacing)

    def check_nyquist(self) -> None:
        """Versao minima de RequiresNyquistSatisfied (secao 7 do documento)."""
        b_half = self.BaudRate*(1 + self.RollOff)/2
        edge = max(abs(f0) for f0 in self.channel_offsets()) + b_half
        if edge > self.fs/2:
            raise ValueError(
                f"Nyquist violado: a grade ocupa +-{edge/1e9:.1f} GHz mas a banda "
                f"de simulacao e +-{self.fs/2e9:.1f} GHz. Aumente SpS ou reduza "
                f"o espacamento de canais.")
        print(f"  [nyquist] grade +-{edge/1e9:.1f} GHz em banda +-{self.fs/2e9:.1f} GHz "
              f"(ocupacao {200*edge/self.fs:.0f}%)")


def transmit_channel(sc: Scenario, ch: int):
    """Um transmissor. C-1: cada canal tem suas proprias seeds, portanto bits e
    ruido de fase independentes. No codigo original os 4 canais WDM usavam
    sync_seed_X=0 e sync_seed_Y=123 e carregavam a MESMA sequencia de bits."""
    seed_x, seed_y = sc.root_seed + 2*ch, sc.root_seed + 2*ch + 1
    f0 = sc.channel_offsets()[ch]
    t, E, E_carrier, s_tx, s_b = DPQAM_transmitter(
        P_laser=sc.P_laser_TX, M=sc.M, SpS=sc.SpS, RollOff=sc.RollOff, ts=sc.ts,
        sync_seed_X=seed_x, sync_seed_Y=seed_y, N_sync=sc.N_sync, N_inf=sc.N_inf,
        N_zeros_init=sc.N_zeros_init, N_zeros_final=sc.N_zeros_final,
        ind_mod=sc.ind_mod, Delta_nu=sc.Delta_nu_TX, Freq_offset=f0, plot_flag=False)
    return t, E, s_b, (seed_x, seed_y), f0


def set_launch_power(E, sc: Scenario, t):
    """Controle de potencia via EDFA, com a potencia medida corretamente (C-5)."""
    p_now = total_power(E) if not LEGACY_PHYSICS else \
            float(np.mean(np.abs(E[0, :])**2)*2)     # bug P-2 preservado no legacy
    gain_db = 10*np.log10(sc.P_LOP/p_now)
    return EDFA(E, t, gain_db, sc.NF_dB, lambda0=LAMBDA0), gain_db


def receive_channel(sc: Scenario, E_RX, t, s_b, seeds, f0, L_acc):
    seed_x, seed_y = seeds
    return DPQAM_receiver(
        E_RX, t, s_b, sc.M, sc.N_inf, sc.SpS, sc.RollOff, sc.ts,
        BW_optical=sc.BW_optical if sc.n_channels > 1 else None,
        P_laser=sc.P_laser_RX, Delta_nu=sc.Delta_nu_RX, N_sync=sc.N_sync,
        sync_seed_X=seed_x, sync_seed_Y=seed_y,
        Freq_offset=f0 + sc.Freq_offset_LO, plot_flag=False, R=sc.R,
        BW=sc.BW_elec, shot_noise=sc.shot_noise, thermal_noise=sc.thermal_noise,
        L=L_acc, D=sc.D)


def propagate(E, t, sc: Scenario, n_spans: int, pmd: PMDConfig):
    """O enlace: SSFM por span, EDFA por span. Campo agregado, um solver so."""
    for _ in range(n_spans):
        E = fiber(E, t, sc.L_span, sc.DeltaL, sc.D, LAMBDA0, sc.alpha_dB,
                  sc.n2, sc.Aeff, sc.SoP_rotation, pmd=pmd)
        E = EDFA(E, t, sc.alpha_dB*sc.L_span, sc.NF_dB, lambda0=LAMBDA0)
    return E


# ============================================================================
def run(sc: Scenario, scenario: str):
    np.random.seed(sc.root_seed)          # legado: RNG global. Ver secao 8.
    sc.check_nyquist()
    pmd = PMDConfig(total_length=sc.L_span*max(1, sc.N_span))
    results: dict = {"scenario": scenario, "legacy_physics": LEGACY_PHYSICS,
                     "config": asdict(sc)}

    # --- transmissores + MUX ---
    chans = [transmit_channel(sc, k) for k in range(sc.n_channels)]
    t = chans[0][0]
    E = sum(c[1] for c in chans)                       # multiplexer = soma
    _rec("tx_field_x", E[0, :])
    E, gain_db = set_launch_power(E, sc, t)
    print(f"  [tx] {sc.n_channels} canal(is), {sc.n_symbols} simbolos, "
          f"{sc.n_samples} amostras, janela {sc.window*1e9:.1f} ns")
    print(f"  [tx] potencia de lancamento {1e3*total_power(E):.3f} mW "
          f"(ganho aplicado {gain_db:+.2f} dB)")

    if scenario == "b2b":
        spans = 0
    else:
        spans = sc.N_span

    ber_per_span = []
    E_link = E
    for s in range(max(1, spans)):
        if spans:
            E_link = propagate(E_link, t, sc, 1, pmd)
        L_acc = sc.L_span*(s + 1) if spans else 0.0
        row = []
        for k, (_, _, s_b, seeds, f0) in enumerate(chans):
            try:
                ber = receive_channel(sc, E_link.copy(), t, s_b, seeds, f0, L_acc)
                row.append([float(ber[0]), float(ber[1])])
                print(f"  [rx] span {s+1 if spans else 0} canal {k}: "
                      f"BER_X={ber[0]:.2e}  BER_Y={ber[1]:.2e}")
            except SyncError as e:
                row.append([float("nan"), float("nan")])
                print(f"  [rx] span {s+1} canal {k}: sincronismo falhou -> {e}")
        ber_per_span.append(row)
        if not spans:
            break

    results["ber"] = ber_per_span
    _rec("ber", np.array(ber_per_span))
    return results


def sweep(sc: Scenario, what: str):
    values = (np.arange(0.001, 0.011, 0.001) if what == "power"
              else np.arange(50e3, 210e3, 25e3))
    out = []
    for v in values:
        if what == "power": sc.P_LOP = float(v)
        else:               sc.L_span, sc.N_span = float(v), 1
        print(f"\n--- {what} = {v:.4g} ---")
        r = run(sc, "span")
        out.append({what: float(v), "ber": r["ber"]})
    return {"scenario": f"sweep-{what}", "legacy_physics": LEGACY_PHYSICS, "points": out}


def main():
    global LEGACY_PHYSICS
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--scenario", default="b2b",
                    choices=["b2b", "span", "wdm", "sweep-power", "sweep-length"])
    ap.add_argument("--fixed", action="store_true",
                    help="aplica as correcoes de fisica (LEGACY_PHYSICS=False)")
    ap.add_argument("--channels", type=int, default=1)
    ap.add_argument("--spacing", type=float, default=None, help="espacamento WDM em Hz")
    ap.add_argument("--sps", type=int, default=16)
    ap.add_argument("--baudrate", type=float, default=14e9)
    ap.add_argument("--spans", type=int, default=1)
    ap.add_argument("--symbols", type=int, default=4*4096)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="resultados")
    a = ap.parse_args()

    LEGACY_PHYSICS = not a.fixed
    sc = Scenario(n_channels=a.channels if a.scenario == "wdm" else 1,
                  N_span=a.spans, N_inf=a.symbols, root_seed=a.seed,
                  SpS=a.sps, BaudRate=a.baudrate)
    if a.spacing:
        sc.channel_spacing = a.spacing
        sc.BW_optical = a.spacing
    if a.scenario == "wdm":
        sc.n_channels = a.channels

    print(f"=== {a.scenario} | LEGACY_PHYSICS={LEGACY_PHYSICS} ===")
    res = (sweep(sc, a.scenario.split("-")[1]) if a.scenario.startswith("sweep")
           else run(sc, a.scenario))

    os.makedirs(a.out, exist_ok=True)
    tag = f"{a.scenario}_{'legacy' if LEGACY_PHYSICS else 'fixed'}"
    with open(os.path.join(a.out, f"{tag}.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    save_record(os.path.join(a.out, f"{tag}.npz"))
    print(f"\nresultados -> {a.out}/{tag}.json  e  {a.out}/{tag}.npz")
    print(f"chaves gravadas para plotagem: {sorted(RECORD)}")


if __name__ == "__main__":
    main()