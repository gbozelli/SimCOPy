"""Mistura do sinal com o oscilador local: hibrida de 90 graus (docs/sistema_coerente.md, secao 3).

Agrawal, Sistemas de Comunicacao por Fibra Optica:
    P(t) = P_s + P_LO + 2*sqrt(P_s*P_LO)*cos(w_IF*t + phi_s - phi_LO)        (Eq. 4.5.3)
num fotodetector so: o batimento vem junto com P_s e P_LO. No receptor balanceado,
    I+- = (1/2)*R_d*(P_s + P_LO) +- R_d*sqrt(P_s*P_LO)*cos(w_IF*t + phi)      (Eq. 10.2.4)
e a diferenca I+ - I- elimina o termo constante e deixa so o batimento. A hibrida de 90 graus com
quatro fotodiodos recupera as partes real e imaginaria de E_s*conj(E_LO) (secao 8.7.3, Fig. 8.39).

Codigo do professor (draft.py, `HybridNetwork`):
    A_1 = (A_s + A_LO)/2      A_2 = (A_s - A_LO)/2
    A_3 = (A_s + i*A_LO)/2    A_4 = (A_s - i*A_LO)/2
A soma das quatro potencias e P_s + P_LO (sem perda). Com fotodiodos de responsividade R,
    R*(|A_1|^2 - |A_2|^2) = R*Re[A_s*conj(A_LO)]
    R*(|A_3|^2 - |A_4|^2) = R*Im[A_s*conj(A_LO)]
entao I + i*Q = R*A_s*conj(A_LO): o campo do sinal, com a fase medida em relacao ao oscilador
local. Vale para envelopes, entao nao depende da convencao da portadora. Um oscilador local em
+df (azul), cujo envelope e exp(-i*2*pi*df*t) na convencao do projeto, faz a constelacao girar
como exp(+i*2*pi*df*t).
"""
from __future__ import annotations
import numpy as np

def hybrid_90(signal, lo) -> np.ndarray:
    """Quatro saidas (4, N) da hibrida de 90 graus para os campos `signal` e `lo` (N,), na mesma
    polarizacao. As saidas 0 e 1 formam o par em fase; as saidas 2 e 3, o par em quadratura."""
    s = np.asarray(signal, dtype=complex)
    l = np.asarray(lo, dtype=complex)
    return 0.5 * np.stack([s + l, s - l, s + 1j * l, s - 1j * l])
