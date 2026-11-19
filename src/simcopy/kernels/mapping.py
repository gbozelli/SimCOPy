"""Mapeamento de bits em simbolos e decisao: BPSK, QPSK, 16QAM e 64QAM (docs/sistema_coerente.md,
secoes 2 e 12).

Agrawal, Sistemas de Comunicacao por Fibra Optica, secao 10.1: no BPSK a fase da portadora vale 0
ou pi; no QPSK o modulador em quadratura da os quatro pontos (+-1 +- i)/sqrt(2) (Fig. 10.1 e 10.3).

Convencao do codigo do professor (draft.py, `QAM_dem`, M = 4): pontos normalizados pelo maior nivel
(+-1 em cada quadratura) e, no QPSK, o primeiro bit do simbolo escolhe a quadratura (0 -> +1,
1 -> -1) e o segundo a fase (0 -> +1, 1 -> -1). E um codigo de Gray: pontos vizinhos diferem em um
bit. O BPSK usa a mesma regra em um eixo so: 0 -> +1, 1 -> -1.

16QAM (draft.py, `QAM_dem`, M = 16): os dois primeiros bits do simbolo escolhem o nivel em fase e
os dois ultimos o nivel em quadratura, com o mesmo codigo de Gray nos dois eixos: 00 -> -3,
01 -> -1, 11 -> +1, 10 -> +3. Os niveis sao divididos por 3 (o professor divide pelo maior).

64QAM: o draft.py nao tem. Segue a mesma regra do 16QAM com tres bits por eixo: a primeira metade
dos bits escolhe o nivel em fase, a segunda o nivel em quadratura, com o codigo de Gray refletido
(000 -> -7, 001 -> -5, 011 -> -3, 010 -> -1, 110 -> +1, 111 -> +3, 101 -> +5, 100 -> +7), niveis
divididos por 7. O 00, 01, 11, 10 do 16QAM e o mesmo codigo com dois bits.

Os niveis saem em [-1, 1], como a pre-distorcao do modulador espera (kernels/modulator/mzm.py).
O 32QAM (constelacao em cruz) ainda nao entra.
"""
from __future__ import annotations
import numpy as np

BITS_PER_SYMBOL = {"BPSK": 1, "QPSK": 2, "16QAM": 4, "64QAM": 6}
BINARY_PER_QUADRATURE = ("BPSK", "QPSK")     # cada bit e uma decisao entre dois niveis
SQUARE_QAM = ("16QAM", "64QAM")              # grade quadrada com Gray por eixo
ROTATIONAL_SYMMETRY = {"BPSK": 2, "QPSK": 4, "16QAM": 4, "64QAM": 4}   # giros que repetem a constelacao

def gray_to_index(bits) -> np.ndarray:
    """Indice do nivel (0 = mais baixo) de cada palavra de Gray refletido. `bits`: (n, k), bit mais
    significativo primeiro. Com k = 2: 00 -> 0, 01 -> 1, 11 -> 2, 10 -> 3."""
    b = np.asarray(bits).astype(int)
    idx = np.zeros(b.shape[0], dtype=int)
    acc = np.zeros(b.shape[0], dtype=int)
    for j in range(b.shape[1]):               # binario: cada bit e o xor dos bits de Gray ate ele
        acc = acc ^ b[:, j]
        idx = (idx << 1) | acc
    return idx

def index_to_gray(index, k: int) -> np.ndarray:
    """Palavra de Gray refletido (n, k) do indice do nivel: o inverso de `gray_to_index`."""
    idx = np.asarray(index).astype(int)
    g = idx ^ (idx >> 1)
    return ((g[:, None] >> np.arange(k - 1, -1, -1)) & 1).astype(np.uint8)

def pam_levels(bits, k: int) -> np.ndarray:
    """Niveis de um eixo com 2^k niveis de Gray, de -1 a 1: (2*indice - (2^k - 1))/(2^k - 1)."""
    m = 2 ** k
    return (2.0 * gray_to_index(bits) - (m - 1)) / (m - 1)

def pam_decide(values, k: int) -> np.ndarray:
    """Indice do nivel mais proximo num eixo com 2^k niveis de -1 a 1."""
    m = 2 ** k
    return np.clip(np.round((np.asarray(values) * (m - 1) + (m - 1)) / 2.0), 0, m - 1).astype(int)

def _check(name: str) -> int:
    if name not in BITS_PER_SYMBOL:
        raise ValueError(f"formato {name!r} sem mapeamento; disponiveis: {sorted(BITS_PER_SYMBOL)}")
    return BITS_PER_SYMBOL[name]

def constellation(name: str) -> np.ndarray:
    """Pontos da constelacao, na ordem do valor inteiro dos bits de cada simbolo."""
    k = _check(name)
    bits = (np.arange(2 ** k)[:, None] >> np.arange(k - 1, -1, -1)) & 1
    return map_bits(bits.ravel(), name)

def map_bits(bits, name: str) -> np.ndarray:
    """Bits (0/1), `bits_per_symbol` por simbolo em sequencia -> simbolos complexos."""
    k = _check(name)
    b = np.asarray(bits).astype(int)
    if b.size % k:
        raise ValueError(f"{b.size} bits nao formam simbolos de {k} bits ({name})")
    if name == "BPSK":
        return (1 - 2 * b).astype(complex)
    if name in SQUARE_QAM:
        b = b.reshape(-1, k)
        return pam_levels(b[:, :k // 2], k // 2) + 1j * pam_levels(b[:, k // 2:], k // 2)
    b = b.reshape(-1, 2)
    return (1 - 2 * b[:, 1]) + 1j * (1 - 2 * b[:, 0])

def demap(symbols, name: str) -> np.ndarray:
    """Decisao pelo ponto mais proximo: simbolos complexos -> bits, na ordem de `map_bits`.
    Os pontos formam uma grade, entao o ponto mais proximo sai de cada quadratura em separado."""
    k = _check(name)
    s = np.asarray(symbols)
    if name == "BPSK":
        return (s.real < 0).astype(np.uint8)
    if name in SQUARE_QAM:
        h = k // 2
        return np.hstack([index_to_gray(pam_decide(s.real, h), h),
                          index_to_gray(pam_decide(s.imag, h), h)]).astype(np.uint8).ravel()
    return np.column_stack([s.imag < 0, s.real < 0]).astype(np.uint8).ravel()

def nearest_points(symbols, name: str) -> np.ndarray:
    """Ponto da constelacao mais proximo de cada simbolo recebido. Sai direto dos niveis, sem
    passar pelos bits, e por isso aceita qualquer forma de entrada, por exemplo (n_pol, n).
    (`demap` achata a saida em bits por simbolo, e so serve para uma sequencia.)"""
    k = _check(name)
    s = np.asarray(symbols)
    if name in SQUARE_QAM:
        m = 2 ** (k // 2)
        level = lambda v: (2.0 * pam_decide(v, k // 2) - (m - 1)) / (m - 1)
        return level(s.real) + 1j * level(s.imag)
    sign = lambda v: np.where(v < 0, -1.0, 1.0)
    if name == "BPSK":
        return sign(s.real).astype(complex)
    return sign(s.real) + 1j * sign(s.imag)

def decision_components(symbols, bits, name: str) -> tuple[np.ndarray, np.ndarray]:
    """(valores, nivel alto): a componente real que decide cada bit e se o bit transmitido esta no
    nivel alto (+1). No BPSK e no QPSK cada bit e uma decisao binaria numa quadratura, o que
    permite o fator Q do Agrawal (Eq. 10.3.2 e 10.3.7)."""
    if name not in BINARY_PER_QUADRATURE:
        raise ValueError(f"{name}: o fator Q de decisao binaria so vale para {BINARY_PER_QUADRATURE}")
    s = np.asarray(symbols)
    b = np.asarray(bits).astype(int)
    if name == "BPSK":
        return s.real, b == 0
    b = b.reshape(-1, 2)
    return np.concatenate([s.imag, s.real]), np.concatenate([b[:, 0] == 0, b[:, 1] == 0])
