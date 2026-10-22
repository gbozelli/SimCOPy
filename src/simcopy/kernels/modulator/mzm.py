"""Modulador de Mach-Zehnder (MZM) e modulador em quadratura (docs/moduladores.md).

Agrawal, Sistemas de Comunicacao por Fibra Optica, secao 10.1.2:
    phi_j = pi*V_j/V_pi                                   (texto antes da Eq. 10.1.3)
    A_t = (1/2)*A_i*[exp(i*phi_1) + exp(i*phi_2)]         (Eq. 10.1.3)
    t_m = cos[(phi_1 - phi_2)/2]*exp[i*(phi_1 + phi_2)/2]  (Eq. 10.1.4)
O modulador em quadratura (Fig. 10.3) poe um MZM em push-pull em cada braco de um
interferometro externo, com defasagem entre os bracos.

Os bracos aceitam amplitudes diferentes (a_1, a_2), como no codigo do professor (perda de
cada braco em dB); com a_1 = a_2 = 1/2 volta a Eq. 10.1.3. Aqui so ha fisica: nada de
configuracao, perda de insercao ou filtro S21.
"""
from __future__ import annotations
import numpy as np

def mzm_field(field, v1, v2, vpi, a1=0.5, a2=0.5):
    """Campo na saida do MZM com tensoes v1, v2 [V] nos bracos (Eq. 10.1.3)."""
    return field * (a1 * np.exp(1j * np.pi * v1 / vpi) + a2 * np.exp(1j * np.pi * v2 / vpi))

def push_pull(v, v_bias):
    """Tensoes dos bracos em push-pull com polarizacao simetrica: (v - V_b/2, -v + V_b/2).

    phi_1 + phi_2 = 0, entao a saida nao tem chirp (texto apos a Eq. 10.1.4), e a
    transmissao e cos^2[pi*(2v - V_b)/(2*V_pi)] (Eq. 10.1.5). O Agrawal poe a polarizacao
    num braco so (V_2 = -V_1 + V_b), o que so acrescenta a fase constante pi*V_b/(2*V_pi);
    a forma simetrica e a do codigo do professor.
    """
    return v - 0.5 * v_bias, -v + 0.5 * v_bias

def iq_field(field, v_i, v_q, vpi, v_bias, quadrature_phase, a1=0.5, a2=0.5):
    """Modulador em quadratura (Fig. 10.3): (1/2)*[t(v_i) + exp(i*theta)*t(v_q)]*field.

    t e o MZM em push-pull de cada braco e theta = `quadrature_phase` [rad], a defasagem
    entre os bracos (pi/2 no ideal). O 1/2 e a divisao e a recombinacao do interferometro
    externo, como na Eq. 10.1.3; o codigo do professor da o mesmo fator (1/sqrt(2) duas vezes).
    """
    t_i = mzm_field(1.0, *push_pull(v_i, v_bias), vpi, a1, a2)
    t_q = mzm_field(1.0, *push_pull(v_q, v_bias), vpi, a1, a2)
    return 0.5 * field * (t_i + np.exp(1j * quadrature_phase) * t_q)
