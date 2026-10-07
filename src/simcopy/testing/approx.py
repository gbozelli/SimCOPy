"""pytest.approx com tolerancia absoluta ZERO por padrao.

pytest.approx tem tolerancia absoluta padrao de 1e-12. Em unidades SI, varias
grandezas deste projeto sao dessa ordem ou menores: DGD ~1e-13 s, beta2 ~2e-26
s^2/m, FWHM ~3e-11 s. Com o padrao, qualquer valor dentro de +-1e-12 "passa" e o
teste vira decorativo: um DGD 3x errado passava. Aqui so a tolerancia relativa
vale, a menos que o teste peca `abs` explicitamente.
"""
import pytest

def approx(expected, rel=1e-6, abs=0.0):
    return pytest.approx(expected, rel=rel, abs=abs)
