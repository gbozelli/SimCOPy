"""Seeds: reprodutivel quando deve ser igual, independente quando deve ser diferente.
Substitui o np.random.seed global do rascunho, que fazia os 4 canais WDM
carregarem a mesma sequencia de bits."""
import numpy as np
from simcopy.core.grid import TimeGrid
from simcopy.core.context import SimContext

CTX = SimContext(TimeGrid(1e12, 16), root_seed=0)

def test_mesma_chave_mesma_sequencia():
    a = CTX.rng_for(2, "payload").integers(0, 2, 1000)
    b = CTX.rng_for(2, "payload").integers(0, 2, 1000)
    assert np.array_equal(a, b)

def test_canais_diferentes_sao_independentes():
    a = CTX.rng_for(0, "payload").integers(0, 2, 10000)
    b = CTX.rng_for(1, "payload").integers(0, 2, 10000)
    assert not np.array_equal(a, b)
    assert abs(np.corrcoef(a, b)[0, 1]) < 0.05

def test_payload_e_sincronismo_sao_independentes():
    a = CTX.rng_for(0, "payload").integers(0, 2, 10000)
    b = CTX.rng_for(0, "sync").integers(0, 2, 10000)
    assert not np.array_equal(a, b)

def test_seed_raiz_diferente_muda_tudo():
    outro = SimContext(TimeGrid(1e12, 16), root_seed=1)
    assert not np.array_equal(CTX.rng_for(0, "payload").integers(0, 2, 100),
                              outro.rng_for(0, "payload").integers(0, 2, 100))
