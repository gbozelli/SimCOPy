import sys
import os
import numpy as np
import pytest

# Garante que o pytest encontre os arquivos na pasta raiz
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Assumindo que você salvou a função otimizada em um arquivo chamado qam_tx.py
from qam_tx import qam_tx 

def test_qam_tx_dimensions():
    """
    Testa se as dimensões dos arrays de saída correspondem à matemática teórica dos parâmetros.
    """
    N_inf = 1024
    M = 16
    
    # Executa a função desabilitando o gráfico para o teste não travar a tela
    s_cont, t, s_b = qam_tx(N_inf=N_inf, M=M, plot_flag=False)
    
    # Verificação 1: O número total de bits gerados deve ser N_inf * log2(M)
    # Ex: 1024 * 4 (para 16-QAM) = 4096 bits
    expected_bits = int(N_inf * np.log2(M))
    assert len(s_b) == expected_bits, f"Esperado {expected_bits} bits, mas obteve {len(s_b)}"
    
    # Verificação 2: O vetor de tempo deve ter o mesmo tamanho do sinal contínuo
    assert len(s_cont) == len(t), "Vetor de tempo (t) incompatível com sinal contínuo (s_cont)"

def test_qam_tx_determinism():
    """
    Testa se a função é reprodutível ao fixarmos as seeds de randomização.
    """
    # Fixamos a seed do numpy antes da chamada para a geração de bits (s_b)
    np.random.seed(42)
    s_cont1, t1, s_b1 = qam_tx(sync_seed=10, plot_flag=False)
    
    # Resetamos a seed global para o exato mesmo estado inicial
    np.random.seed(42)
    s_cont2, t2, s_b2 = qam_tx(sync_seed=10, plot_flag=False)
    
    # Verifica se todos os arrays resultantes são estritamente iguais
    np.testing.assert_array_equal(s_b1, s_b2)
    np.testing.assert_array_equal(s_cont1, s_cont2)
    np.testing.assert_array_equal(t1, t2)

@pytest.mark.parametrize("MIMO_i", [0, 1])
def test_qam_tx_mimo_execution(MIMO_i):
    """
    Garante que as duas lógicas de construção da sequência MIMO funcionam sem estourar erros.
    """
    try:
        s_cont, t, s_b = qam_tx(MIMO_i=MIMO_i, plot_flag=False)
        assert s_cont is not None
    except Exception as e:
        pytest.fail(f"A função falhou ao processar MIMO_i={MIMO_i}. Erro: {e}")