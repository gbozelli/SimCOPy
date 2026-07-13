import sys
import os
import numpy as np
import pytest

# Garante que o pytest encontre os arquivos na pasta raiz e/ou na pasta deprecated
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Ajuste os imports abaixo conforme o nome exato dos seus arquivos e pastas
from QAM_transmitter import QAM_transmitter # Importa a versão antiga
from qam_tx import qam_tx                   # Importa a versão nova/otimizada

def test_deprecated_qam_transmitter_dimensions():
    """
    Testa se as dimensões geradas pela função antiga estão de acordo com o esperado 
    (lembrando que ela omite a sequência MIMO na concatenação final).
    """
    N_inf = 1024
    M = 16
    SpS = 16
    
    # Executa a função antiga desabilitando o gráfico
    s_cont, t, s_b = QAM_transmitter(N_inf=N_inf, M=M, SpS=SpS, plot_flag=False)
    
    # Verificação dos bits gerados
    expected_bits = int(N_inf * np.log2(M))
    assert len(s_b) == expected_bits, f"Esperado {expected_bits} bits, mas obteve {len(s_b)}"
    
    # Verificação da compatibilidade de tamanho entre sinal contínuo e vetor de tempo
    assert len(s_cont) == len(t), "Vetor de tempo incompatível com o sinal contínuo na versão deprecated"

def test_deprecated_determinism():
    """
    Garante que a versão antiga sempre gera a mesma saída dado o mesmo estado inicial (seed).
    """
    np.random.seed(100)
    s_cont1, t1, s_b1 = QAM_transmitter(sync_seed=5, plot_flag=False)
    
    np.random.seed(100)
    s_cont2, t2, s_b2 = QAM_transmitter(sync_seed=5, plot_flag=False)
    
    np.testing.assert_array_equal(s_b1, s_b2)
    np.testing.assert_array_equal(s_cont1, s_cont2)

def test_compare_bits_new_vs_deprecated():
    """
    Compara a versão NOVA com a versão ANTIGA.
    Como a versão antiga omitia o MIMO no sinal contínuo (s_cont), comparamos 
    apenas o vetor de bits da informação (s_b), que deve ser estritamente igual 
    se o ambiente randômico for o mesmo.
    """
    N_inf = 512
    M = 64
    sync_seed = 42
    
    # Gera resultados com a versão antiga
    np.random.seed(99)
    _, _, s_b_old = QAM_transmitter(N_inf=N_inf, M=M, sync_seed=sync_seed, plot_flag=False)
    
    # Gera resultados com a versão nova
    np.random.seed(99)
    _, _, s_b_new = qam_tx(N_inf=N_inf, M=M, sync_seed=sync_seed, plot_flag=False)
    
    # O vetor de bits de informação deve ser exatamente o mesmo
    np.testing.assert_array_equal(s_b_new, s_b_old, err_msg="Os bits gerados divergem entre a versão nova e a deprecated!")