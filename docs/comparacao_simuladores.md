# O simcopy e os outros simuladores abertos (2026-11-15)

O simcopy já tem três coisas que nenhuma das páginas lidas mostra junta: a comparação ponto a ponto
com um simulador comercial, com os dados no repositório, o transmissor IMDD por equações de taxa e o
coerente no mesmo pacote, e o material de ensino em português. Em compensação está atrás em quatro
frentes: testes, DSP cego, métricas e código corretor de erros, e velocidade. Ser melhor que todos
em tudo não é um alvo realista. O alvo defensável é ser o simulador aberto com a validação mais
rastreável e o melhor para ensinar, e chegar à paridade no que os outros já tratam como básico. A
seção 5 traz a ordem proposta.

---

## 1. Como esta comparação foi feita

- **Fontes:** as páginas oficiais de cada projeto (repositório, documentação e PyPI), lidas em
  2026-11-15. A lista está na seção 7.
- **Não instalei nem rodei nenhum deles.** O que está aqui é o que as páginas declaram, não o que eu
  medi.
- **"Não consta"** quer dizer que a página lida não fala do assunto. Não quer dizer que o projeto
  não tem.
- **Maiman Studio:** parte da página do PyPI foi lida por resumo, não palavra por palavra. Os itens
  de código corretor de erros e a data da versão vêm desse resumo.
- **O simcopy** está descrito pelo estado do repositório nesta data, incluindo os defeitos de
  `docs/revisao_2026-11-14.md`.

---

## 2. Os sete projetos

| Projeto | O que é, nas palavras da página | Linguagem | Licença | Atividade |
|---|---|---|---|---|
| OptiCommPy | Simula sistemas de fibra "dos bits no transmissor às métricas no receptor" | Python, com Numba e CuPy opcional | GPL-3.0 | 2351 commits, 225 estrelas, artigo no JOSS |
| QAMpy | Cadeia de DSP para simular e equalizar sinais ópticos | Python, com pythran | GPL-3.0 | 2109 commits, 89 estrelas, "ainda em alfa" |
| Robochameleon | Estrutura de código e biblioteca de componentes para simulação e análise experimental | MATLAB | GPL-3.0 | 137 commits, 70 estrelas |
| PyOFSS | Monta um sistema de fibra a partir de módulos | Python 2 | GPL-3.0 | Última versão em 2012 |
| comnum | Simulação de enlace de sistemas digitais, sem fio e ópticos | Rust, com interface Python | MIT ou Apache-2.0 | 64 commits, 2 estrelas |
| Maiman Studio | Simulador modular de enlaces ópticos e sistemas fotônicos | Python, com CuPy opcional | Apache-2.0 | 260 commits, versão 0.17.0 de outubro de 2026 |
| OpticalCommunications | Código de apoio às aulas de Comunicações Ópticas da UFCG | Python (notebooks) | GPL-3.0 | 295 commits, 51 estrelas |

O `OpticalCommunications` não é um simulador: são os notebooks do curso da UFCG, em inglês e em
português, que usam o OptiCommPy. Entra na comparação como material de ensino.

O PyOFSS resolve a equação não linear de Schrödinger generalizada para pulsos (dispersão de ordem
superior, automodulação de fase, auto-inclinação e Raman). Não tem transmissor, receptor nem DSP.
Serve de referência para a física da fibra, não para o sistema.

---

## 3. Recurso por recurso

"sim" e "não" são o que a página declara. "?" é "não consta".

### 3.1 Transmissor e formatos

| | simcopy | OptiCommPy | QAMpy | Robochameleon | comnum | Maiman |
|---|---|---|---|---|---|---|
| OOK com detecção direta | sim | sim | ? | sim, sem detalhe | ? | sim |
| PAM4 | não | sim | ? | ? | PAM | sim |
| Laser por equações de taxa | sim, com chirp e ruído de Langevin | ? (laser com ruído de fase e RIN) | não | ? | ? | sim, com modelo térmico e Fabry-Perot multimodo |
| Formatos coerentes | BPSK, QPSK, 16QAM | M-QAM, PSK, APSK, OFDM | BPSK a QAM de ordem alta | ? | QAM, PSK, PAM, CPM | M-QAM até 256 |
| Formatação probabilística | não | sim | ? | ? | sim | ? |
| Defeitos do modulador IQ | extinção, desbalanço de fase e de ganho, skew, banda S21 | desbalanço e skew no receptor | não linearidade, resposta do DAC, ENOB | ? | desbalanço IQ | ? |
| DAC e ADC com quantização | não | sim, com ENOB e jitter | ENOB | ? | sim | ? |

### 3.2 Enlace

| | simcopy | OptiCommPy | QAMpy | Robochameleon | comnum | Maiman |
|---|---|---|---|---|---|---|
| Split-step com passo adaptativo | sim | sim | ? | pelo complemento `ssprop` | sim | sim |
| Dupla polarização (Manakov) | sim | sim | ? | ? | sim | ? |
| PMD | sim | ? | sim, como defeito | ? | sim | ? |
| Raman, dispersão de ordem superior | não | ? | ? | ? | amplificação Raman | ? |
| EDFA | ganho, ASE e saturação | simples e por equações de taxa | ? | ? | sim | ASE, saturação e transientes |
| Modelos GN e EGN | não | modelos de perturbação | ? | ? | sim | ? |
| WDM | mux e demux com filtro | sim | ? | ? | ? | sim, com nós ROADM |

### 3.3 Receptor, DSP e métricas

| | simcopy | OptiCommPy | QAMpy | Robochameleon | comnum | Maiman |
|---|---|---|---|---|---|---|
| Receptor coerente com híbrida e detecção balanceada | sim | sim | só o DSP | sim, sem detalhe | sim | sim |
| Compensação de dispersão | sim | sim | ? | ? | sim | sim |
| Equalizador adaptativo cego | não: mínimos quadrados com os símbolos conhecidos | CMA, RDE, NLMS, DD-LMS, RLS | CMA, RDE, dirigido por decisão | "blocos de DSP padrão" | sim | sim, em dupla polarização |
| Recuperação de portadora cega | não: com os símbolos conhecidos | busca cega de fase, PLL, Viterbi-Viterbi | busca cega de fase, Viterbi-Viterbi | ? | sim | sim |
| Recuperação de relógio | não | Gardner | sincronismo por pilotos | ? | sim | sim |
| Retropropagação digital | não | sim | ? | ? | sim | ? |
| Código corretor de erros | não | LDPC, Hamming | ? | ? | LDPC, polar, turbo, RS, BCH | RS, staircase, oFEC |
| Métricas | Q, BER, SER, EVM e SNR | BER, SER, SNR, Q, EVM, MI, GMI, NGMI | SER, BER, EVM, GMI | ? | BER, capacidades | Q, BER, OSNR, SNR, EVM |

### 3.4 Uso, validação e ensino

| | simcopy | OptiCommPy | QAMpy | Robochameleon | comnum | Maiman |
|---|---|---|---|---|---|---|
| Como se define uma simulação | YAML ou JSON, validado antes de rodar, ou Python | script Python | script Python | código MATLAB | campanhas em YAML, ou Python | diagrama de blocos no navegador, JSON ou Python |
| Linha de comando | validar, resolver, rodar, varrer | ? | ? | ? | sim | sim |
| Interface gráfica | não | ? | não | não | página de teste em WebAssembly | sim |
| Aceleração | não | Numba e GPU | pythran | ? | Rust, com Monte Carlo paralelo | GPU opcional |
| Testes declarados | 73, no núcleo e na fibra | suíte automática, sem número | pasta de testes | ? | 833 de unidade, 219 de validação, 293 em Python | mais de 2250 |
| Validação declarada | constelações do VPI em 125 e 140 km, e formas fechadas nos notebooks | revisão por pares no JOSS | ? | "validado em dados reais" | cada bloco contra uma referência nomeada | 161 formas fechadas, 159 a cada commit |
| Material de ensino | 9 notebooks em português | cerca de 20 notebooks, com Colab | documentação da API | manual e guia rápido | tutoriais e galeria | 48 lições |
| Publicado no PyPI | não | sim | sim | não se aplica | sim | sim |

---

## 4. Onde o simcopy está

### 4.1 O que ele já tem de diferente

1. **Comparação com um simulador comercial, com os dados no repositório.** No cenário de 4 canais
   DP-16QAM, a SNR fica de 0,03 a 0,16 dB acima da do VPI em 15 pontos, e os três termos do ajuste
   fecham em 0,25 dB (`docs/comparacao_vpi.md`). Nas páginas lidas, o OptiCommPy e o QAMpy não
   declaram comparação com simulador comercial. O comnum valida contra referências abertas
   (OptiCommPy, GNPy, Optilux). O Maiman tem uma página de comparação com OptiSystem e VPI, que não
   li.
2. **Diagnóstico da referência, não só reprodução.** A mesma análise achou o provável erro nas
   tabelas de BER do grupo e separou o piso do VPI em imagem de quadratura, vazamento de portadora
   e resto.
3. **IMDD com laser por equações de taxa e coerente no mesmo pacote e na mesma configuração.** Só o
   Maiman declara os dois.
4. **Cada número com fonte.** Agrawal, código do professor ou decisão registrada, com a dúvida
   aberta quando não há fonte (`docs/duvidas_abertas.md`).
5. **Ensino em português, na ordem das aulas.** Os notebooks montam cada cadeia à mão antes de usar
   o pacote. O `OpticalCommunications` da UFCG tem notebooks em português, apoiados no OptiCommPy.
6. **Licença MIT.** OptiCommPy, QAMpy, Robochameleon e PyOFSS são GPL-3.0. O comnum e o Maiman
   também são permissivos.
7. **Camadas com contrato verificado** (import-linter) e saída de todo estágio guardada para
   inspeção.

### 4.2 Onde ele está atrás

1. **Testes.** 73 testes, nenhum para blocos, montagem, medida, DSP ou CLI. O comnum declara 1345 e
   o Maiman mais de 2250. A revisão de 2026-11-14 listou 37 pontos, 6 deles com resultado errado sem
   aviso.
2. **DSP cego.** O receptor usa os símbolos transmitidos. OptiCommPy, QAMpy, comnum e Maiman têm
   equalização e recuperação de portadora cegas. Sem isso o simcopy não reproduz o limiar de
   convergência do VPI, nem serve para dados de laboratório.
3. **Métricas e código corretor.** Sem MI, GMI nem OSNR, e sem nenhum código. OptiCommPy, comnum e
   Maiman têm.
4. **Velocidade.** Python puro. Um ponto do cenário de 4 canais leva de 30 a 200 s com $2^{14}$
   símbolos, e os $2^{18}$ do VPI custariam 16 vezes mais. OptiCommPy tem Numba e GPU, o comnum é
   Rust.
5. **Formatos e dispositivos.** Sem PAM4, 32QAM, 64QAM, DAC e ADC, APD, Raman e inclinação da
   dispersão.
6. **Distribuição.** Não está no PyPI, não tem site de documentação e os testes novos não rodam na
   integração contínua.
7. **Interface.** Só o Maiman tem interface gráfica. A proposta de API e interface do simcopy está
   em `docs/proposta_api_interface.md`, sem implementação.

---

## 5. Como deixar o simcopy melhor

**[Opinião minha, para decisão do Gabriel.]** Três projetos puxam a régua, cada um num eixo:
o OptiCommPy em amplitude de DSP e comunidade, o Maiman em testes e interface, o comnum em
rastreabilidade da validação e velocidade. Alcançar os três em tudo é trabalho de anos. A ordem
abaixo fecha primeiro o que impede o simcopy de ser levado a sério, e depois aprofunda o que ele já
faz melhor.

### 5.1 Primeiro: confiança (o que os outros já tratam como básico)

| # | O que fazer | Por quê | Quem puxa a régua |
|---|---|---|---|
| 1 | Corrigir os 6 defeitos de resultado errado sem aviso da revisão de 2026-11-14 | Um número errado sem aviso anula a validação | todos |
| 2 | Testes dos blocos, da montagem, da medida e do DSP, cada um contra forma fechada, e rodando a cada commit | Hoje a validação do coerente está em notebooks e num script | Maiman, comnum |
| 3 | Transformar a comparação com o VPI em teste de regressão: os termos A, B e C dentro da tolerância | É o diferencial, e hoje nada impede que uma mudança o quebre | comnum |
| 4 | Publicar no PyPI, com versão e registro de mudanças | Todos os ativos estão lá | OptiCommPy |

### 5.2 Depois: paridade de função

| # | O que fazer | Por quê |
|---|---|---|
| 5 | DSP cego: equalizador adaptativo 2×2, recuperação de portadora e de relógio | É o que mais falta. Permite reproduzir o limiar de convergência do VPI e processar dados medidos |
| 6 | MI, GMI e OSNR | São as métricas que os artigos usam |
| 7 | PAM4, 32QAM e 64QAM | Estão na configuração, e a montagem recusa |
| 8 | DAC e ADC com quantização | Faz parte do cenário do VPI |
| 9 | Velocidade: FFT com vários núcleos, laço do laser de equações de taxa, e depois Numba | Sem isso os $2^{18}$ símbolos do VPI não cabem numa sessão de trabalho |

### 5.3 Em paralelo: aprofundar o que já é melhor

| # | O que fazer | Por quê |
|---|---|---|
| 10 | Estender a comparação com o VPI: 130 e 135 km, várias sementes, espectros ópticos, e o DSP cego contra o limiar de baixa potência | Nenhuma página lida mostra uma comparação desse tamanho com dados no repositório |
| 11 | Comparar bloco a bloco com o OptiCommPy, como o comnum faz | Segunda referência independente, e aberta |
| 12 | Fechar os notebooks de ensino: um por aula, com a forma fechada e o erro em painel próprio | É o que diferencia o simcopy para o grupo e para as aulas |
| 13 | Gravar a configuração e a versão em todo resultado | Reprodutibilidade, que o comnum já anuncia como princípio |

### 5.4 Por último

Código corretor de erros, interface gráfica, Raman e retropropagação digital. O OptiCommPy, o comnum
e o Maiman já cobrem bem essas partes, e nenhuma delas é necessária para os cenários do grupo hoje.

### 5.5 O que não copiar

- **Código.** OptiCommPy, QAMpy, Robochameleon e PyOFSS são GPL-3.0, e o simcopy é MIT. A regra do
  projeto para o BIFROST (GPL) é nunca importar. Vale decidir se a mesma regra se aplica a usar o
  OptiCommPy só como referência em testes (item 11).
- **Amplitude sem validação.** O QAMpy avisa que está em alfa e que nem todos os arquivos estão
  atualizados. O simcopy cresce melhor um bloco validado por vez.

---

## 6. Decisões que ficam com o Gabriel

1. O alvo: validação e ensino primeiro, como proposto, ou amplitude de recursos.
2. A ordem da seção 5, em especial se o DSP cego vem antes ou depois dos testes.
3. Usar o OptiCommPy como referência em testes, dado que é GPL.
4. Publicar no PyPI agora, com o que existe, ou só depois dos itens 1 a 3.

---

## 7. Fontes

Todas lidas em 2026-11-15.

- OptiCommPy: <https://github.com/edsonportosilva/OptiCommPy> e
  <https://opticommpy.readthedocs.io/en/latest/>
- QAMpy: <https://github.com/ChalmersPhotonicsLab/QAMpy>
- Robochameleon: <https://github.com/dtu-dsp/Robochameleon>
- PyOFSS: <https://pythonhosted.org/pyofss/overview.html>, <https://github.com/galilley/pyofss> e
  <https://pypi.org/project/pyofss/>
- comnum: <https://github.com/com-rs/comnum>
- Maiman Studio: <https://github.com/ehsun-sh/maiman-studio>, <https://maimanstudio.com/docs/> e
  <https://pypi.org/project/maiman/>
- OpticalCommunications: <https://github.com/edsonportosilva/OpticalCommunications>
