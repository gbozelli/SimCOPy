Todas as entradas precisam ser vetoriais para PMD

Fibra
Todos isolados
1 - fibra propagação de um pulso gaussiano/rrc
2 - atenuação (1 pulso, como ele varia ao longo da distancia)
3 - dispersão cromática (1 pulso, depois 2, mesmo lambda, propagação, ver fase em termos do tempo, e depois em frequência) - espectro em função da frequencia
4 - pmd -> esfera de poincaré (ver como o estado de plarização vai mudando)
5 - efeitos não lineares , só kerr por enquanto (n simulamos biriwan e raman pois n afeta muito) (1 pullso, depois 1 pulso com dobro de amplitude, tirar phi(t), sem considerar dispersão, e depois considerando ela. fazer espectograma na saída, F x D, criar as componentes de four =wave mixing, e depois isso vai se criando em cascata, e depois n se gera, pois não há sobreposição dos pulsos)
6 - teste de todos os efeitos juntos

Notebook pra isso (revisão dos conceitos, depois código, escrever a equação, e zerar os termos conforme formos testando)

Lasers
- Também tem as equações descritas matematicamente, e os modelos de rate equations (Mod. direta) e o modelo quase-monocromático (Mod. externa), para diferentes correntes de injeção, existem diferentes fluxos de fotons. gráfico N_c x t para um Ibias, embaixo um fluxo de fótons phi x t
 e embaixo Pout X I (), e um de saida pela entrada (entrada um pulso quadrado I(t), sai uma potência(t))

Moduladores
- Podem ser descritos por monocromático pois se usam ondas continuas (geralmente não há modulação no laser, exceto para sistemas específicos).

Amplificadores (EDFA, não usaremos Raman, )
- Como se agrega ruido no EDFA, como se calcula snr na saída e na entrada

Fotodetectores
- MUX E DEMUX

A partir daqui podemos fazer TX, enlace e RX

Passivos (Não se pode testar sozinhos, pois precisa da cadeia de transmissão tx/rx)

IMDD < 5w
Criar sequência pseudoaleatória, converter para sinal elétrico (n se usa Nyquist, e sim um filtro de Bessel)
Laser - Fibra - Fotodetector
Teste: Colocar o laser e ver a saída óptica, entrada elétrica Colocar Laser e fibra diretos, sem a fotodetector, e depois laser, fibra e fotodetector, e monitorar a saída óptica em todos.
Laser de onda contínua - MZM - fibra - fotodetector
Colocar duas fibras ao inves de uma (SMF, DCF)
Colocar um aplificador depois da fibra DCF
Múltiplos spans
Múltiplos canais (mux e demux, enlace compartilhado)

ja tinha modulação no IMDD
Coerente digital
Inicialmente BPSK
Laser - MZM (modulação em amplitude, phase = 0) - 
Laser - MZM (modulação QAM, I(t) e Q(t)), se i e q são OOK, é constante, não há recuperação de fase) - photodetector
Laser - MZM (modulação QAM, I(t) e Q(t)), se i e q são OOK, é constante, não há recuperação de fase) mas no final de mzm tem um laser - photodetector
Laser - MZM (modulação QAM, I(t) e Q(t)), se i e q são OOK, é constante, não há recuperação de fase) mas no final de mzm tem um laser, -n balanceador - photodetector (ai vai conseguir recuperar info)
