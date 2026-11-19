# Potência de lançamento: o PowerControl basta? Proposta para discussão (2026-11-18)

Nasceu de uma pergunta feita olhando o `plot_chain` do notebook 08: o `PowerControl` é a melhor
forma de impor a potência de lançamento numa simulação realista?

**Resposta em uma frase:** é uma boa escolha como atenuador ideal (medir a BER contra a potência
recebida) e uma escolha ruim como amplificador, porque amplifica sem ruído. A proposta é manter o
bloco, separar bem os dois papéis e avisar quando ele estiver amplificando.

Nada aqui foi implementado. A seção 5 lista o que precisa ser decidido.

---

## 1. O que existe hoje

**[Fato: código em `src/simcopy`]**

- **`PowerControl`** (`components/optical/power.py`): multiplica o campo por
  `sqrt(potência_alvo / potência_atual)`, com a potência média somando X e Y. Não tem ruído, perda
  de inserção nem limite. É montado no fim do transmissor, quando o canal tem `launch_power`
  (`pipelines/imdd.py:75` e `pipelines/coherent.py:95`).
- **`Amplifier`** (`components/optical/amplifier.py`): um EDFA com ganho fixo (com saturação) ou
  controlado por potência (`target_output_power`), e com ASE (figura de ruído, Eq. 7.2.15). Fica em
  `link.booster` (antes das fibras) e em `span.amplifier` (depois de cada fibra).
- **Os dois juntos.** O `launch_power` é a potência do canal na saída do transmissor. Se há booster,
  o que entra na fibra é a saída dele: a regra do limiar de Brillouin
  (`validation/rules.py:200`) já trata assim.

**[Fato medido, 2026-11-18, rodando os exemplos]** O que cada exemplo faz com a potência:

| Exemplo | O que o `PowerControl` faz | Observação |
|---|---|---|
| `coerente_qpsk.yaml`, `coerente_bpsk.yaml` | atenua 29,8 dB (−11,2 para −41 dBm) | sem fibra: é a potência na entrada do receptor |
| `coerente_dp16qam.yaml`, `coerente_dp64qam.yaml` | impõe −20 dBm | o YAML diz "potência na entrada do receptor (não há enlace)" |
| `imdd_wdm_4ch.yaml` | **amplifica** 3 dB (−3 para 0 dBm), sem ruído | depois há um amplificador de ganho fixo de 22,4 dB na DCF |
| `imdd_ook_mzm.yaml`, `imdd_ook_dml.yaml` | não existe (sem `launch_power`) | o lançamento é a saída do modulador ou do laser (+7 dBm) |
| `400zr_4ch_75ghz.yaml` | não existe (sem `launch_power`) | booster de 13 dBm controlado por potência, com NF de 5 dB (não rodei este) |

O último é o caso do VPI: o `AmpSysOpt` do documento do VPI é um amplificador `PowerControlled`
com NF de 5 dB que define a potência de lançamento (`docs/amplificadores.md` §4 e
`docs/arquitetura.md` §14).

## 2. Onde o PowerControl é bom

Atenuar é passivo e sem ruído: tira potência e mais nada. É o que se faz na bancada com um
atenuador variável antes do receptor para medir a BER contra a potência recebida, como no
notebook 08. Nesse uso o `PowerControl` e o equipamento real fazem a mesma coisa. Também serve
para comparar dois casos com a mesma potência e isolar um efeito.

## 3. Onde ele esconde um defeito

1. **Amplifica sem ruído.** Se o alvo for maior que a potência que sai do transmissor, o
   `PowerControl` age como um amplificador perfeito. Na vida real esse ganho vem de um EDFA, que
   acrescenta ASE. **[Dedução]** A OSNR (em 0,1 nm) na saída de um booster com figura de ruído NF é
   aproximadamente `58 + P_entrada[dBm] − NF[dB]`. Conferi com o kernel `amplifier.py`
   (NF = 5 dB, saída de 0 dBm), e a diferença é menor que 0,1 dB:

   | P na entrada do booster | ganho necessário | OSNR na saída |
   |---|---|---|
   | −10 dBm | 10 dB | 43 dB |
   | −20 dBm | 20 dB | 33 dB |
   | −30 dBm | 30 dB | 23 dB |
   | −40 dBm | 40 dB | 13 dB |

   Com ganhos até uns 10 dB (o +3 dB do `imdd_wdm_4ch`), a OSNR fica muito acima da do enlace e
   ignorar o ruído quase não muda nada. Com 30 dB ou mais, o ruído que falta é o que domina.
2. **Anula parâmetros do transmissor.** Com `launch_power`, a potência do laser e a perda de
   inserção do MZM deixam de importar, porque a saída é renormalizada (a decisão M9 pediu isso;
   `docs/sistema_coerente.md` diz que o `PowerControl` "repõe" os 6 dB do MZM). Quem muda o
   `laser.power` e não vê efeito pode achar que é um erro.
3. **O nome engana.** `launch_power` é a potência do canal na saída do transmissor, mas, com
   booster, a potência lançada na fibra é a saída do booster. E é por canal (decisão da dúvida 3),
   enquanto o mux soma os canais (no `imdd_wdm_4ch`, 0 dBm por canal vira 6,02 dBm).
4. **Não limita nada.** Um EDFA real satura. O `PowerControl` aceita qualquer alvo sem aviso, e o
   `Amplifier` só tem `saturation_power` no modo de ganho fixo. A única conferência é o aviso de
   Brillouin do IMDD (acima de 5 mW por canal em mais de 10 km).
5. **Fixa a média da janela simulada.** **[Dedução]** Numa janela curta, a potência média real
   varia um pouco com os bits. É um efeito pequeno, e a especificação em dBm médios pede mesmo uma
   média fixa.

## 4. Opções

- **A. Deixar como está.** Simples, e os exemplos de back-to-back continuam certos. Mantém o risco do
  item 1 da seção 3, sem aviso.
- **B. Separar os papéis.** O `PowerControl` só atenua: se o alvo pedir ganho, dá erro e aponta para
  o booster. Quem precisa amplificar usa o `Amplifier`, com ASE. A cadeia realista fica: laser,
  modulador (com perda), atenuador opcional, booster, fibra. É o que o exemplo 400ZR já faz.
- **C. Aviso na validação.** O `PowerControl` continua amplificando, mas uma regra avisa quando o
  ganho passa de um limite (por exemplo 10 dB) e diz a OSNR que o EDFA real daria, pela conta
  `58 + P − NF`. Barato e não quebra nenhum exemplo.

**Recomendação.** C agora, e B se o grupo quiser que o padrão seja realista. Para simulações de
enlace, o lançamento deve vir do booster, como no VPI. O `PowerControl` fica para o back-to-back e
para as comparações.

## 5. O que decidir

1. **Papel do `PowerControl`.** Só atenuador (erro se precisar de ganho, opção B), ou como hoje com
   aviso (opção C)?
2. **Limite do aviso.** Quantos dB de ganho sem ruído toleramos antes de avisar? Ou o aviso deve
   sair pela OSNR equivalente (por exemplo, abaixo de 30 dB)?
3. **Nome.** `launch_power` passa a se chamar algo como `tx_output_power`, e o lançamento fica com
   `link.booster.output_power`? Ou mantemos o nome e só documentamos?
4. **Limite do booster.** O modo de potência alvo deve respeitar uma potência de saturação e um
   ganho máximo? O VPI tem `SaturationPower = 10e-3` W e `GainMax = 100` dB no `AmpSysOpt`.
5. **O exemplo `imdd_wdm_4ch`.** Troca o `PowerControl` de +3 dB por um booster, ou mantemos?

## 6. Como testar (insumo para a issue, que fica com o grupo)

**[Proposta]**
- O `PowerControl` com alvo menor que a entrada: a potência final é a do alvo e a forma do sinal não
  muda.
- Com alvo maior que a entrada: aparece o aviso (ou o erro, se for a opção B).
- O `Amplifier` de potência alvo: a OSNR na saída bate com `58 + P − NF` dentro de 0,1 dB.
- Com `launch_power` definido, mudar o `laser.power` não muda a saída (é o comportamento atual, e o
  teste o documenta).
