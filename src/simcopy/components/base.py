"""Contrato de bloco (arquitetura, secao 5.3 e anexo A.1; docs/mapa_arquitetura.md, 4.5).

Todo bloco e um dataclass congelado com a configuracao do seu dispositivo no construtor, e o
processo e o __call__(x, ctx). Os blocos novos tambem tem:
- `id`: chave do ctx.rng_for (sementes independentes por bloco) e rotulo do estagio no resultado;
- `output_unit` (atributo de classe): unidade da saida, para os graficos (simcopy.plotting).
Blocos nao plotam, nao imprimem e nao mutam a entrada.
"""
from __future__ import annotations
from typing import Any, Protocol, runtime_checkable
from ..core.context import SimContext

@runtime_checkable
class Block(Protocol):
    def __call__(self, x: Any, ctx: SimContext) -> Any: ...
