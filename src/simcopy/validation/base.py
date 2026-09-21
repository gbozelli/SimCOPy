from __future__ import annotations
from dataclasses import dataclass
from enum import Enum

class Severity(Enum):
    ERROR = "ERROR"
    WARNING = "WARNING"

@dataclass(frozen=True)
class Violation:
    validator: str
    severity: Severity
    message: str
    suggestion: str | None = None

    def __str__(self):
        s = f"{self.severity.value}  [{self.validator}]\n  {self.message}"
        return s + (f"\n  Sugestao: {self.suggestion}" if self.suggestion else "")

class ValidationFailed(RuntimeError):
    def __init__(self, violations):
        self.violations = violations
        super().__init__("\n".join(str(v) for v in violations))
