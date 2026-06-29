from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Severity(Enum):
    ERROR = "ERROR"
    WARNING = "WARNING"
    INFO = "INFO"


@dataclass(frozen=True)
class Diagnostic:
    severity: Severity
    line_no: int | None
    code: str
    message: str


class GffParseError(Exception):
    def __init__(self, diagnostic: Diagnostic):
        self.diagnostic = diagnostic
        line = "?" if diagnostic.line_no is None else diagnostic.line_no
        super().__init__(
            f"{diagnostic.severity.value} (line {line}) "
            f"[{diagnostic.code}] {diagnostic.message}"
        )
