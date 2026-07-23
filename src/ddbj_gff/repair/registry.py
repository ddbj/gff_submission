from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from ..model import GffDocument
from ..normalize.report import Change
from .context import RepairContext
from .report import Candidate


@dataclass
class Operation:
    name: str
    summary: str
    requires_sequence: bool
    detect: Callable[[GffDocument, RepairContext], list[Candidate]]
    apply: Callable[[GffDocument, RepairContext, "list[Candidate] | None"], list[Change]]


REGISTRY: dict[str, Operation] = {}


def register(op: Operation) -> Operation:
    REGISTRY[op.name] = op
    return op


def get_operation(name: str) -> Operation:
    try:
        return REGISTRY[name]
    except KeyError:
        raise KeyError(f"unknown repair operation {name!r}; "
                       f"available: {sorted(REGISTRY)}") from None


def list_operations() -> list[Operation]:
    return list(REGISTRY.values())
