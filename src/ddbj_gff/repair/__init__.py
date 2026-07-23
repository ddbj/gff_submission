"""GFF repair / curation layer.

Modular, individually-invokable GFF->GFF operations. Each Operation is a
two-phase detect(doc, ctx)->list[Candidate] (non-destructive) + apply(doc, ctx,
selection)->list[Change] (mutating) unit registered in REGISTRY. Add a new
operation by writing its detect/apply and calling register() in operations.py.
"""
from __future__ import annotations

from .context import RepairContext
from .registry import Operation, REGISTRY, register, get_operation, list_operations
from .report import Candidate, candidates_to_json, render_candidates

__all__ = ["RepairContext", "Operation", "REGISTRY", "register", "get_operation",
           "list_operations", "Candidate", "candidates_to_json", "render_candidates"]

from . import operations as _operations  # noqa: E402,F401  (populates REGISTRY)
