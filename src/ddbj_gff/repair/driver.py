from __future__ import annotations

from Bio import SeqIO

from ..io import open_text
from .context import RepairContext
from .registry import get_operation, list_operations
from .report import Candidate

DEFAULT_ORDER = ["internal-stop-to-misc", "utr-absent-to-partial-mrna",
                 "missing-start-stop-to-partial-cds"]


def load_sequences(path: str) -> dict:
    with open_text(path) as fh:
        return {rec.id: rec.seq for rec in SeqIO.parse(fh, "fasta")}


def run_detect(doc, ctx: RepairContext, names=None) -> list[Candidate]:
    ops = [get_operation(n) for n in names] if names else list_operations()
    out: list[Candidate] = []
    for op in ops:
        out.extend(op.detect(doc, ctx))
    return out


def run_apply(doc, ctx: RepairContext, names) -> list:
    changes: list = []
    for n in names:
        changes.extend(get_operation(n).apply(doc, ctx, None))
    return changes
