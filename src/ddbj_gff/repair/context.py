from __future__ import annotations

from dataclasses import dataclass


@dataclass
class RepairContext:
    sequences: dict | None = None   # seqid -> Bio.Seq.Seq (nucleotide)
    transl_table: int = 1           # default table when a CDS omits transl_table
