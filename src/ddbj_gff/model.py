from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Span:
    seqid: str
    start: int  # 1-based inclusive
    end: int
    strand: str = "."  # one of: + - . ?
    phase: int | None = None  # 0/1/2 for CDS
    score: float | None = None
    part: int | None = None

    def sort_key(self) -> tuple:
        return (
            self.seqid,
            self.start,
            self.end,
            self.strand,
            -1 if self.phase is None else self.phase,
            float("-inf") if self.score is None else self.score,
            -1 if self.part is None else self.part,
        )


@dataclass
class Directive:
    raw: str
    kind: str
    value: object = None
