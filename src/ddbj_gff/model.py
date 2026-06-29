from __future__ import annotations

from dataclasses import dataclass, field

from Bio.SeqFeature import CompoundLocation, FeatureLocation


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


@dataclass
class Feature:
    id: str | None
    source: str
    type: str
    spans: list[Span] = field(default_factory=list)
    attributes: dict[str, list[str]] = field(default_factory=dict)
    parent_ids: list[str] = field(default_factory=list)
    children: list["Feature"] = field(default_factory=list, repr=False)
    parents: list["Feature"] = field(default_factory=list, repr=False)

    _STRAND_MAP = {"+": 1, "-": -1}

    def _first(self, key: str) -> str | None:
        vals = self.attributes.get(key)
        return vals[0] if vals else None

    def _int(self, key: str) -> int | None:
        v = self._first(key)
        return int(v) if v is not None and v != "" else None

    @property
    def name(self) -> str | None:
        return self._first("Name")

    @property
    def locus_tag(self) -> str | None:
        return self._first("locus_tag")

    @property
    def gene(self) -> str | None:
        return self._first("gene")

    @property
    def product(self) -> str | None:
        return self._first("product")

    @property
    def protein_id(self) -> str | None:
        return self._first("protein_id")

    @property
    def transl_table(self) -> int | None:
        return self._int("transl_table")

    @property
    def number(self) -> int | None:
        return self._int("number")

    @property
    def exon_number(self) -> int | None:
        return self._int("exon_number")

    @property
    def note(self) -> list[str]:
        return self.attributes.get("Note", [])

    @property
    def dbxref(self) -> list[str]:
        return self.attributes.get("Dbxref", [])

    @property
    def is_circular(self) -> bool:
        return self._first("Is_circular") == "true"

    @property
    def is_ordered(self) -> bool:
        return self._first("is_ordered") == "true"

    @property
    def is_trans_spliced(self) -> bool:
        v = self._first("exception")
        return v is not None and v.replace("_", "-") == "trans-splicing"

    def ordered_spans(self) -> list[Span]:
        spans = list(self.spans)
        if not spans:
            return spans
        if any(s.part is not None for s in spans):
            return sorted(spans, key=lambda s: (s.part is None, s.part if s.part is not None else 0))
        if all(s.strand == "-" for s in spans):
            return sorted(spans, key=lambda s: s.start, reverse=True)
        return sorted(spans, key=lambda s: s.start)

    @property
    def codon_start(self) -> int | None:
        if self.type != "CDS" or not self.spans:
            return None
        phase = self.ordered_spans()[0].phase
        return None if phase is None else phase + 1

    def to_biopython_location(self):
        parts = []
        for s in self.ordered_spans():
            strand = self._STRAND_MAP.get(s.strand, 0)
            parts.append(FeatureLocation(s.start - 1, s.end, strand=strand))
        if len(parts) == 1:
            return parts[0]
        return CompoundLocation(parts)
