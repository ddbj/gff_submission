from __future__ import annotations

from dataclasses import dataclass, field

from Bio.SeqFeature import CompoundLocation, FeatureLocation

from .errors import Diagnostic


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


@dataclass
class GffDocument:
    directives: list[Directive] = field(default_factory=list)
    features: list[Feature] = field(default_factory=list)
    feature_index: dict[str, Feature] = field(default_factory=dict)
    roots: list[Feature] = field(default_factory=list)
    fasta: dict | None = None
    sequences: dict | None = None
    diagnostics: list[Diagnostic] = field(default_factory=list)

    def _directive(self, kind: str) -> Directive | None:
        for d in self.directives:
            if d.kind == kind:
                return d
        return None

    @property
    def gff_version(self) -> str | None:
        d = self._directive("gff-version")
        return d.value if d else None

    @property
    def insdc_gff_version(self) -> str | None:
        d = self._directive("insdc-gff-version")
        return d.value if d else None

    @property
    def species(self) -> int | None:
        d = self._directive("species")
        return d.value if d else None

    @property
    def sequence_regions(self) -> dict[str, tuple[int, int]]:
        out: dict[str, tuple[int, int]] = {}
        for d in self.directives:
            if d.kind == "sequence-region" and d.value:
                seqid, start, end = d.value
                out[seqid] = (start, end)
        return out

    @property
    def transl_table_map(self) -> dict | None:
        d = self._directive("transl_table")
        return d.value if d else None

    def get(self, feature_id: str) -> Feature | None:
        return self.feature_index.get(feature_id)

    @staticmethod
    def _directive_key(d: "Directive"):
        v = d.value
        if isinstance(v, dict):
            v = tuple(sorted(v.items()))
        elif isinstance(v, list):
            v = tuple(v)
        return (d.kind, v)

    @staticmethod
    def _feature_key(f: "Feature"):
        return (
            f.id,
            f.type,
            f.source,
            tuple(sorted(s.sort_key() for s in f.spans)),
            tuple(sorted((k, tuple(v)) for k, v in f.attributes.items())),
            tuple(sorted(f.parent_ids)),
            tuple(sorted(c.id for c in f.children if c.id)),
            tuple(sorted(p.id for p in f.parents if p.id)),
        )

    def semantically_equals(self, other: "GffDocument") -> bool:
        """Return True if other is semantically equal (order/whitespace ignored).

        Intentional limitations (Phase 1 round-trip oracle):
        - Compares the SET of directive keys and feature keys, so exact-duplicate
          directives and exact-duplicate ID-less features collapse and are not counted.
        - Does NOT compare the FASTA payload (peptide sequences); test FASTA fidelity separately.
        - Attribute values are compared as an ORDERED tuple (stricter than a multiset);
          this is safe because the writer preserves value order on round-trip.
        """
        if {self._directive_key(d) for d in self.directives} != {
            other._directive_key(d) for d in other.directives
        }:
            return False
        return {self._feature_key(f) for f in self.features} == {
            other._feature_key(f) for f in other.features
        }
