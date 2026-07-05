from __future__ import annotations

from ..model import Span

_STRAND = {1: "+", -1: "-", 0: "?", None: "."}

# INSDC qualifier -> canonical GFF attribute key
_QUAL_MAP = {
    "gene": "gene", "locus_tag": "locus_tag", "product": "product",
    "note": "Note", "protein_id": "protein_id", "gene_synonym": "gene_synonym",
    "transl_table": "transl_table", "db_xref": "Dbxref", "pseudo": "pseudo",
    "ncRNA_class": "ncRNA_class",
}
_DROP_QUALS = {"translation", "codon_start"}


def bio_location_to_spans(location, seqid, *, is_cds, codon_start=1):
    """BioPython location -> 1-based Spans in biological (5'->3') order.
    parts come from BioPython in transcription order. For CDS, per-segment phase
    is derived from codon_start; else phase is None."""
    parts = list(location.parts)
    spans = []
    phase = (codon_start - 1) if is_cds else None
    for p in parts:
        start = int(p.start) + 1
        end = int(p.end)
        strand = _STRAND.get(p.strand, ".")
        spans.append(Span(seqid, start, end, strand, phase=phase))
        if is_cds:
            seg_len = end - start + 1
            phase = (3 - ((seg_len - phase) % 3)) % 3
    return spans


def qualifiers_to_attrs(feature) -> dict:
    attrs: dict[str, list[str]] = {}
    for k, vals in feature.qualifiers.items():
        if k in _DROP_QUALS:
            continue
        gk = _QUAL_MAP.get(k, k)
        attrs[gk] = list(vals)
    return attrs
