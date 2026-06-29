from __future__ import annotations

from .attributes import serialize_attributes
from .model import Feature, GffDocument, Span


def _fmt_score(score: float | None) -> str:
    return "." if score is None else repr(score)


def _format_row(feat: Feature, span: Span) -> str:
    attrs = dict(feat.attributes)
    if span.part is not None:
        attrs["part"] = [str(span.part)]
    col9 = serialize_attributes(attrs)
    phase = "." if span.phase is None else str(span.phase)
    return "\t".join(
        [
            span.seqid,
            feat.source,
            feat.type,
            str(span.start),
            str(span.end),
            _fmt_score(span.score),
            span.strand,
            phase,
            col9,
        ]
    )


def _format_fasta(fasta: dict) -> str:
    out: list[str] = []
    for fid, seq in fasta.items():
        out.append(f">{fid}")
        s = str(seq)
        for i in range(0, len(s), 60):
            out.append(s[i : i + 60])
    return "\n".join(out) + "\n"


def _canonical_sort(features: list[Feature]) -> list[Feature]:
    def key(f: Feature):
        s = f.ordered_spans()
        first = s[0] if s else Span("", 0, 0)
        return (first.seqid, first.start, first.end)

    return sorted(features, key=key)


def write(doc: GffDocument, *, canonical_sort: bool = False) -> str:
    lines: list[str] = []
    for d in doc.directives:
        if d.kind == "FASTA":
            continue
        lines.append(d.raw)

    features = _canonical_sort(doc.features) if canonical_sort else doc.features
    for feat in features:
        for span in feat.ordered_spans():
            lines.append(_format_row(feat, span))

    text = "\n".join(lines)
    if text:
        text += "\n"
    if doc.fasta:
        text += "##FASTA\n" + _format_fasta(doc.fasta)
    return text
