from __future__ import annotations

import re
from dataclasses import dataclass
from io import StringIO

from Bio import SeqIO

from .attributes import parse_attributes
from .errors import Diagnostic, GffParseError, Severity
from .model import Directive, Feature, GffDocument, Span

_TAXID_RE = re.compile(r"id=(\d+)")


def parse_directive(line: str) -> Directive:
    raw = line.rstrip("\r\n")
    if raw.strip() == "###":
        return Directive(raw, "resolution-boundary", None)

    content = raw[2:].strip() if raw.startswith(("##", "#!")) else raw.lstrip("#").strip()
    parts = content.split(None, 1)
    name = parts[0] if parts else ""
    rest = parts[1] if len(parts) > 1 else ""

    if name in ("gff-version", "gff-spec-version", "insdc-gff-version"):
        return Directive(raw, name, rest.strip())
    if name == "sequence-region":
        fields = rest.split()
        if len(fields) >= 3:
            try:
                return Directive(raw, "sequence-region", (fields[0], int(fields[1]), int(fields[2])))
            except ValueError:
                return Directive(raw, "sequence-region", None)
        return Directive(raw, "sequence-region", None)
    if name == "species":
        m = _TAXID_RE.search(rest)
        return Directive(raw, "species", int(m.group(1)) if m else rest.strip())
    if name == "transl_table":
        table: dict[str, int] = {}
        for item in re.split(r"[,\s]+", rest.strip()):
            if ":" in item:
                k, v = item.split(":", 1)
                table[k] = int(v)
        return Directive(raw, "transl_table", table)
    if name == "FASTA":
        return Directive(raw, "FASTA", None)
    return Directive(raw, "unknown", rest if rest else None)


@dataclass
class ParsedRow:
    id: str | None
    source: str
    type: str
    span: Span
    attributes: dict[str, list[str]]
    parent_ids: list[str]
    line_no: int


def parse_feature_line(
    line: str, line_no: int, diagnostics: list[Diagnostic]
) -> ParsedRow | None:
    cols = line.rstrip("\r\n").split("\t")
    if len(cols) != 9:
        diagnostics.append(
            Diagnostic(Severity.ERROR, line_no, "col-count", f"expected 9 columns, got {len(cols)}")
        )
        return None
    seqid, source, ftype, start_s, end_s, score_s, strand, phase_s, attr_s = cols
    try:
        start = int(start_s)
        end = int(end_s)
    except ValueError:
        diagnostics.append(
            Diagnostic(Severity.ERROR, line_no, "coord", f"non-integer start/end: {start_s!r},{end_s!r}")
        )
        return None
    score = None if score_s == "." else float(score_s)
    phase = None if phase_s == "." else int(phase_s)

    if start > end:
        diagnostics.append(
            Diagnostic(Severity.WARNING, line_no, "start-gt-end",
                       f"start>end ({start}>{end}); possible origin-spanning feature")
        )
    if not attr_s.isascii():
        diagnostics.append(
            Diagnostic(Severity.WARNING, line_no, "non-ascii", "non-ASCII characters in attributes")
        )

    attrs = parse_attributes(attr_s)
    part = None
    if attrs.get("part"):
        part = int(attrs["part"][0])
        attrs.pop("part", None)

    span = Span(seqid, start, end, strand, phase, score, part)
    fid = attrs.get("ID", [None])[0]
    parent_ids = list(attrs.get("Parent", []))
    return ParsedRow(fid, source, ftype, span, attrs, parent_ids, line_no)


def _add_row(doc: GffDocument, row: ParsedRow) -> None:
    if row.id is not None and row.id in doc.feature_index:
        feat = doc.feature_index[row.id]
        if feat.type != row.type:
            doc.diagnostics.append(
                Diagnostic(
                    Severity.ERROR, row.line_no, "id-type-mismatch",
                    f"ID {row.id!r} reused with different type ({feat.type} vs {row.type})",
                )
            )
            doc.features.append(
                Feature(row.id, row.source, row.type, [row.span], dict(row.attributes), list(row.parent_ids))
            )
            return
        if row.attributes != feat.attributes:
            doc.diagnostics.append(
                Diagnostic(
                    Severity.WARNING, row.line_no, "attr-mismatch",
                    f"row for ID {row.id!r} has attributes differing from the first row; keeping the first row's attributes",
                )
            )
        feat.spans.append(row.span)
        return
    feat = Feature(row.id, row.source, row.type, [row.span], dict(row.attributes), list(row.parent_ids))
    doc.features.append(feat)
    if row.id is not None:
        doc.feature_index[row.id] = feat


def _resolve_graph(doc: GffDocument) -> None:
    for feat in doc.features:
        for pid in feat.parent_ids:
            parent = doc.feature_index.get(pid)
            if parent is None:
                doc.diagnostics.append(
                    Diagnostic(Severity.WARNING, None, "dangling-parent",
                               f"Parent {pid!r} not found for feature {feat.id!r}")
                )
                continue
            parent.children.append(feat)
            feat.parents.append(parent)
    doc.roots = [f for f in doc.features if not f.parent_ids]


def _parse_fasta(lines: list[str]) -> dict:
    handle = StringIO("\n".join(lines))
    return {rec.id: rec.seq for rec in SeqIO.parse(handle, "fasta")}


def parse(text: str, *, strict: bool = False) -> GffDocument:
    doc = GffDocument()
    in_fasta = False
    fasta_lines: list[str] = []
    dropped_comments = 0

    for line_no, line in enumerate(text.splitlines(), start=1):
        if in_fasta:
            fasta_lines.append(line)
            continue
        if line == "":
            continue
        if line.startswith("#"):
            if line.startswith(("##", "#!")) or line.strip() == "###":
                directive = parse_directive(line)
                doc.directives.append(directive)
                if directive.kind == "FASTA":
                    in_fasta = True
            else:
                dropped_comments += 1
            continue
        row = parse_feature_line(line, line_no, doc.diagnostics)
        if row is not None:
            _add_row(doc, row)

    if fasta_lines:
        doc.fasta = _parse_fasta(fasta_lines)
    if dropped_comments:
        doc.diagnostics.append(
            Diagnostic(Severity.INFO, None, "dropped-comments", f"{dropped_comments} bare comment line(s) ignored")
        )

    _resolve_graph(doc)

    if strict:
        for d in doc.diagnostics:
            if d.severity == Severity.ERROR:
                raise GffParseError(d)
    return doc
