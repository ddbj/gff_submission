from __future__ import annotations

import re
from dataclasses import dataclass

from .attributes import parse_attributes
from .errors import Diagnostic, Severity
from .model import Directive, Span

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
            return Directive(raw, "sequence-region", (fields[0], int(fields[1]), int(fields[2])))
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
