from __future__ import annotations

import re

from .model import Directive

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
