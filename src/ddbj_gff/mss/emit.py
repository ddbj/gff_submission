from __future__ import annotations

from .model import MssDocument, MssQualifier


def emit_ann(doc: MssDocument) -> str:
    lines: list[str] = list(doc.common_rows)
    for entry in doc.entries:
        first_row_of_entry = True
        for feat in entry.features:
            quals = feat.qualifiers or [MssQualifier("", "")]
            for i, q in enumerate(quals):
                col1 = entry.name if first_row_of_entry else ""
                col2 = feat.key if i == 0 else ""
                col3 = feat.location if i == 0 else ""
                lines.append("\t".join([col1, col2, col3, q.key, q.value]))
                first_row_of_entry = False
    return "\n".join(lines) + "\n"


def emit_fasta(seqs: dict[str, object]) -> str:
    out: list[str] = []
    for name, seq in seqs.items():
        out.append(f">{name}")
        s = str(seq)
        for i in range(0, len(s), 60):
            out.append(s[i : i + 60])
        out.append("//")
    return "\n".join(out) + "\n"
