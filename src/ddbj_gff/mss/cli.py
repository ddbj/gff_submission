from __future__ import annotations

import argparse
import sys

from Bio import SeqIO

from .. import parse
from .config import load_common, load_config
from .convert import convert
from .emit import emit_ann, emit_fasta
from ..errors import Severity


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="ddbj_gff.mss",
                                 description="Convert GFF3 + FASTA to DDBJ MSS (.ann + .fasta)")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--common", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args(argv)

    with open(args.gff, encoding="ascii", errors="replace") as fh:
        doc = parse(fh.read())
    seqs = {rec.id: rec.seq for rec in SeqIO.parse(args.fasta, "fasta")}
    cfg, cfg_diags = load_config(args.config)
    common_rows = load_common(args.common)

    mss, diags = convert(doc, seqs, cfg, common_rows, strict=args.strict)
    diags = cfg_diags + diags

    with open(f"{args.out}.ann", "w", encoding="ascii") as fh:
        fh.write(emit_ann(mss))
    with open(f"{args.out}.fasta", "w", encoding="ascii") as fh:
        fh.write(emit_fasta(seqs))

    counts: dict[str, int] = {}
    for d in diags:
        counts[d.severity.value] = counts.get(d.severity.value, 0) + 1
    if counts:
        sys.stderr.write("diagnostics: "
                         + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())) + "\n")
    return 1 if counts.get(Severity.ERROR.value) else 0
