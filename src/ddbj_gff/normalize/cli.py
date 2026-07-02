from __future__ import annotations

import argparse
import sys

from .. import parse
from ..io import open_text
from ..writer import write
from .config import NormalizeConfig, load_normalize_config
from .normalize import normalize


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="ddbj_gff.normalize",
                                 description="Normalize a GFF3 toward the INSDC profile")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta")
    ap.add_argument("--config")
    ap.add_argument("--taxid", type=int)
    ap.add_argument("--transl-table", type=int, dest="transl_table")
    ap.add_argument("--insdc-gff-version", dest="insdc_gff_version")
    ap.add_argument("--out")
    ap.add_argument("--report")
    args = ap.parse_args(argv)

    cfg = load_normalize_config(args.config) if args.config else NormalizeConfig()
    if args.taxid is not None:
        cfg.taxid = args.taxid
    if args.transl_table is not None:
        cfg.transl_table = args.transl_table
    if args.insdc_gff_version is not None:
        cfg.insdc_gff_version = args.insdc_gff_version

    seq_lengths = None
    if args.fasta:
        from Bio import SeqIO
        with open_text(args.fasta) as fh:
            seq_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(fh, "fasta")}

    with open(args.gff, encoding="ascii", errors="replace") as fh:
        doc = parse(fh.read())

    norm, report = normalize(doc, seq_lengths=seq_lengths, config=cfg)

    out_text = write(norm)
    if args.out:
        with open(args.out, "w", encoding="ascii") as fh:
            fh.write(out_text)
    else:
        sys.stdout.write(out_text)

    report_text = report.render()
    if args.report:
        with open(args.report, "w", encoding="ascii") as fh:
            fh.write(report_text)
    else:
        sys.stderr.write(report_text)

    return 0
