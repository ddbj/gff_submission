from __future__ import annotations

import argparse
import sys

from .. import parse
from ..writer import write
from ..normalize.report import NormalizationReport
from .context import RepairContext
from .registry import list_operations, REGISTRY
from .report import candidates_to_json, render_candidates
from .driver import run_detect, run_apply, load_sequences, DEFAULT_ORDER


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="ddbj_gff.repair",
                                 description="Apply modular GFF curation operations")
    ap.add_argument("--list", action="store_true", help="list available operations")
    ap.add_argument("--gff")
    ap.add_argument("--fasta")
    ap.add_argument("--transl-table", type=int, default=1, dest="transl_table")
    ap.add_argument("--detect", action="store_true", help="preview candidates; write nothing")
    ap.add_argument("--json", action="store_true", help="detect output as JSON")
    ap.add_argument("--only", help="comma-separated operation names (detect)")
    ap.add_argument("--apply", help="comma-separated operation names, or 'all'")
    ap.add_argument("--out")
    ap.add_argument("--report")
    args = ap.parse_args(argv)

    if args.list:
        for op in list_operations():
            seq = " (needs FASTA)" if op.requires_sequence else ""
            sys.stdout.write(f"{op.name}{seq}: {op.summary}\n")
        return 0

    if not args.gff:
        ap.error("--gff is required unless --list")

    with open(args.gff, encoding="ascii", errors="replace") as fh:
        doc = parse(fh.read())
    sequences = load_sequences(args.fasta) if args.fasta else None
    ctx = RepairContext(sequences=sequences, transl_table=args.transl_table)

    if args.detect:
        names = args.only.split(",") if args.only else None
        cands = run_detect(doc, ctx, names)
        sys.stdout.write(candidates_to_json(cands) + "\n" if args.json
                         else render_candidates(cands))
        return 0

    if args.apply:
        if args.apply == "all":
            names = [n for n in DEFAULT_ORDER if n in REGISTRY]
        else:
            names = args.apply.split(",")
        changes = run_apply(doc, ctx, names)
        out_text = write(doc)
        if args.out:
            with open(args.out, "w", encoding="ascii") as fh:
                fh.write(out_text)
        else:
            sys.stdout.write(out_text)
        report = NormalizationReport(applied=changes, unresolved=[]).render()
        if args.report:
            with open(args.report, "w", encoding="ascii") as fh:
                fh.write(report)
        else:
            sys.stderr.write(report)
        return 0

    ap.error("nothing to do: pass --list, --detect, or --apply")
    return 2
