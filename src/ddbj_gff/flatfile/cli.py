from __future__ import annotations

import argparse

from Bio import SeqIO

from .convert import flatfile_to_gff
from ..writer import write


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="flatfile2gff",
                                 description="DDBJ flatfile -> canonical INSDC GFF3")
    ap.add_argument("--in", dest="infile", required=True, help="input DDBJ flatfile (.gbk)")
    ap.add_argument("--out", dest="outfile", required=True, help="output GFF3")
    args = ap.parse_args(argv)
    rec = SeqIO.read(args.infile, "genbank")
    doc = flatfile_to_gff(rec)
    with open(args.outfile, "w", encoding="utf-8") as fh:
        fh.write(write(doc))
    print(f"[flatfile2gff] -> {args.outfile}")


if __name__ == "__main__":
    main()
