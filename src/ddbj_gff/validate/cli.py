from __future__ import annotations

import argparse
import sys

from .. import parse
from ..errors import Severity
from .severities import resolve_level
from .validate import validate


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="ddbj_gff.validate",
                                 description="Validate a GFF3 file against the INSDC profile")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--severity", action="append", default=[],
                    metavar="CODE=LEVEL", help="override a rule severity (error/warning/info/off)")
    args = ap.parse_args(argv)

    overrides: dict[str, str] = {}
    for item in args.severity:
        if "=" not in item:
            ap.error(f"--severity expects CODE=LEVEL, got {item!r}")
        code, level = item.split("=", 1)
        try:
            resolve_level(level.strip())
        except ValueError as e:
            ap.error(str(e))
        overrides[code.strip()] = level.strip()

    with open(args.gff, encoding="ascii", errors="replace") as fh:
        doc = parse(fh.read())

    diags = list(doc.diagnostics) + validate(doc, severity_overrides=overrides)

    counts: dict[str, int] = {}
    for d in diags:
        counts[d.severity.value] = counts.get(d.severity.value, 0) + 1
        sys.stderr.write(f"{d.severity.value}\t{d.code}\t{d.message}\n")
    sys.stderr.write("summary: "
                     + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())) + "\n")
    return 1 if counts.get(Severity.ERROR.value) else 0
