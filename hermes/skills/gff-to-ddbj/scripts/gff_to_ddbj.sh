#!/usr/bin/env bash
# gff_to_ddbj.sh — happy-path driver: normalize -> validate -> repair(apply all) -> gff2mss.
# Env-specific paths come from GFF_TO_DDBJ_ENV_BIN (or --env-bin); default = NIG/DDBJ cluster.
set -euo pipefail

ENV_BIN="${GFF_TO_DDBJ_ENV_BIN:-/lustre9/open/home/yt/micromamba/envs/mss_tools/bin}"
GFF="" FASTA="" MSSCFG="" COMMON="" OUTP="" WORK="."

usage() {
  echo "usage: $0 --gff IN.gff --fasta IN.fa --mss-config MSS.toml \\"
  echo "          --common COMMON.json --out-prefix DIR/NAME [--workdir DIR] [--env-bin DIR]"
  exit 2
}

while [ $# -gt 0 ]; do
  case "$1" in
    --gff)         GFF="$2"; shift 2;;
    --fasta)       FASTA="$2"; shift 2;;
    --mss-config)  MSSCFG="$2"; shift 2;;
    --common)      COMMON="$2"; shift 2;;
    --out-prefix)  OUTP="$2"; shift 2;;
    --workdir)     WORK="$2"; shift 2;;
    --env-bin)     ENV_BIN="$2"; shift 2;;
    -h|--help)     usage;;
    *) echo "unknown argument: $1" >&2; usage;;
  esac
done
[ -n "$GFF" ] && [ -n "$FASTA" ] && [ -n "$MSSCFG" ] && [ -n "$COMMON" ] && [ -n "$OUTP" ] || usage

PY="$ENV_BIN/python"
GFF2MSS="$ENV_BIN/gff2mss"
mkdir -p "$WORK" "$(dirname "$OUTP")"

echo "[1/4] normalize"
"$PY" -m ddbj_gff.normalize --gff "$GFF" --fasta "$FASTA" \
      --out "$WORK/norm.gff" --report "$WORK/normalize.txt"

echo "[2/4] validate (detect-only; review any ERR: lines)"
"$PY" -m ddbj_gff.validate --gff "$WORK/norm.gff" || \
  echo "  (validate reported issues — see output above; continuing happy path)"

echo "[3/4] repair (apply all)"
"$PY" -m ddbj_gff.repair --gff "$WORK/norm.gff" --fasta "$FASTA" \
      --apply all --out "$WORK/repaired.gff" --report "$WORK/repair.txt"

echo "[4/4] gff2mss"
"$GFF2MSS" --gff "$WORK/repaired.gff" --fasta "$FASTA" \
      --mss-config "$MSSCFG" --common "$COMMON" --out "$OUTP"

echo "done: ${OUTP}.ann / ${OUTP}.fasta"
echo "next: validate the MSS with ddbj-validator (see references/validator.md) —"
echo "      large genomes in the background, always with -f in a non-TTY shell."
