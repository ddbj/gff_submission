# ddbj_gff.repair (reference)

Two-phase curation: `--detect` (non-destructive preview) then `--apply` (writes a new GFF).

- List operations: `$PY -m ddbj_gff.repair --list`
- Preview: `$PY -m ddbj_gff.repair --gff G.gff --fasta FA --detect --json [--only OPS]`
- Apply: `$PY -m ddbj_gff.repair --gff G.gff --fasta FA --apply OPS|all --out repaired.gff [--report R]`

Operations:
- `internal-stop-to-misc` — a CDS whose translation has an internal stop → retype the CDS to
  `misc_feature` + Note (gene/mRNA kept). Sequence-based (needs `--fasta`).
- `utr-absent-to-partial-mrna` — mRNA missing a UTR on a genomic end → mark that end partial.
  Structural (no FASTA needed).
- `missing-start-stop-to-partial-cds` — CDS lacking a start/stop codon → mark the end partial.
  Sequence-based (needs `--fasta`).

Partiality is written as INSDC `partial=true` + `start_range`/`end_range`. Sequence-based ops
silently find nothing without `--fasta`. Default `--apply all` order runs
`internal-stop-to-misc` first (so retyped features are excluded from the partial-CDS op).
`--detect --json` prints machine-readable candidates for an agent to select from.
