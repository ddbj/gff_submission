# ddbj_gff.normalize (reference)

`$PY -m ddbj_gff.normalize --gff GFF [--fasta FA] [--config C.toml] [--taxid N] [--transl-table N] [--insdc-gff-version V] [--out OUT] [--report R]`

Canonicalizes a GFF3 toward the INSDC profile (adds directives; coerces coding
`transcript`→`mRNA`; wraps/reparents CDS under mRNA; SO term mapping; trans-splicing
location; recoded/anticodon children). Pass `--fasta` so `##sequence-region` lengths are
exact (otherwise approximated from max feature end). `--report` writes an applied /
needs-attention summary — review `needs-manual` and `unmapped-type` lines. Without `--out`
the normalized GFF goes to stdout; without `--report` the report goes to stderr.
