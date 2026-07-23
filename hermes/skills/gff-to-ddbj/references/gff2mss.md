# gff2mss (reference)

`$GFF2MSS --gff GFF --fasta FA --mss-config MSS.toml --common COMMON.json [--sequence-roles roles.tsv] [--submission-category CAT] [--locus-tag-start N] --out OUTPREFIX`

Writes `OUTPREFIX.ann` and `OUTPREFIX.fasta`. `--mss-config` and `--common` are REQUIRED.

Minimal `MSS.toml`:
```toml
[source]
organism = "Genus species"
mol_type = "genomic DNA"

[locus_tag]
prefix = "ABCD"   # official BioSample locus_tag prefix
width = 6
start = 10
step = 10

[cds]
transl_table = 1

[product]
default = "hypothetical protein"
map = "product_map.tsv"   # optional 2-col TSV: id<TAB>product
```

`--common` supplies the COMMON block (DBLINK BioProject/BioSample, SUBMITTER, REFERENCE,
DATE.hold_date, ASSEMBLY_GAP, SOURCE …) as JSON, validated by pydantic (the `gff2mss`
console script requires JSON). `--sequence-roles roles.tsv`
(`#seq_id  type  seq_name  status  topology`) sets organelle topology (e.g. circular) and
`/organelle`. Use `--locus-tag-start` to continue numbering across companion submissions
(e.g. organelle after nuclear).
