---
name: gff-to-ddbj
description: "Use when you have an INSDC/SO GFF3 gene annotation and its genome FASTA and need to produce or validate a DDBJ MSS submission (.ann + .fasta) on the NIG/DDBJ cluster."
version: 0.1.0
platforms: [linux]
metadata:
  hermes:
    tags: [bioinformatics, ddbj, gff3, mss, submission, insdc]
    requires_toolsets: [terminal]
    config:
      - key: gff_to_ddbj.env_bin
        description: bin dir of the Python env holding ddbj_gff + gff2mss (provides python and gff2mss)
        default: /lustre9/open/home/yt/micromamba/envs/mss_tools/bin
        prompt: "Path to the mss_tools env bin directory"
      - key: gff_to_ddbj.validator_dir
        description: ddbj-validator production install directory
        default: /home/w3const/ddbj-validator-production
        prompt: "Path to ddbj-validator-production"
---

# GFF3 → DDBJ MSS submission (ddbj_gff)

## When to use

You have an INSDC/SO GFF3 gene annotation plus its genome FASTA and need a DDBJ MSS
submission (`.ann` + `.fasta`): to canonicalize a GFF3 toward the INSDC profile, curate
untranslatable/partial features, convert to MSS, and validate the result before submission.
Not for building the `gff2mss` metadata config itself, nor for MetaboBank/BioSample/DRA records.

## Quick reference

Resolve executables from Hermes config (defaults are for the NIG/DDBJ cluster):

- `PY = ${skills.config.gff_to_ddbj.env_bin}/python` — run module CLIs as `$PY -m ddbj_gff.<tool>`
- `GFF2MSS = ${skills.config.gff_to_ddbj.env_bin}/gff2mss`
- `VALIDATOR = ${skills.config.gff_to_ddbj.validator_dir}/ddbj-validator`

Inputs you must have: the GFF3, the genome FASTA, and — for `gff2mss` — an mss-config TOML
and a common metadata JSON file; optionally a `sequence_roles.tsv` (organelle
topology) and a `product_map.tsv`. See `references/gff2mss.md`.

Outputs: `norm.gff` (normalized), `repaired.gff` (curated), `submission/NAME.ann` +
`submission/NAME.fasta` (MSS), and the validator's `reports/`, `fixed/`, `aa/`.

Per-tool detail lives in `references/` — load the one you need:
`normalize.md`, `validate.md`, `repair.md`, `gff2mss.md`, `validator.md`.

A happy-path driver chains steps 1–4: `scripts/gff_to_ddbj.sh` (see `references/` for the
detect-only preview and the validator step, which stay manual).

## Procedure

1. **Normalize** — canonicalize the GFF3 toward the INSDC profile:
   `$PY -m ddbj_gff.normalize --gff IN.gff --fasta GENOME.fa --out norm.gff --report normalize.txt`
   Review `normalize.txt` for `needs-manual` / `unmapped-type` lines.
2. **Validate** — detect-only INSDC profile check:
   `$PY -m ddbj_gff.validate --gff norm.gff`
   Inspect `ERR:` lines; fix inputs and re-run until no ERROR (warnings may remain).
3. **Repair** — preview, then apply curation (sequence-based ops need `--fasta`):
   - Preview: `$PY -m ddbj_gff.repair --gff norm.gff --fasta GENOME.fa --detect --json`
   - Apply chosen ops (or all): `$PY -m ddbj_gff.repair --gff norm.gff --fasta GENOME.fa --apply all --out repaired.gff --report repair.txt`
   Ops: `internal-stop-to-misc` (CDS with internal stop → `misc_feature`), `utr-absent-to-partial-mrna`, `missing-start-stop-to-partial-cds`. `internal-stop-to-misc` applies first.
4. **gff2mss** — convert to MSS:
   `$GFF2MSS --gff repaired.gff --fasta GENOME.fa --mss-config MSS.toml --common COMMON.json --out submission/NAME`
   (add `--sequence-roles roles.tsv` for organelle topology; `--locus-tag-start N` to continue numbering).
5. **Validate the MSS** — final check (see `references/validator.md`):
   `$VALIDATOR ddbj submission -o submission/out -f -j 1`
   For a large (GB-scale) genome, run this in the background, not a foreground timeout, and always pass `-f` in a non-TTY shell.

Steps 1–4 for the happy path:
`scripts/gff_to_ddbj.sh --gff IN.gff --fasta GENOME.fa --mss-config MSS.toml --common COMMON.json --out-prefix submission/NAME --workdir work`

## Pitfalls

- **Repair sequence ops require `--fasta`** — the sequence-based ops (`internal-stop-to-misc`, `missing-start-stop-to-partial-cds`) error (exit 2) without `--fasta`. Apply order matters: run `internal-stop-to-misc` before the partial ops (the default `--apply all` order does this).
- **`gff2mss` requires `--mss-config` and `--common`** — it errors without them. The `.ann` entry (seqid) names must match the FASTA.
- **`ddbj-validator` without `-f` crashes with `EOFError`** at the interactive auto-fix prompt in a non-TTY shell. `-f` applies fixes non-interactively to `OUT/fixed/` (non-destructive).
- **Large genome + foreground timeout → truncated `fixed/`** — background the validator; use `-j 1` on a single-core node.
- **Validator reports describe the INPUT**, not the fixed output — re-run the validator on `OUT/fixed/` to confirm fixes cleared warnings.

## Verification

- `ddbj_gff.validate` on the normalized/repaired GFF reports **0 ERROR**.
- `submission/NAME.ann` and `submission/NAME.fasta` have **matching entry (seqid) names**.
- `ddbj-validator` summary shows **Error 0** (warnings triaged).
- If fixes were applied, `OUT/fixed/NAME.ann` is complete (its last locus_tag/gene matches the input's last — not truncated).
