# ddbj-validator (reference, MSS check)

`$VALIDATOR ddbj <target-dir-or-files> -o OUT [-f] [-j N] [-n|-l]`

Validates DDBJ `.ann` + FASTA (the `gff2mss` output). Pass the submission directory (or the
files) positionally — the `-d/-a/-s` flags are deprecated.

- `-f` / `--force-fix` — apply all auto-fixes without the interactive prompt. **Required in a
  non-TTY / agent shell**: without it the tool opens an `Action: [a]/[i]/[q]?` prompt and
  crashes with `EOFError`. Fixes are written to `OUT/fixed/` (originals untouched).
- **Large (GB-scale) genome → run in the background**, never a foreground timeout (a killed
  run leaves a truncated `OUT/fixed/`). Use `-j 1` on a single-core node (higher `-j` only
  multiplies ~15–18 GB/process memory).
- `-n` NCBI API / `-l` fully local. Outputs: `reports/{validation_report_summary,details}.txt`,
  `autofix_confirmation_summary.txt`, `fixed/`, `aa/` (CDS→protein).
- Reports describe the INPUT; re-run on `OUT/fixed/` to confirm fixes cleared warnings.

Canonical install: `${skills.config.gff_to_ddbj.validator_dir}` (default
`/home/w3const/ddbj-validator-production`).
