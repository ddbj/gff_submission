# Retrieval test — gff-to-ddbj skill

Verifies that a fresh agent, given ONLY this skill (`SKILL.md` + `references/`), produces the
correct command sequence and avoids the named pitfalls. Mirrors the RED/GREEN method used for
the `ddbj-validator` skill.

## Scenario (given to a fresh agent)

> On the NIG/DDBJ cluster (non-interactive shell, no TTY), take `ann.gff3` + `genome.fa`
> (a large eukaryotic genome) to a DDBJ MSS submission and validate it. You have `mss.toml`
> and `common.json`. Give the exact command sequence (no execution).

The agent is told to read only `hermes/skills/gff-to-ddbj/SKILL.md` and the files in
`hermes/skills/gff-to-ddbj/references/`, and nothing else.

## Answer key (must appear, in order)

Executables resolve from config: `PY=${gff_to_ddbj.env_bin}/python`,
`GFF2MSS=${gff_to_ddbj.env_bin}/gff2mss`, `VALIDATOR=${gff_to_ddbj.validator_dir}/ddbj-validator`.

1. `$PY -m ddbj_gff.normalize --gff ann.gff3 --fasta genome.fa --out norm.gff --report normalize.txt`
2. `$PY -m ddbj_gff.validate --gff norm.gff` (detect-only; inspect ERROR lines)
3. `$PY -m ddbj_gff.repair --gff norm.gff --fasta genome.fa --apply all --out repaired.gff`
   — **`--fasta` present** (sequence ops require it; CLI errors exit 2 without it)
4. `$GFF2MSS --gff repaired.gff --fasta genome.fa --mss-config mss.toml --common common.json --out submission/NAME`
5. `$VALIDATOR ddbj submission -o submission/out -f -j 1` — **`-f`** and **run in the background**
   (large genome; never a foreground timeout)

## Pitfalls that MUST be avoided

- Running `repair` without `--fasta`.
- Running `ddbj-validator` without `-f` (would hit the interactive prompt / `EOFError` in non-TTY).
- Running the validator in the foreground for a large genome (risks a truncated `fixed/`).

## Result (2026-07-23)

GREEN — see the SDD task-4 report. A fresh subagent given only the skill produced all five
steps in order with `--fasta` on repair, `-f` on the validator, and backgrounded the validator
for the large genome; no pitfalls triggered.
