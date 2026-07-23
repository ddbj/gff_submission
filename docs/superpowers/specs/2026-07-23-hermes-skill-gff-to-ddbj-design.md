# Hermes skill `gff-to-ddbj` — design

**Date:** 2026-07-23
**Status:** design (approved, for spec review)
**Scope:** Package this repository's `ddbj_gff` tools (plus the sibling `gff2mss` and the external
`ddbj-validator`) as a single end-to-end **Hermes agent skill** so a Hermes agent can drive the
GFF3 → DDBJ MSS submission workflow. Follows the Hermes skill format
(agentskills.io base + `metadata.hermes` extensions), per
<https://www.glukhov.org/ja/ai-systems/hermes/authoring-hermes-skill/>.

## Problem / motivation

`gff_submission` exposes its functionality as Python module CLIs (`python -m ddbj_gff.normalize`
/ `.validate` / `.repair` / `.flatfile.cli`) — a library, not console scripts. A Hermes agent (or
any agent-AI) that wants to take a GFF3 + genome FASTA to a DDBJ MSS submission currently has to
rediscover: which env has the tools, the correct step order, the non-obvious traps (repair apply
order, `gff2mss` required config, `ddbj-validator`'s `-f`/EOF and background-for-large-genomes
behavior). A baseline test (a capable fresh agent, no skill) confirmed these are not discoverable:
it could not locate the canonical installs, guessed at flags, and had no correct step sequence.

A Hermes skill captures the workflow once as an on-demand knowledge document with progressive
disclosure, so the agent loads a compact index, then the full procedure only when the task matches.

## Decisions (confirmed)

1. **One end-to-end skill** (`gff-to-ddbj`), not per-tool skills. It orchestrates the full pipeline:
   normalize → validate → repair → gff2mss → ddbj-validator.
2. **Hermes skill format** (not MCP, not a JSON tool schema): `SKILL.md` with YAML frontmatter +
   Markdown body, `metadata.hermes` block, `references/` for per-tool detail, optional `scripts/`.
3. **This-cluster defaults, overridable via `metadata.hermes.config`.** Env-specific paths default
   to the known NIG/DDBJ cluster locations and can be overridden per install via Hermes config.
4. **Structure = SKILL.md (workflow) + references/ (per-tool detail) + scripts/ (happy-path driver)**
   — the progressive-disclosure layout Hermes prescribes.
5. **Scope note:** "this repo's tools" = `ddbj_gff`; the end-to-end choice means `gff2mss` (sibling
   `ddbj_mss_tools`) and `ddbj-validator` (external) are referenced as pipeline **steps**. The skill
   centers on this repo's tools and links out to the others' invocation detail.

## Skill location & discovery

```
~/.hermes/skills/bioinformatics/gff-to-ddbj/
├── SKILL.md
├── references/
│   ├── normalize.md
│   ├── validate.md
│   ├── repair.md
│   ├── gff2mss.md
│   └── validator.md
└── scripts/
    └── gff_to_ddbj.sh
```

- Category folder `bioinformatics/`; skill folder name = `name` field = `gff-to-ddbj`; slash route
  `/gff-to-ddbj`.
- **Source of truth lives in this repo** under `hermes/skills/gff-to-ddbj/` (version-controlled with
  the tools it documents); it is installed/symlinked into `~/.hermes/skills/bioinformatics/` at
  deploy time. This keeps the skill in lockstep with `ddbj_gff` changes.

## Frontmatter (SKILL.md)

```yaml
name: gff-to-ddbj
description: "Use when converting an INSDC/SO GFF3 + genome FASTA into a DDBJ MSS submission
  (.ann/.fasta) on the NIG/DDBJ cluster — normalize, validate, repair (internal-stop→misc_feature /
  partial), gff2mss, and the final ddbj-validator check."
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
```

- `description` states ONLY when-to-use (no workflow summary), per skill-authoring guidance.
- Executables resolve as: `${gff_to_ddbj.env_bin}/python -m ddbj_gff.*`, `${gff_to_ddbj.env_bin}/gff2mss`,
  `${gff_to_ddbj.validator_dir}/ddbj-validator`.

## SKILL.md body sections

1. **When to use** — plain-language triggers (have a GFF3 + genome FASTA, want a DDBJ MSS submission).
2. **Quick reference** — executable resolution from config; the additional inputs `gff2mss` needs
   (mss-config TOML, common JSON, and optional product_map TSV / sequence_roles TSV); output layout
   (`norm.gff`, `repaired.gff`, `submission/NAME.{ann,fasta}`, validator `reports/`/`fixed/`/`aa/`).
3. **Procedure** (ordered, ready to run):
   1. `python -m ddbj_gff.normalize --gff IN.gff [--fasta FA] --out norm.gff --report norm.txt`
   2. `python -m ddbj_gff.validate --gff norm.gff` — detect-only; inspect ERROR lines
   3. `python -m ddbj_gff.repair --gff norm.gff [--fasta FA] --detect --json` → review, then
      `python -m ddbj_gff.repair --gff norm.gff --fasta FA --apply all --out repaired.gff`
      (default op order internal-stop→misc first; requires FASTA for sequence-based ops)
   4. `gff2mss --gff repaired.gff --fasta FA --mss-config X.toml --common X.json --out submission/NAME`
   5. `ddbj-validator ddbj <submission_dir> -o <out> -f -j 1` — **run in the background for large
      genomes; always `-f` in a non-TTY shell**
4. **Pitfalls** — repair apply order; sequence ops need `--fasta`; `gff2mss` required config;
   `ddbj-validator` `-f`/EOFError-in-non-TTY; foreground-timeout truncates `fixed/`; reports describe
   the INPUT so re-validate `fixed/` to confirm; `-j 1` on a single-core node.
5. **Verification** — `validate` reports 0 ERROR; `.ann`/`.fasta` entry (seqid) names match;
   `ddbj-validator` summary shows Error 0; `fixed/.ann` complete (last locus_tag/gene matches input).

## references/ (partial-load detail)

Each is a self-contained page (Hermes loads on demand via `skill_view`):

- **normalize.md** — full flag list, what each pass does, `--config`/`--taxid`/`--transl-table`.
- **validate.md** — detect-only nature, severity codes, that unknown attributes pass.
- **repair.md** — detect/apply two-phase model, the three operations, `--list`/`--only`/`--apply`,
  INSDC partial encoding, misc_feature retype semantics.
- **gff2mss.md** — required args, minimal `mss-config`/`common` examples, sequence_roles/product_map.
- **validator.md** — the `ddbj-validator` essentials (subcommands, `-f`, background for large
  genomes, re-validate `fixed/`), self-contained so Hermes needs nothing external.

## scripts/gff_to_ddbj.sh (optional happy-path driver)

A thin bash driver chaining steps 1→4 (normalize → validate → repair --apply all → gff2mss) with
paths injected from config/env. It uses `${HERMES_SKILL_DIR}` for its own location and accepts
`--env-bin` / input/output arguments. The **detect-only preview** (step 3 review) and the
**large-genome validator run** (step 5) stay as explicit body-Procedure steps rather than being
buried in the driver, because both involve judgment (which repairs to apply) or long/background
execution. Standard-library / existing-tooling only (no new deps), per Hermes guidance.

## Testing (Hermes not installed in this environment)

- **Frontmatter validity:** parses as YAML; required `name`/`description` present; `metadata.hermes`
  structure matches the Hermes spec (config keys have key/description/default/prompt).
- **Executable resolution:** each command's `--help` succeeds on this cluster using the default
  config paths (`${env_bin}/python -m ddbj_gff.{normalize,validate,repair}`, `${env_bin}/gff2mss`,
  `${validator_dir}/ddbj-validator`).
- **Driver dry-run:** `scripts/gff_to_ddbj.sh` runs on a tiny fixture GFF+FASTA through step 4 and
  produces a `.ann`/`.fasta`.
- **Retrieval test (skill efficacy):** a fresh subagent given only `SKILL.md` + `references/` must
  produce the correct command sequence and avoid the named pitfalls for a GFF→MSS task — mirroring
  the RED/GREEN method used for the `ddbj-validator` skill.

## Non-goals / out of scope

- MCP server or JSON function-schema variants (a different consumer; not this task).
- Auto-generating `gff2mss` config (mss-config/common/product_map) — the skill documents the
  required inputs; producing them is the submitter's step.
- Installing/configuring Hermes itself, or a Hermes runtime on this cluster.
- The organelle/AGAT-specific preprocessing in `dev/heterosigma` (project-specific, not generic).

## Deployment

Source in-repo at `hermes/skills/gff-to-ddbj/`; install by copying/symlinking into
`~/.hermes/skills/bioinformatics/gff-to-ddbj/`. A short `hermes/README.md` documents the install
step and the config keys.
