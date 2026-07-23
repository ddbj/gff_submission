# Hermes skill `gff-to-ddbj` — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Package this repo's `ddbj_gff` tools (plus sibling `gff2mss` and external `ddbj-validator`) as one end-to-end **Hermes agent skill** `gff-to-ddbj`, version-controlled in-repo and installable into `~/.hermes/skills/`.

**Architecture:** A Hermes skill is a directory with `SKILL.md` (YAML frontmatter + Markdown body), optional `references/` (partial-load detail), and `scripts/`. Source lives in the repo at `hermes/skills/gff-to-ddbj/`; it is copied/symlinked to `~/.hermes/skills/bioinformatics/gff-to-ddbj/` at deploy time. The skill documents the workflow normalize → validate → repair → gff2mss → ddbj-validator; a thin bash driver chains the first four steps.

**Tech Stack:** Markdown + YAML frontmatter (agentskills.io + `metadata.hermes` extensions), bash (standard tooling only), the `ddbj_gff` module CLIs and `gff2mss` console script (Python), `ddbj-validator`. Verification uses PyYAML (present in the mss_tools env) and the cluster tool installs.

## Global Constraints

- Hermes skill format per <https://www.glukhov.org/ja/ai-systems/hermes/authoring-hermes-skill/>: `SKILL.md` = YAML frontmatter (`---`…`---`) + Markdown body; required frontmatter fields `name`, `description`; Hermes extensions under `metadata.hermes`.
- `description` states ONLY when-to-use — NO workflow/process summary.
- `name: gff-to-ddbj` (lowercase, hyphens); body sections exactly, in order: `## When to use`, `## Quick reference`, `## Procedure`, `## Pitfalls`, `## Verification`.
- Source of truth in repo: `hermes/skills/gff-to-ddbj/`. Install target: `~/.hermes/skills/bioinformatics/gff-to-ddbj/`.
- Two `metadata.hermes.config` keys, with these exact defaults (this NIG/DDBJ cluster), overridable:
  - `gff_to_ddbj.env_bin` default `/lustre9/open/home/yt/micromamba/envs/mss_tools/bin` (provides `python`, `gff2mss`)
  - `gff_to_ddbj.validator_dir` default `/home/w3const/ddbj-validator-production`
- Executable resolution in all docs/scripts: `${gff_to_ddbj.env_bin}/python -m ddbj_gff.{normalize,validate,repair}`, `${gff_to_ddbj.env_bin}/gff2mss`, `${gff_to_ddbj.validator_dir}/ddbj-validator`.
- Exact tool syntax (verified):
  - `python -m ddbj_gff.normalize --gff GFF [--fasta FA] [--config C] [--taxid N] [--transl-table N] [--out OUT] [--report R]`
  - `python -m ddbj_gff.validate --gff GFF [--severity CODE=LEVEL]` (detect-only)
  - `python -m ddbj_gff.repair [--list] [--gff GFF] [--fasta FA] [--transl-table N] [--detect] [--json] [--only OPS] [--apply OPS|all] [--out OUT] [--report R]`; ops: `internal-stop-to-misc`, `utr-absent-to-partial-mrna`, `missing-start-stop-to-partial-cds`; sequence ops need `--fasta`; default apply order runs `internal-stop-to-misc` first.
  - `gff2mss --gff GFF --fasta FA --mss-config TOML --common FILE [--sequence-roles TSV] [--submission-category CAT] [--locus-tag-start N] --out OUTPREFIX` → writes `OUTPREFIX.ann` + `OUTPREFIX.fasta`.
  - `ddbj-validator ddbj <dir|targets> -o OUT [-f] [-j N] [-n|-l]` (positional target; `-f` applies fixes non-interactively → `OUT/fixed/`; large genome → background + `-j 1`).
- No new runtime dependencies; scripts use only standard tooling already present.
- Run all verification with the mss_tools env binaries (absolute paths above). Do NOT rely on `ddbj-validator`/`gff2mss` being on `$PATH`.
- Commit trailer on every commit:
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>` and the `Claude-Session:` line from repo history.

---

### Task 1: SKILL.md (frontmatter + body) and install README

**Files:**
- Create: `hermes/skills/gff-to-ddbj/SKILL.md`
- Create: `hermes/README.md`
- Test: `hermes/skills/gff-to-ddbj/tests/test_skill_md.py`

**Interfaces:**
- Produces: the skill entrypoint `SKILL.md` with frontmatter fields `name`, `description`, `version`, `platforms`, and `metadata.hermes` (`tags`, `requires_toolsets`, `config` list of `{key, description, default, prompt}`); body with the five required `##` sections. Later tasks (references, driver) rely on the config key names `gff_to_ddbj.env_bin` and `gff_to_ddbj.validator_dir`.

- [ ] **Step 1: Write the failing test**

```python
# hermes/skills/gff-to-ddbj/tests/test_skill_md.py
import re
from pathlib import Path

import yaml  # PyYAML, present in the mss_tools env

SKILL = Path(__file__).resolve().parents[1] / "SKILL.md"


def _frontmatter_and_body(text):
    m = re.match(r"^---\n(.*?)\n---\n(.*)$", text, re.DOTALL)
    assert m, "SKILL.md must start with a YAML frontmatter block delimited by ---"
    return yaml.safe_load(m.group(1)), m.group(2)


def test_frontmatter_required_and_hermes_config():
    fm, _ = _frontmatter_and_body(SKILL.read_text())
    assert fm["name"] == "gff-to-ddbj"
    assert isinstance(fm["description"], str) and fm["description"].strip()
    assert fm["platforms"] == ["linux"]
    cfg = {c["key"]: c for c in fm["metadata"]["hermes"]["config"]}
    assert set(cfg) == {"gff_to_ddbj.env_bin", "gff_to_ddbj.validator_dir"}
    assert cfg["gff_to_ddbj.env_bin"]["default"] == \
        "/lustre9/open/home/yt/micromamba/envs/mss_tools/bin"
    assert cfg["gff_to_ddbj.validator_dir"]["default"] == \
        "/home/w3const/ddbj-validator-production"
    for c in cfg.values():
        assert {"key", "description", "default", "prompt"} <= set(c)
    assert "terminal" in fm["metadata"]["hermes"]["requires_toolsets"]


def test_description_has_no_workflow_summary():
    fm, _ = _frontmatter_and_body(SKILL.read_text())
    # description is triggering-conditions only; must start with "Use when"
    assert fm["description"].lstrip().startswith("Use when")


def test_body_has_five_sections_in_order():
    _, body = _frontmatter_and_body(SKILL.read_text())
    heads = re.findall(r"^## (.+)$", body, re.MULTILINE)
    assert heads[:5] == ["When to use", "Quick reference", "Procedure",
                         "Pitfalls", "Verification"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/lustre9/open/home/yt/micromamba/envs/mss_tools/bin/python -m pytest hermes/skills/gff-to-ddbj/tests/test_skill_md.py -v`
Expected: FAIL — `SKILL.md` does not exist (FileNotFoundError).

- [ ] **Step 3: Write `SKILL.md`**

```markdown
---
name: gff-to-ddbj
description: "Use when converting an INSDC/SO GFF3 + genome FASTA into a DDBJ MSS submission (.ann/.fasta) on the NIG/DDBJ cluster — covers normalize, validate, repair (internal-stop→misc_feature and partial-CDS/mRNA), gff2mss, and the final ddbj-validator check."
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
and a common metadata file (JSON or TSV); optionally a `sequence_roles.tsv` (organelle
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
   `$GFF2MSS --gff repaired.gff --fasta GENOME.fa --mss-config MSS.toml --common COMMON.(json|tsv) --out submission/NAME`
   (add `--sequence-roles roles.tsv` for organelle topology; `--locus-tag-start N` to continue numbering).
5. **Validate the MSS** — final check (see `references/validator.md`):
   `$VALIDATOR ddbj submission -o submission/out -f -j 1`
   For a large (GB-scale) genome, run this in the background, not a foreground timeout, and always pass `-f` in a non-TTY shell.

Steps 1–4 for the happy path:
`scripts/gff_to_ddbj.sh --gff IN.gff --fasta GENOME.fa --mss-config MSS.toml --common COMMON.tsv --out-prefix submission/NAME --workdir work`

## Pitfalls

- **Repair sequence ops silently no-op without `--fasta`** — always pass `--fasta` to `repair` when you want internal-stop / start-stop detection. Apply order matters: run `internal-stop-to-misc` before the partial ops (the default `--apply all` order does this).
- **`gff2mss` requires `--mss-config` and `--common`** — it errors without them. The `.ann` entry (seqid) names must match the FASTA.
- **`ddbj-validator` without `-f` crashes with `EOFError`** at the interactive auto-fix prompt in a non-TTY shell. `-f` applies fixes non-interactively to `OUT/fixed/` (non-destructive).
- **Large genome + foreground timeout → truncated `fixed/`** — background the validator; use `-j 1` on a single-core node.
- **Validator reports describe the INPUT**, not the fixed output — re-run the validator on `OUT/fixed/` to confirm fixes cleared warnings.

## Verification

- `ddbj_gff.validate` on the normalized/repaired GFF reports **0 ERROR**.
- `submission/NAME.ann` and `submission/NAME.fasta` have **matching entry (seqid) names**.
- `ddbj-validator` summary shows **Error 0** (warnings triaged).
- If fixes were applied, `OUT/fixed/NAME.ann` is complete (its last locus_tag/gene matches the input's last — not truncated).
```

- [ ] **Step 4: Write `hermes/README.md`**

```markdown
# Hermes skills (in-repo source)

Source for Hermes agent skills that expose this repo's tools. Skills are authored per
the Hermes format (agentskills.io + `metadata.hermes`).

## gff-to-ddbj

End-to-end GFF3 → DDBJ MSS workflow (normalize → validate → repair → gff2mss →
ddbj-validator).

### Install

Copy or symlink the skill into the Hermes skills tree (category `bioinformatics`):

```bash
mkdir -p ~/.hermes/skills/bioinformatics
ln -s "$(pwd)/hermes/skills/gff-to-ddbj" ~/.hermes/skills/bioinformatics/gff-to-ddbj
```

### Config keys (Hermes `config.yaml` → `skills.config`)

| key | default (NIG/DDBJ cluster) | purpose |
|---|---|---|
| `gff_to_ddbj.env_bin` | `/lustre9/open/home/yt/micromamba/envs/mss_tools/bin` | env with `python` (ddbj_gff) + `gff2mss` |
| `gff_to_ddbj.validator_dir` | `/home/w3const/ddbj-validator-production` | `ddbj-validator` install |

Override per install if the tools live elsewhere.
```

- [ ] **Step 5: Run test to verify it passes**

Run: `/lustre9/open/home/yt/micromamba/envs/mss_tools/bin/python -m pytest hermes/skills/gff-to-ddbj/tests/test_skill_md.py -v`
Expected: PASS (3 tests).

- [ ] **Step 6: Commit**

```bash
git add hermes/skills/gff-to-ddbj/SKILL.md hermes/README.md hermes/skills/gff-to-ddbj/tests/test_skill_md.py
git commit -m "feat(hermes): gff-to-ddbj SKILL.md + install README"
```

---

### Task 2: references/ per-tool detail pages

**Files:**
- Create: `hermes/skills/gff-to-ddbj/references/normalize.md`
- Create: `hermes/skills/gff-to-ddbj/references/validate.md`
- Create: `hermes/skills/gff-to-ddbj/references/repair.md`
- Create: `hermes/skills/gff-to-ddbj/references/gff2mss.md`
- Create: `hermes/skills/gff-to-ddbj/references/validator.md`
- Test: `hermes/skills/gff-to-ddbj/tests/test_references.py`

**Interfaces:**
- Consumes: config key names from Task 1 (`gff_to_ddbj.env_bin`, `gff_to_ddbj.validator_dir`).
- Produces: five self-contained reference pages that document exact flags for each step.

- [ ] **Step 1: Write the failing test**

```python
# hermes/skills/gff-to-ddbj/tests/test_references.py
from pathlib import Path

REF = Path(__file__).resolve().parents[1] / "references"

# each reference file must exist and mention these key tokens (exact flags/ops)
REQUIRED = {
    "normalize.md": ["ddbj_gff.normalize", "--gff", "--out", "--report", "--transl-table"],
    "validate.md":  ["ddbj_gff.validate", "detect-only", "ERR:"],
    "repair.md":    ["ddbj_gff.repair", "--detect", "--apply", "internal-stop-to-misc",
                     "utr-absent-to-partial-mrna", "missing-start-stop-to-partial-cds"],
    "gff2mss.md":   ["gff2mss", "--mss-config", "--common", "--sequence-roles", ".ann"],
    "validator.md": ["ddbj-validator", "-f", "-j 1", "fixed/", "background"],
}


def test_reference_pages_exist_and_document_flags():
    for fname, tokens in REQUIRED.items():
        text = (REF / fname).read_text()
        for tok in tokens:
            assert tok in text, f"{fname} missing {tok!r}"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/lustre9/open/home/yt/micromamba/envs/mss_tools/bin/python -m pytest hermes/skills/gff-to-ddbj/tests/test_references.py -v`
Expected: FAIL — references directory / files do not exist.

- [ ] **Step 3: Write the five reference pages**

`references/normalize.md`:
```markdown
# ddbj_gff.normalize (reference)

`$PY -m ddbj_gff.normalize --gff GFF [--fasta FA] [--config C.toml] [--taxid N] [--transl-table N] [--insdc-gff-version V] [--out OUT] [--report R]`

Canonicalizes a GFF3 toward the INSDC profile (adds directives; coerces coding
`transcript`→`mRNA`; wraps/reparents CDS under mRNA; SO term mapping; trans-splicing
location; recoded/anticodon children). Pass `--fasta` so `##sequence-region` lengths are
exact (otherwise approximated from max feature end). `--report` writes an applied /
needs-attention summary — review `needs-manual` and `unmapped-type` lines. Without `--out`
the normalized GFF goes to stdout; without `--report` the report goes to stderr.
```

`references/validate.md`:
```markdown
# ddbj_gff.validate (reference)

`$PY -m ddbj_gff.validate --gff GFF [--severity CODE=LEVEL]`

Detect-only INSDC profile validator (does not modify the file). Emits diagnostics with a
severity (`ERR:` / `WAR:` / `INFO:`) and a code. Fix inputs until there are no `ERR:`
lines; warnings may legitimately remain. `--severity CODE=off|error|warning|info` overrides
a rule's level. Unknown/extra GFF attributes are left untouched (not flagged). Run after
`normalize` and again after `repair`.
```

`references/repair.md`:
```markdown
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
```

`references/gff2mss.md`:
```markdown
# gff2mss (reference)

`$GFF2MSS --gff GFF --fasta FA --mss-config MSS.toml --common COMMON.(json|tsv) [--sequence-roles roles.tsv] [--submission-category CAT] [--locus-tag-start N] --out OUTPREFIX`

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
DATE.hold_date, ASSEMBLY_GAP, SOURCE …) as JSON or TSV. `--sequence-roles roles.tsv`
(`#seq_id  type  seq_name  status  topology`) sets organelle topology (e.g. circular) and
`/organelle`. Use `--locus-tag-start` to continue numbering across companion submissions
(e.g. organelle after nuclear).
```

`references/validator.md`:
```markdown
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `/lustre9/open/home/yt/micromamba/envs/mss_tools/bin/python -m pytest hermes/skills/gff-to-ddbj/tests/test_references.py -v`
Expected: PASS.

- [ ] **Step 5: Verify the documented commands resolve on this cluster**

Run each and confirm exit 0 / usage text (proves the resolution paths in the docs are real):
```bash
EB=/lustre9/open/home/yt/micromamba/envs/mss_tools/bin
$EB/python -m ddbj_gff.normalize --help >/dev/null && echo normalize-ok
$EB/python -m ddbj_gff.validate --help >/dev/null && echo validate-ok
$EB/python -m ddbj_gff.repair --help >/dev/null && echo repair-ok
$EB/gff2mss --help >/dev/null && echo gff2mss-ok
/home/w3const/ddbj-validator-production/ddbj-validator ddbj --help >/dev/null 2>&1; echo "validator-exit=$?"
```
Expected: `normalize-ok`, `validate-ok`, `repair-ok`, `gff2mss-ok` printed; validator prints a usage line (exit 0 or 2, not 127).

- [ ] **Step 6: Commit**

```bash
git add hermes/skills/gff-to-ddbj/references hermes/skills/gff-to-ddbj/tests/test_references.py
git commit -m "feat(hermes): gff-to-ddbj per-tool reference pages"
```

---

### Task 3: scripts/gff_to_ddbj.sh happy-path driver + dry-run

**Files:**
- Create: `hermes/skills/gff-to-ddbj/scripts/gff_to_ddbj.sh`
- Test: `hermes/skills/gff-to-ddbj/tests/test_driver.py`

**Interfaces:**
- Consumes: `${gff_to_ddbj.env_bin}` (via env var `GFF_TO_DDBJ_ENV_BIN` or `--env-bin`), the tool CLIs from the Global Constraints.
- Produces: an executable `gff_to_ddbj.sh` chaining normalize → validate → repair(apply all) → gff2mss, writing `<out-prefix>.ann` and `<out-prefix>.fasta`.

- [ ] **Step 1: Write the failing test**

```python
# hermes/skills/gff-to-ddbj/tests/test_driver.py
import os
import shutil
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "gff_to_ddbj.sh"
ENV_BIN = "/lustre9/open/home/yt/micromamba/envs/mss_tools/bin"
FIX = Path("/lustre9/open/home/yt/ddbj/ddbj_mss_tools/tests/mss_fixtures")


def test_driver_produces_ann_and_fasta(tmp_path):
    # stage a known-good minimal gff2mss fixture set
    for f in ("mini.gff3", "mini.fa", "config.toml", "common.metadata.tsv"):
        shutil.copy(FIX / f, tmp_path / f)
    out_prefix = tmp_path / "submission" / "mini"
    r = subprocess.run(
        ["bash", str(SCRIPT),
         "--gff", str(tmp_path / "mini.gff3"),
         "--fasta", str(tmp_path / "mini.fa"),
         "--mss-config", str(tmp_path / "config.toml"),
         "--common", str(tmp_path / "common.metadata.tsv"),
         "--out-prefix", str(out_prefix),
         "--workdir", str(tmp_path / "work")],
        env={**os.environ, "GFF_TO_DDBJ_ENV_BIN": ENV_BIN},
        capture_output=True, text=True)
    assert r.returncode == 0, f"driver failed:\nSTDOUT{r.stdout}\nSTDERR{r.stderr}"
    ann = out_prefix.with_suffix(".ann")
    fasta = out_prefix.with_suffix(".fasta")
    assert ann.exists() and ann.stat().st_size > 0
    assert fasta.exists() and fasta.stat().st_size > 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/lustre9/open/home/yt/micromamba/envs/mss_tools/bin/python -m pytest hermes/skills/gff-to-ddbj/tests/test_driver.py -v`
Expected: FAIL — script does not exist.

- [ ] **Step 3: Write `scripts/gff_to_ddbj.sh`**

```bash
#!/usr/bin/env bash
# gff_to_ddbj.sh — happy-path driver: normalize -> validate -> repair(apply all) -> gff2mss.
# Env-specific paths come from GFF_TO_DDBJ_ENV_BIN (or --env-bin); default = NIG/DDBJ cluster.
set -euo pipefail

ENV_BIN="${GFF_TO_DDBJ_ENV_BIN:-/lustre9/open/home/yt/micromamba/envs/mss_tools/bin}"
GFF="" FASTA="" MSSCFG="" COMMON="" OUTP="" WORK="."

usage() {
  echo "usage: $0 --gff IN.gff --fasta IN.fa --mss-config MSS.toml \\"
  echo "          --common COMMON.(json|tsv) --out-prefix DIR/NAME [--workdir DIR] [--env-bin DIR]"
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
```

Make it executable:
```bash
chmod +x hermes/skills/gff-to-ddbj/scripts/gff_to_ddbj.sh
```

- [ ] **Step 4: Run test to verify it passes**

Run: `/lustre9/open/home/yt/micromamba/envs/mss_tools/bin/python -m pytest hermes/skills/gff-to-ddbj/tests/test_driver.py -v`
Expected: PASS. If `gff2mss` rejects the `mini` set through the normalize→repair chain, first reproduce by running the four commands by hand on the staged fixture to find which step fails, fix the driver's argument passing (not the assertions), and re-run. Do not weaken the assertions.

- [ ] **Step 5: Commit**

```bash
git add hermes/skills/gff-to-ddbj/scripts/gff_to_ddbj.sh hermes/skills/gff-to-ddbj/tests/test_driver.py
git commit -m "feat(hermes): gff-to-ddbj happy-path driver + dry-run test"
```

---

### Task 4: Retrieval verification (skill efficacy) + gap fixes

**Files:**
- Create: `hermes/skills/gff-to-ddbj/tests/RETRIEVAL_TEST.md` (records the scenario + expected answer key)
- Modify (only if gaps found): `hermes/skills/gff-to-ddbj/SKILL.md` and/or `references/*.md`

**Interfaces:**
- Consumes: the finished `SKILL.md` + `references/`.
- Produces: evidence that a fresh agent, given only the skill, produces the correct command sequence and avoids the named pitfalls.

- [ ] **Step 1: Write the retrieval scenario + answer key**

Create `hermes/skills/gff-to-ddbj/tests/RETRIEVAL_TEST.md` with:
- **Scenario:** "On the NIG/DDBJ cluster (non-interactive shell), take `ann.gff3` + `genome.fa` to a DDBJ MSS submission and validate it. You have `mss.toml` and `common.json`."
- **Answer key (must appear):** (1) `${env_bin}/python -m ddbj_gff.normalize … --out norm.gff`; (2) `… ddbj_gff.validate --gff norm.gff`; (3) `… ddbj_gff.repair --gff norm.gff --fasta genome.fa --apply all --out repaired.gff` (with `--fasta`); (4) `${env_bin}/gff2mss --gff repaired.gff --fasta genome.fa --mss-config mss.toml --common common.json --out submission/NAME`; (5) `${validator_dir}/ddbj-validator ddbj submission -o … -f -j 1`, **backgrounded**.
- **Pitfalls that must be avoided:** running `repair` without `--fasta`; running `ddbj-validator` without `-f`; running the validator in the foreground for a large genome.

- [ ] **Step 2: Run the retrieval test (fresh subagent)**

Dispatch a fresh general-purpose subagent whose ONLY provided knowledge is the contents of
`hermes/skills/gff-to-ddbj/SKILL.md` (tell it to also read the `references/` files in that
directory). Give it the Scenario from Step 1. Ask for the exact command sequence (no
execution). Compare its answer against the answer key.

Expected: the subagent lists the five steps in order with `--fasta` on repair, `-f` on the
validator, and backgrounds the validator for the large genome.

- [ ] **Step 3: Close gaps (only if the subagent missed something)**

If the subagent omitted `--fasta`, `-f`, backgrounding, or got the order wrong, edit
`SKILL.md` / the relevant `references/*.md` to make that guidance more prominent (e.g. move
it into the Pitfalls list or the Procedure step), then re-run Step 2 until the answer key is met.

- [ ] **Step 4: Run the full skill test suite**

Run: `/lustre9/open/home/yt/micromamba/envs/mss_tools/bin/python -m pytest hermes/skills/gff-to-ddbj/tests -v`
Expected: all tests pass (SKILL.md structure, references, driver).

- [ ] **Step 5: Commit**

```bash
git add hermes/skills/gff-to-ddbj/tests/RETRIEVAL_TEST.md hermes/skills/gff-to-ddbj
git commit -m "test(hermes): gff-to-ddbj retrieval verification + any gap fixes"
```

---

## Self-Review

**1. Spec coverage:**
- One end-to-end skill `gff-to-ddbj` → Tasks 1–3. ✓
- Hermes format (frontmatter + `metadata.hermes` + body sections) → Task 1 (+ test). ✓
- This-cluster defaults, config-overridable → Task 1 config keys + README (Task 1). ✓
- SKILL.md + references/ + scripts/ structure → Tasks 1, 2, 3. ✓
- references: normalize/validate/repair/gff2mss/validator → Task 2. ✓
- scripts happy-path driver (steps 1–4; detect-preview + validator stay manual) → Task 3 + SKILL Procedure. ✓
- In-repo source + install to ~/.hermes/skills/bioinformatics/ → Task 1 README. ✓
- Testing: frontmatter validity, executable resolution, driver dry-run, retrieval test → Tasks 1, 2, 3, 4. ✓
- Non-goals (MCP, config auto-gen, Hermes install, heterosigma preprocessing) → not built. ✓

**2. Placeholder scan:** No TBD/TODO; every file's full content is provided; every command has an expected result. ✓

**3. Type/name consistency:** Config keys `gff_to_ddbj.env_bin` / `gff_to_ddbj.validator_dir` identical across SKILL.md, README, references, driver (`GFF_TO_DDBJ_ENV_BIN`), and tests. Body section names match the test's expected list. Operation names match `ddbj_gff.repair`. ✓

## Notes for the implementer

- Run every command with the mss_tools env binaries (absolute paths); the tools are NOT on `$PATH`.
- Task 3's fixture (`ddbj_mss_tools/tests/mss_fixtures/{mini.gff3,mini.fa,config.toml,common.metadata.tsv}`) is a known-good gff2mss input; if the normalize→repair→gff2mss chain fails on it, debug the driver by running the four steps by hand on the staged copy — do not relax the test.
- The skill's source lives in-repo; installing to `~/.hermes/skills/` is a manual deploy step documented in `hermes/README.md`, not part of these tasks.
