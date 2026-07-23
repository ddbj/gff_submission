# pass_reparent_gene_children_to_mrna Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a default-on `ddbj-gff` normalize pass that reparents a gene's direct structural sub-features (CDS/exon/intron/UTR/codon) onto the gene's single mRNA, fixing AUGUSTUS-dialect GFFs where those sub-features are parented to the gene instead of the mRNA.

**Architecture:** One new pure function `pass_reparent_gene_children_to_mrna(doc, ctx) -> list[Change]` in `src/ddbj_gff/normalize/passes.py`, gated by a new `NormalizeConfig.reparent_gene_children` bool (default `True`). Registered in `ALL_PASSES` between `pass_coerce_transcript_to_mrna` and `pass_wrap_cds_in_mrna`, with its applied action `"reparent-to-mrna"` added to the `_APPLIED` set. It mutates only feature links (`parent_ids` / `attributes["Parent"]` / `parents` / `children`) and the mRNA span; it creates and removes no features, so `doc.features` / `doc.feature_index` / `doc.roots` are untouched.

**Tech Stack:** Python 3.11+ (tomllib, dataclasses), pytest, BioPython (already a dep). Tests drive GFF text through `ddbj_gff.parse` → `ddbj_gff.normalize.normalize`.

## Global Constraints

- Single repository: `gff_submission` (package `ddbj_gff`). No `ddbj_mss_tools` change.
- Default **on**: `reparent_gene_children: bool = True`; on well-formed 3-level input the pass is a no-op (structure unchanged).
- Reparent only when the gene has **exactly one transcript child and it is an `mRNA`**, and the gene has ≥1 direct structural child. Otherwise: `≥2 transcripts` or `sole transcript non-mRNA` → record a `needs-manual` attention and leave untouched; `zero transcripts` → skip silently (leave for `pass_wrap_cds_in_mrna`).
- Structural sub-feature types eligible for reparenting: exactly `{"CDS", "exon", "intron", "three_prime_UTR", "five_prime_UTR", "start_codon", "stop_codon"}`.
- mRNA span is **extend-only** to a single `Span(seqid, lo, hi, strand)` covering the mRNA's existing extent plus all reparented children; never shrinks.
- Applied action string is exactly `"reparent-to-mrna"` (must be added to `_APPLIED`). Ambiguity action string is exactly `"needs-manual"` (already treated as attention; NOT in `_APPLIED`).
- Never use `in` / `not in` / `==` on `Feature` objects or lists of them — `Feature` is a dataclass with cyclic `parents`/`children`, so value-equality can recurse. Filter by `id()` identity.
- Test container command (run from repo root inside the dev container):
  `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest <path> -v'`

---

### Task 1: NormalizeConfig.reparent_gene_children field

**Files:**
- Modify: `src/ddbj_gff/normalize/config.py`
- Test: `tests/test_normalize_reparent.py` (new file; Task 2 appends to it)

**Interfaces:**
- Consumes: nothing (first task).
- Produces: `NormalizeConfig.reparent_gene_children: bool` (default `True`); `load_normalize_config` reads `[normalize].reparent_gene_children` (default `True`).

- [ ] **Step 1: Write the failing tests**

Create `tests/test_normalize_reparent.py` with exactly:

```python
from ddbj_gff.normalize.config import NormalizeConfig, load_normalize_config


def test_reparent_config_default_on():
    assert NormalizeConfig().reparent_gene_children is True


def test_reparent_config_loader_reads_flag(tmp_path):
    p = tmp_path / "n.toml"
    p.write_text("[normalize]\nreparent_gene_children = false\n")
    assert load_normalize_config(str(p)).reparent_gene_children is False


def test_reparent_config_loader_defaults_on(tmp_path):
    p = tmp_path / "n.toml"
    p.write_text("[normalize]\ntaxid = 1\n")
    assert load_normalize_config(str(p)).reparent_gene_children is True
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_normalize_reparent.py -v'`
Expected: FAIL — `AttributeError: 'NormalizeConfig' object has no attribute 'reparent_gene_children'` (and loader test fails too).

- [ ] **Step 3: Add the field and loader line**

In `src/ddbj_gff/normalize/config.py`, add the field to the dataclass (after `wrap_cds_in_mrna`):

```python
@dataclass
class NormalizeConfig:
    taxid: int | None = None
    transl_table: int = 1
    insdc_gff_version: str = "1.0.0"
    coerce_transcript_to_mrna: bool = True
    wrap_cds_in_mrna: bool = True
    reparent_gene_children: bool = True
    merge_overlapping_loci: bool = False
    merge_overlap_min_fraction: float = 0.0
```

And add the loader line inside `load_normalize_config`'s `NormalizeConfig(...)` call (after the `wrap_cds_in_mrna=` line):

```python
        wrap_cds_in_mrna=n.get("wrap_cds_in_mrna", True),
        reparent_gene_children=n.get("reparent_gene_children", True),
        merge_overlapping_loci=n.get("merge_overlapping_loci", False),
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_normalize_reparent.py -v'`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/normalize/config.py tests/test_normalize_reparent.py
git commit -m "feat(normalize): NormalizeConfig.reparent_gene_children (default on)"
```

---

### Task 2: pass_reparent_gene_children_to_mrna + register + behaviour tests

**Files:**
- Modify: `src/ddbj_gff/normalize/passes.py` (add the pass function)
- Modify: `src/ddbj_gff/normalize/normalize.py` (import, insert into `ALL_PASSES`, extend `_APPLIED`)
- Test: `tests/test_normalize_reparent.py` (append behaviour tests)

**Interfaces:**
- Consumes: `NormalizeConfig.reparent_gene_children` (Task 1). `Feature` fields `id, source, type, spans, attributes, parent_ids, children, parents` and `Span(seqid, start, end, strand)` from `ddbj_gff.model` (already imported at the top of `passes.py`). `Change(action, target, message)` from `.report` (already imported). `NormalizeContext` has `.config`.
- Produces: `pass_reparent_gene_children_to_mrna(doc, ctx) -> list[Change]`; applied action `"reparent-to-mrna"`, attention action `"needs-manual"`.

- [ ] **Step 1: Write the failing behaviour tests**

Append to `tests/test_normalize_reparent.py`:

```python
from ddbj_gff import parse
from ddbj_gff.normalize.normalize import normalize

HDR = "##gff-version 3\n##sequence-region c 1 100000\n"


def _norm(gff, **cfg):
    doc = parse(gff)
    work, report = normalize(doc, config=NormalizeConfig(taxid=1, **cfg))
    return work, report


def _feat(doc, fid):
    return doc.feature_index.get(fid)


def _validate(doc):
    from ddbj_gff.validate import validate
    return validate(doc)


def test_misparented_cds_exon_reparented_to_mrna():
    gff = HDR + (
        "c\tx\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "c\tx\tmRNA\t100\t500\t.\t+\t.\tID=g1.1;Parent=g1\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g1.cds1;Parent=g1\n"
        "c\tx\texon\t100\t200\t.\t+\t.\tID=g1.ex1;Parent=g1\n"
    )
    doc, report = _norm(gff)
    gene, mrna = _feat(doc, "g1"), _feat(doc, "g1.1")
    cds, exon = _feat(doc, "g1.cds1"), _feat(doc, "g1.ex1")
    assert cds.parent_ids == ["g1.1"] and exon.parent_ids == ["g1.1"]
    assert cds.attributes["Parent"] == ["g1.1"]
    assert {c.id for c in mrna.children} == {"g1.cds1", "g1.ex1"}
    assert {c.id for c in gene.children} == {"g1.1"}
    assert any(ch.action == "reparent-to-mrna" for ch in report.applied)
    assert "dangling-parent" not in {d.code for d in _validate(doc)}


def test_wellformed_three_level_is_noop():
    gff = HDR + (
        "c\tx\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "c\tx\tmRNA\t100\t500\t.\t+\t.\tID=g1.1;Parent=g1\n"
        "c\tx\tCDS\t100\t500\t.\t+\t0\tID=g1.cds1;Parent=g1.1\n"
        "c\tx\texon\t100\t500\t.\t+\t.\tID=g1.ex1;Parent=g1.1\n"
    )
    doc, report = _norm(gff)
    assert not any(ch.action == "reparent-to-mrna" for ch in report.applied)
    assert _feat(doc, "g1.cds1").parent_ids == ["g1.1"]
    assert {c.id for c in _feat(doc, "g1").children} == {"g1.1"}


def test_two_mrnas_gene_level_cds_not_reparented():
    gff = HDR + (
        "c\tx\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "c\tx\tmRNA\t100\t500\t.\t+\t.\tID=g1.1;Parent=g1\n"
        "c\tx\tmRNA\t100\t500\t.\t+\t.\tID=g1.2;Parent=g1\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g1.cds1;Parent=g1\n"
    )
    doc, report = _norm(gff)
    assert _feat(doc, "g1.cds1").parent_ids == ["g1"]
    assert not any(ch.action == "reparent-to-mrna" for ch in report.applied)
    assert any(ch.action == "needs-manual" for ch in report.unresolved)


def test_no_transcript_left_for_wrap():
    gff = HDR + (
        "c\tx\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "c\tx\tCDS\t100\t500\t.\t+\t0\tID=g1.cds1;Parent=g1\n"
        "c\tx\texon\t100\t500\t.\t+\t.\tID=g1.ex1;Parent=g1\n"
    )
    doc, report = _norm(gff)
    mrnas = [f for f in doc.features if f.type == "mRNA"]
    assert len(mrnas) == 1
    assert {c.id for c in mrnas[0].children} == {"g1.cds1", "g1.ex1"}
    assert any(ch.action == "add-child-feature" for ch in report.applied)
    assert not any(ch.action == "reparent-to-mrna" for ch in report.applied)


def test_sole_trna_transcript_not_reparented():
    gff = HDR + (
        "c\tx\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "c\tx\ttRNA\t100\t500\t.\t+\t.\tID=g1.t;Parent=g1\n"
        "c\tx\texon\t100\t500\t.\t+\t.\tID=g1.ex1;Parent=g1\n"
    )
    doc, report = _norm(gff)
    assert _feat(doc, "g1.ex1").parent_ids == ["g1"]
    assert not any(ch.action == "reparent-to-mrna" for ch in report.applied)
    assert any(ch.action == "needs-manual" for ch in report.unresolved)


def test_multiexon_misparent_reparented_span_covers():
    gff = HDR + (
        "c\tx\tgene\t100\t900\t.\t+\t.\tID=g1\n"
        "c\tx\tmRNA\t100\t900\t.\t+\t.\tID=g1.1;Parent=g1\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g1.c1;Parent=g1\n"
        "c\tx\tCDS\t800\t900\t.\t+\t2\tID=g1.c2;Parent=g1\n"
        "c\tx\texon\t100\t200\t.\t+\t.\tID=g1.e1;Parent=g1\n"
        "c\tx\texon\t800\t900\t.\t+\t.\tID=g1.e2;Parent=g1\n"
    )
    doc, report = _norm(gff)
    mrna = _feat(doc, "g1.1")
    assert {c.id for c in mrna.children} == {"g1.c1", "g1.c2", "g1.e1", "g1.e2"}
    assert mrna.spans[0].start == 100 and mrna.spans[0].end == 900
    for cid in ("g1.c1", "g1.c2", "g1.e1", "g1.e2"):
        assert _feat(doc, cid).parent_ids == ["g1.1"]


def test_flag_off_no_reparent():
    gff = HDR + (
        "c\tx\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "c\tx\tmRNA\t100\t500\t.\t+\t.\tID=g1.1;Parent=g1\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g1.cds1;Parent=g1\n"
    )
    doc, report = _norm(gff, reparent_gene_children=False)
    assert _feat(doc, "g1.cds1").parent_ids == ["g1"]
    assert not any(ch.action == "reparent-to-mrna" for ch in report.applied)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_normalize_reparent.py -v'`
Expected: the 3 config tests PASS; the 7 new behaviour tests FAIL — `test_misparented_cds_exon_reparented_to_mrna`, `test_two_mrnas_gene_level_cds_not_reparented`, `test_sole_trna_transcript_not_reparented`, `test_multiexon_misparent_reparented_span_covers` fail on assertions (CDS still parented to `g1`, no `reparent-to-mrna` change); `test_wellformed_three_level_is_noop`, `test_no_transcript_left_for_wrap`, `test_flag_off_no_reparent` may already pass (no-op cases) — that is fine.

- [ ] **Step 3: Implement the pass**

In `src/ddbj_gff/normalize/passes.py`, add (place it immediately before `pass_wrap_cds_in_mrna`):

```python
_REPARENT_STRUCTURAL = {"CDS", "exon", "intron", "three_prime_UTR",
                        "five_prime_UTR", "start_codon", "stop_codon"}


def pass_reparent_gene_children_to_mrna(doc, ctx) -> list:
    """Reparent a gene's direct structural sub-features onto the gene's mRNA.

    Fixes AUGUSTUS-dialect GFFs where CDS/exon/etc. are parented to the gene
    instead of the mRNA (so the mRNA is empty and gff2mss would drop the gene).
    Complement of pass_wrap_cds_in_mrna, which handles genes with structural
    children and NO transcript. Acts only when the reparent target is
    unambiguous: exactly one transcript child that is an mRNA.
    """
    changes: list = []
    if not getattr(ctx.config, "reparent_gene_children", True):
        return changes
    for gene in doc.features:
        if gene.type != "gene":
            continue
        structural = [c for c in gene.children if c.type in _REPARENT_STRUCTURAL]
        if not structural:
            continue
        transcripts = [c for c in gene.children if c.type not in _REPARENT_STRUCTURAL]
        if not transcripts:
            continue  # no transcript at all -> pass_wrap_cds_in_mrna's job
        if len(transcripts) != 1 or transcripts[0].type != "mRNA":
            reason = ("multiple transcript children" if len(transcripts) > 1
                      else f"sole transcript is {transcripts[0].type!r}, not mRNA")
            changes.append(Change("needs-manual", gene.id or "?",
                                  f"gene {gene.id!r} has {len(structural)} gene-level "
                                  f"sub-feature(s) but {reason}; left for manual review"))
            continue
        mrna = transcripts[0]
        struct_ids = {id(c) for c in structural}
        for c in structural:
            c.parent_ids = [mrna.id]
            c.parents = [mrna]
            if "Parent" in c.attributes:
                c.attributes["Parent"] = [mrna.id]
            mrna.children.append(c)
        gene.children = [c for c in gene.children if id(c) not in struct_ids]
        seqid = mrna.spans[0].seqid if mrna.spans else structural[0].spans[0].seqid
        strand = mrna.spans[0].strand if mrna.spans else structural[0].spans[0].strand
        lo = min([s.start for s in mrna.spans] + [s.start for c in structural for s in c.spans])
        hi = max([s.end for s in mrna.spans] + [s.end for c in structural for s in c.spans])
        mrna.spans = [Span(seqid, lo, hi, strand)]
        changes.append(Change("reparent-to-mrna", mrna.id or "?",
                              f"reparented {len(structural)} gene-level sub-feature(s) of "
                              f"gene {gene.id!r} to mRNA {mrna.id!r}"))
    return changes
```

- [ ] **Step 4: Register the pass**

In `src/ddbj_gff/normalize/normalize.py`, add `pass_reparent_gene_children_to_mrna` to the import from `.passes`, insert it into `ALL_PASSES` between `pass_coerce_transcript_to_mrna` and `pass_wrap_cds_in_mrna`, and add `"reparent-to-mrna"` to `_APPLIED`:

```python
from .passes import (NormalizeContext, pass_directives, pass_coerce_transcript_to_mrna,
                     pass_reparent_gene_children_to_mrna,
                     pass_wrap_cds_in_mrna, pass_merge_overlapping_loci, pass_circular_origin,
                     pass_trans_splicing_location, pass_so_terms, pass_transl_except, pass_anticodon)
from .report import NormalizationReport

ALL_PASSES = [pass_directives, pass_coerce_transcript_to_mrna,
              pass_reparent_gene_children_to_mrna, pass_wrap_cds_in_mrna,
              pass_merge_overlapping_loci, pass_circular_origin, pass_trans_splicing_location,
              pass_so_terms, pass_transl_except, pass_anticodon]

# actions that represent a clean applied change; everything else needs human attention
_APPLIED = {"add-directive", "rename-type", "add-qualifier", "add-child-feature",
            "merge-loci", "reparent-to-mrna"}
```

- [ ] **Step 5: Run the reparent tests to verify they pass**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_normalize_reparent.py -v'`
Expected: PASS (10 passed — 3 config + 7 behaviour).

- [ ] **Step 6: Run the full not-slow suite to verify no regression**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/ -q -m "not slow"'`
Expected: all pass, 0 failed (the previously-green count plus the 10 new reparent tests). Default-on pass must not change any existing test outcome (well-formed inputs are no-ops).

- [ ] **Step 7: Commit**

```bash
git add src/ddbj_gff/normalize/passes.py src/ddbj_gff/normalize/normalize.py tests/test_normalize_reparent.py
git commit -m "feat(normalize): pass_reparent_gene_children_to_mrna (gene-level CDS/exon -> mRNA, default on)"
```

---

## Notes for the executor

- Branch off `main` first (do not implement on `main` directly). Suggested branch: `feat/reparent-gene-children`.
- This is a `gff_submission`-only change; do not touch `ddbj_mss_tools`.
- Do not `git push`; local commits only.
- The dev container `ddbj-gff-dev` bind-mounts the repo `src/` and `tests/` at `/workspace`, so edits on the host are live inside the container — no copy step needed.
