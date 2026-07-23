# Trans-splicing Support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Canonicalize trans-spliced features into a `location=join(...)` attribute in `ddbj-gff` normalize, and make `gff2mss` emit the correct trans-spliced CDS (location + per-part-strand translation + `/trans_splicing`) and intron features.

**Architecture:** Two repos, one-way dependency (`gff2mss → ddbj-gff`). `ddbj-gff` gains a normalize pass `pass_trans_splicing_location` that builds the INSDC `location=` string from a feature's part-ordered, per-part-strand spans. `gff2mss` honors that `location=` verbatim for the MSS location, translates a trans-spliced CDS via a per-part-strand `CompoundLocation`, and emits intron features (currently skipped).

**Tech Stack:** Python 3.10+, BioPython 1.87 (`_insdc_location_string`, `CompoundLocation`, `FeatureLocation`), pydantic 2, pytest 9. All tests run inside the `ddbj-gff-dev` Docker container.

## Global Constraints

- **Two repos, one-way dependency:** `gff2mss` imports `ddbj_gff`; `ddbj_gff` MUST NOT import `gff2mss`/`common`.
- **Per-part strand mapping (the rule that makes BioPython emit the right INSDC string):** GFF strand `-` → BioPython `strand=-1` (emitted as a per-part `complement(...)`); GFF strand `+`, `?`, `.` → `strand=+1`. Build the `CompoundLocation` from `feature.ordered_spans()` (part order, honoring `is_ordered`) — do NOT sort by coordinate.
- **Verified location strings (exact, from the fixture — do not alter):** for the rps12 fixture (seqlen 1754):
  - trans CDS → `join(complement(1641..1754),93..324,829..854)`
  - trans intron 1 → `join(complement(855..1640),1..92)`
  - cis intron 2 → `325..828`
- **Verified translation (exact):** the trans CDS translates (table 11, codon_start 1) to
  `MPTIQQLIRNKRQPIENRTKSPALKGCPQRRGVCTRVYTTTPKKPNSALRKIARVRLTSGFEITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVRYHIIRGTLDAVGVKDRQQGRSKYGVKKSK`
  (123 aa + terminal stop; M-start; 0 internal stops).
- **Preserve existing `location=`:** if a feature already has a `location=` attribute (authoritative; may contain remote seqids), do NOT overwrite it.
- **Commit policy:**
  - `gff_submission` (Task 1): branch `feat/trans-splicing` off `main` (@ current HEAD, includes the circular-origin work). Normal commits of the specific files.
  - `ddbj_mss_tools` (Tasks 2–4): branch `feat/trans-splicing` off `main`. Commit ONLY `src/gff2mss/**` and `tests/**` files you create/modify. NEVER `git add -A` — ~49 pre-existing uncommitted files in `examples/`/`docs/`/`data/` must stay untouched. Do NOT push.
- **Reference spec:** `docs/superpowers/specs/2026-07-05-trans-splicing-design.md`.
- **Fixture (already created, in repo):** `gff_submission/tests/normalize_fixtures/trans_splicing_rps12.{gff3,fasta}` — a re-coordinatized excerpt of *Marchantia polymorpha* chloroplast `AP025455.1` (real 1..854 ++ 65903..66802 = 1754 bp), rps12 trans-spliced gene. Parsed shape (verified): CDS `is_trans_spliced=True`, 3 spans parts `[1,2,3]` strands `[-,+,+]` phases `[0,0,2]`; intron `id-Mp_Cg00010` trans, 2 spans; intron `id-Mp_Cg00010-2` cis, 1 span; genes not trans-spliced.

## Test Environment (`ddbj-gff-dev` container)

- `ddbj_gff` source bind-mounted live at `/workspace/src`; tests at `/workspace/tests` (Task 1 needs NO sync).
- `gff2mss` + `common` at `/opt/mss_src` (docker-cp'd — re-sync after every edit). ddbj_mss_tools tests at `/opt/mss_tests`, fixtures at `/opt/mss_tests/mss_fixtures`.
- venv `/opt/ddbj-venv` (Bio 1.87, pydantic 2.13, pytest 9.1).

**Run ddbj-gff tests (Task 1):**
```bash
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev \
  bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/<file> -v'
```

**Sync + run gff2mss tests (Tasks 2–4):**
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/<file> ddbj-gff-dev:/opt/mss_tests/<file>
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest <file> -v'
```
> `docker cp .../gff2mss/.` (trailing `/.`) copies contents into the existing dir; without it the copy NESTS at `/opt/mss_src/gff2mss/gff2mss`.

## File Structure

- `gff_submission/src/ddbj_gff/normalize/passes.py` — add `pass_trans_splicing_location` (+ Bio imports at module top).
- `gff_submission/src/ddbj_gff/normalize/normalize.py` — register the pass in `ALL_PASSES`.
- `ddbj_mss_tools/src/gff2mss/convert.py` — add `_location_attr` helper; `_build_trans_spliced_cds`; branch in `build_cds_feature`; `build_intron_feature`; intron collection in `build_entry_features`.
- New tests: `gff_submission/tests/test_normalize_trans_splicing.py`; `ddbj_mss_tools/tests/test_mss_trans_splicing.py`; `ddbj_mss_tools/tests/mss_fixtures/trans_splicing_rps12.{gff3,fasta}` (copied).

---

### Task 1: `pass_trans_splicing_location` normalize pass (ddbj-gff)

**Files:**
- Modify: `gff_submission/src/ddbj_gff/normalize/passes.py` (add pass + 2 Bio imports at top)
- Modify: `gff_submission/src/ddbj_gff/normalize/normalize.py` (import + register in `ALL_PASSES` after `pass_circular_origin`)
- Test: `gff_submission/tests/test_normalize_trans_splicing.py`

**Interfaces:**
- Consumes: `Feature.is_trans_spliced`, `Feature.ordered_spans()`, `Feature.attributes`, `Feature.spans` (each `Span` has `.seqid/.start/.end/.strand`), `doc.sequence_regions`, `ctx.seq_lengths`, `Change(action, target, message)` from `.report`.
- Produces: `pass_trans_splicing_location(doc, ctx) -> list[Change]`. For each trans-spliced feature with ≥2 spans and no existing `location=`, sets `f.attributes["location"] = [insdc_string]`. Registered in `ALL_PASSES`.

- [ ] **Step 1: Write the failing test**

Create `gff_submission/tests/test_normalize_trans_splicing.py`:
```python
from ddbj_gff import parse
from ddbj_gff.model import Feature, Span, Directive, GffDocument
from ddbj_gff.normalize.passes import pass_trans_splicing_location, NormalizeContext
from ddbj_gff.normalize.normalize import normalize
from ddbj_gff.validate import validate


def _ctx():
    return NormalizeContext(vocab=None, seq_lengths=None, config=None)


def test_builds_location_for_mixed_strand_trans_cds():
    # 3 parts, part order 1,2,3, strands -,+,+
    spans = [Span("c", 1641, 1754, "-", 0, part=1),
             Span("c", 93, 324, "+", 0, part=2),
             Span("c", 829, 854, "+", 2, part=3)]
    cds = Feature("cds", "S", "CDS", spans, {"exception": ["trans-splicing"]}, [])
    seq_dir = Directive("x", "sequence-region", ("c", 1, 1754))
    doc = GffDocument(directives=[seq_dir], features=[cds])
    changes = pass_trans_splicing_location(doc, _ctx())
    assert cds.attributes["location"] == ["join(complement(1641..1754),93..324,829..854)"]
    assert len(changes) == 1


def test_two_part_trans_intron_and_qmark_strand():
    spans = [Span("c", 855, 1640, "-", None, part=1),
             Span("c", 1, 92, "?", None, part=2)]
    intron = Feature("i", "S", "intron", spans, {"exception": ["trans-splicing"]}, [])
    seq_dir = Directive("x", "sequence-region", ("c", 1, 1754))
    doc = GffDocument(directives=[seq_dir], features=[intron])
    pass_trans_splicing_location(doc, _ctx())
    assert intron.attributes["location"] == ["join(complement(855..1640),1..92)"]


def test_skips_non_trans_and_single_span_and_existing_location():
    non_trans = Feature("g", "S", "gene",
                        [Span("c", 855, 1754, "-", None, part=1), Span("c", 1, 854, "?", None, part=2)],
                        {"is_ordered": ["true"]}, [])           # multi-part but not trans-spliced
    cis = Feature("c1", "S", "intron", [Span("c", 325, 828, "+")], {"exception": ["trans-splicing"]}, [])  # single span
    preset = Feature("p", "S", "CDS",
                     [Span("c", 1641, 1754, "-", 0, part=1), Span("c", 93, 324, "+", 0, part=2)],
                     {"exception": ["trans-splicing"], "location": ["join(remote:1..9,1..2)"]}, [])
    doc = GffDocument(directives=[Directive("x", "sequence-region", ("c", 1, 1754))],
                      features=[non_trans, cis, preset])
    pass_trans_splicing_location(doc, _ctx())
    assert "location" not in non_trans.attributes
    assert "location" not in cis.attributes
    assert preset.attributes["location"] == ["join(remote:1..9,1..2)"]   # preserved


def test_full_normalize_and_validate_clean_on_fixture():
    with open("tests/normalize_fixtures/trans_splicing_rps12.gff3") as fh:
        doc = parse(fh.read())
    work, _report = normalize(doc)
    cds = next(f for f in work.features if f.type == "CDS")
    intron1 = next(f for f in work.features if f.type == "intron" and len(f.spans) == 2)
    assert cds.attributes["location"] == ["join(complement(1641..1754),93..324,829..854)"]
    assert intron1.attributes["location"] == ["join(complement(855..1640),1..92)"]
    # canonical form set -> no noncanonical-special-case for trans-splicing
    codes = {d.code for d in validate(work)}
    assert "noncanonical-special-case" not in codes
```

- [ ] **Step 2: Run test to verify it fails**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_normalize_trans_splicing.py -v'`
Expected: FAIL — `ImportError: cannot import name 'pass_trans_splicing_location'`.

- [ ] **Step 3: Add the Bio imports + the pass**

In `gff_submission/src/ddbj_gff/normalize/passes.py`, add near the top imports (after the existing `from .report import Change`):
```python
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.SeqIO.InsdcIO import _insdc_location_string
```
Append the pass:
```python
def pass_trans_splicing_location(doc, ctx) -> list:
    """Build the canonical INSDC location=join(...) attribute for trans-spliced
    multi-part features. Per-part strand: '-' -> complement, '+'/'?'/'.' -> forward;
    parts kept in feature.ordered_spans() order (honors is_ordered). An existing
    location= is authoritative (may be remote) and is left untouched."""
    changes: list = []
    regions = doc.sequence_regions
    for f in doc.features:
        if not f.is_trans_spliced or len(f.spans) < 2 or f.attributes.get("location"):
            continue
        parts = f.ordered_spans()
        seqlen = regions.get(parts[0].seqid, (None, None))[1]
        if seqlen is None and ctx.seq_lengths:
            seqlen = ctx.seq_lengths.get(parts[0].seqid)
        seqlen = seqlen or max(s.end for s in parts)
        locs = [FeatureLocation(s.start - 1, s.end, strand=(-1 if s.strand == "-" else 1))
                for s in parts]
        loc_str = _insdc_location_string(CompoundLocation(locs), seqlen)
        f.attributes["location"] = [loc_str]
        changes.append(Change("add-qualifier", f.id or "?",
                              f"built location={loc_str} for trans-spliced feature"))
    return changes
```

- [ ] **Step 4: Register the pass**

In `gff_submission/src/ddbj_gff/normalize/normalize.py`, add `pass_trans_splicing_location` to the import from `.passes` and insert it into `ALL_PASSES` immediately after `pass_circular_origin`:
```python
from .passes import (NormalizeContext, pass_directives, pass_coerce_transcript_to_mrna,
                     pass_wrap_cds_in_mrna, pass_circular_origin, pass_trans_splicing_location,
                     pass_so_terms, pass_transl_except, pass_anticodon)

ALL_PASSES = [pass_directives, pass_coerce_transcript_to_mrna, pass_wrap_cds_in_mrna,
              pass_circular_origin, pass_trans_splicing_location, pass_so_terms,
              pass_transl_except, pass_anticodon]
```

- [ ] **Step 5: Run test to verify it passes**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_normalize_trans_splicing.py -v'`
Expected: PASS (4 passed).

- [ ] **Step 6: Full ddbj-gff regression (not-slow)**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests -m "not slow" -q'`
Expected: PASS (prior 154 + 4 new = 158). Zero failures. (The existing `test_roundtrip.py` trans_splicing test asserts structure preservation; the new pass only ADDS a `location=` attribute, which round-trips through the writer — if roundtrip breaks, the pass is being applied where it shouldn't; investigate. Note: `test_roundtrip` parses without normalize, so it is unaffected.)

- [ ] **Step 7: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/gff_submission
git add src/ddbj_gff/normalize/passes.py src/ddbj_gff/normalize/normalize.py tests/test_normalize_trans_splicing.py tests/normalize_fixtures/trans_splicing_rps12.gff3 tests/normalize_fixtures/trans_splicing_rps12.fasta
git commit -m "feat(normalize): pass_trans_splicing_location builds INSDC location=join for trans-spliced features"
```

---

### Task 2: trans-spliced CDS in gff2mss (location + translation + `/trans_splicing`)

**Files:**
- Modify: `ddbj_mss_tools/src/gff2mss/convert.py` (add `_location_attr`, `_build_trans_spliced_cds`; branch in `build_cds_feature`)
- Create: `ddbj_mss_tools/tests/mss_fixtures/trans_splicing_rps12.gff3` + `.fasta` (copies)
- Test: `ddbj_mss_tools/tests/test_mss_trans_splicing.py`

**Interfaces:**
- Consumes: `ddbj_gff.model.{Feature,Span}`, existing `_product`, `_submitter_note`, `MssFeature`, `MssQualifier`, `Diagnostic`, `Severity`, `Seq`, `CompoundLocation`, `FeatureLocation`, `CodonTable`. `Feature.is_trans_spliced`, `Feature.ordered_spans()`, `Feature.transl_table`, `Feature.attributes`, `Feature.gene`.
- Produces: `_location_attr(feature) -> str | None`; `_build_trans_spliced_cds(mrna, gene, locus_tag, genome_seq, cfg, diagnostics, cds_feat) -> MssFeature`. `build_cds_feature` returns a trans-spliced CDS (or `misc_feature` on internal stop) when `cds_feat.is_trans_spliced`.

- [ ] **Step 1: Copy the fixture into ddbj_mss_tools**

```bash
cp /Users/tanizawa/projects/ddbj/gff_submission/tests/normalize_fixtures/trans_splicing_rps12.gff3 \
   /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/mss_fixtures/trans_splicing_rps12.gff3
cp /Users/tanizawa/projects/ddbj/gff_submission/tests/normalize_fixtures/trans_splicing_rps12.fasta \
   /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/mss_fixtures/trans_splicing_rps12.fasta
```

- [ ] **Step 2: Write the failing test**

Create `ddbj_mss_tools/tests/test_mss_trans_splicing.py`:
```python
import os
from Bio import SeqIO
from ddbj_gff import parse
from ddbj_gff.normalize.normalize import normalize
from gff2mss.convert import build_entry_features
from gff2mss.config import MssConfig

FIX = os.path.join(os.path.dirname(__file__), "mss_fixtures")
PROT = ("MPTIQQLIRNKRQPIENRTKSPALKGCPQRRGVCTRVYTTTPKKPNSALRKIARVRLTSGF"
        "EITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVRYHIIRGTLDAVGVKDRQQGRSKYGVKKSK")


def _load():
    with open(os.path.join(FIX, "trans_splicing_rps12.gff3")) as fh:
        doc, _ = normalize(parse(fh.read()))
    seqs = {r.id: r.seq for r in SeqIO.parse(os.path.join(FIX, "trans_splicing_rps12.fasta"), "fasta")}
    cfg = MssConfig(source={}, transl_table=11, product_default="hypothetical protein")
    cfg.emit_mrna = False
    return doc, seqs, cfg


def test_trans_spliced_cds_location_and_translation():
    doc, seqs, cfg = _load()
    feats = build_entry_features(doc, seqs, cfg, [])["AP025455.1"]
    cds = [f for f in feats if f.key == "CDS"]
    assert len(cds) == 1
    assert cds[0].location == "join(complement(1641..1754),93..324,829..854)"
    # clean translation -> CDS (not misc_feature)
    assert not any(f.key == "misc_feature" for f in feats)
    # /trans_splicing qualifier present (valueless)
    assert any(q.key == "trans_splicing" for q in cds[0].qualifiers)
    # correct product + codon_start
    assert any(q.key == "product" and q.value == "ribosomal protein S12" for q in cds[0].qualifiers)
    assert any(q.key == "codon_start" and q.value == "1" for q in cds[0].qualifiers)
```

- [ ] **Step 3: Sync + run to verify failure**

Run:
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/mss_fixtures/trans_splicing_rps12.gff3 ddbj-gff-dev:/opt/mss_tests/mss_fixtures/trans_splicing_rps12.gff3
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/mss_fixtures/trans_splicing_rps12.fasta ddbj-gff-dev:/opt/mss_tests/mss_fixtures/trans_splicing_rps12.fasta
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/test_mss_trans_splicing.py ddbj-gff-dev:/opt/mss_tests/test_mss_trans_splicing.py
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_mss_trans_splicing.py -v'
```
Expected: FAIL — the CDS location is wrong (single-strand `_ordered`/`build_insdc_location` produce something like `complement(join(...))` over all spans) and/or `/trans_splicing` absent.

- [ ] **Step 4: Add `_location_attr` and `_build_trans_spliced_cds`; branch in `build_cds_feature`**

In `ddbj_mss_tools/src/gff2mss/convert.py`, add the helper near `_ordered`:
```python
def _location_attr(feature):
    """Return the verbatim INSDC location string from a location= attribute, or None."""
    v = feature.attributes.get("location")
    return v[0] if v else None


def _trans_compound(parts):
    """CompoundLocation from part-ordered spans, per-part strand ('-' -> complement)."""
    return CompoundLocation([FeatureLocation(s.start - 1, s.end,
                             strand=(-1 if s.strand == "-" else 1)) for s in parts])
```
Add the trans-spliced CDS builder (place it just before `build_cds_feature`):
```python
def _build_trans_spliced_cds(mrna, gene, locus_tag, genome_seq, cfg, diagnostics, cds_feat):
    parts = cds_feat.ordered_spans()
    location = _location_attr(cds_feat)
    if location is None:                       # normalize should have set it; be defensive
        location = _insdc_location_string(_trans_compound(parts), len(genome_seq))
    table_id = cds_feat.transl_table or cfg.transl_table
    codon_start = 1 if parts[0].phase is None else parts[0].phase + 1
    quals = [
        MssQualifier("locus_tag", locus_tag),
        MssQualifier("transl_table", str(table_id)),
        MssQualifier("codon_start", str(codon_start)),
        MssQualifier("product", _product(mrna, gene, cfg)),
    ]
    if gene.gene or mrna.gene:
        quals.append(MssQualifier("gene", gene.gene or mrna.gene))
    quals.append(_submitter_note(gene, mrna))
    quals.append(MssQualifier("trans_splicing", ""))

    entry_seqid = parts[0].seqid
    if any(s.seqid != entry_seqid for s in parts) or ":" in location:
        diagnostics.append(Diagnostic(Severity.WARNING, None, "trans-splicing-remote",
                                      f"CDS {mrna.id!r} trans-splicing references a remote seqid; "
                                      f"emitted without translation validation"))
        return MssFeature("CDS", location, quals)

    coding = str(_trans_compound(parts).extract(genome_seq))[codon_start - 1:].upper()
    coding_full = coding[: len(coding) - len(coding) % 3]
    protein = str(Seq(coding_full).translate(table=table_id))
    body = protein[:-1] if protein.endswith("*") else protein
    if "*" in body:
        diagnostics.append(Diagnostic(Severity.WARNING, None, "translation-internal-stop",
                                      f"trans-spliced CDS {mrna.id!r} has an internal stop codon"))
        note = f"internal stop codon(s) detected in CDS {mrna.id}; not translated"
        return MssFeature("misc_feature", location,
                          [MssQualifier("locus_tag", locus_tag), MssQualifier("note", note),
                           MssQualifier("trans_splicing", "")])
    return MssFeature("CDS", location, quals)
```
In `build_cds_feature`, right after `cds_feat = next((c for c in mrna.children if c.type == "CDS"), None)` is computed (it currently appears later — move the `cds_feat` lookup up to just after the `spans`/`if not spans` block, OR add a fresh lookup), insert the branch BEFORE the `_ordered(spans)` line:
```python
    ts_feat = next((c for c in mrna.children if c.type == "CDS" and c.is_trans_spliced), None)
    if ts_feat is not None:
        return _build_trans_spliced_cds(mrna, gene, locus_tag, genome_seq, cfg, diagnostics, ts_feat)
```
(Keep the existing non-trans path untouched below this branch.)

- [ ] **Step 5: Sync + run to verify pass**

Run:
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_mss_trans_splicing.py -v'
```
Expected: PASS.

- [ ] **Step 6: Commit (gff2mss + test + fixtures)**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add src/gff2mss/convert.py tests/test_mss_trans_splicing.py tests/mss_fixtures/trans_splicing_rps12.gff3 tests/mss_fixtures/trans_splicing_rps12.fasta
git commit -m "feat(gff2mss): trans-spliced CDS location + per-part-strand translation + /trans_splicing"
```

---

### Task 3: emit intron features in gff2mss

**Files:**
- Modify: `ddbj_mss_tools/src/gff2mss/convert.py` (add `build_intron_feature`; collect introns in `build_entry_features`)
- Test: `ddbj_mss_tools/tests/test_mss_trans_splicing.py` (extend)

**Interfaces:**
- Consumes: `_location_attr` (Task 2), `build_insdc_location`, `Feature.gene`, `Feature._first`, `Feature.note`, `Feature.is_trans_spliced`, `Feature.spans`, `_span_start`.
- Produces: `build_intron_feature(intron, seqlen) -> MssFeature`; `build_entry_features` now emits `intron` features (in coordinate order with the other per-seqid features).

- [ ] **Step 1: Write the failing test (append)**

Append to `ddbj_mss_tools/tests/test_mss_trans_splicing.py`:
```python
def test_intron_features_emitted():
    doc, seqs, cfg = _load()
    feats = build_entry_features(doc, seqs, cfg, [])["AP025455.1"]
    introns = [f for f in feats if f.key == "intron"]
    locs = {f.location for f in introns}
    assert "join(complement(855..1640),1..92)" in locs   # trans intron 1
    assert "325..828" in locs                             # cis intron 2
    trans_intron = next(f for f in introns if f.location == "join(complement(855..1640),1..92)")
    assert any(q.key == "trans_splicing" for q in trans_intron.qualifiers)
    assert any(q.key == "number" and q.value == "1" for q in trans_intron.qualifiers)
    assert any(q.key == "gene" and q.value == "rps12" for q in trans_intron.qualifiers)
    cis_intron = next(f for f in introns if f.location == "325..828")
    assert not any(q.key == "trans_splicing" for q in cis_intron.qualifiers)
```

- [ ] **Step 2: Sync + run to verify failure**

Run:
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/test_mss_trans_splicing.py ddbj-gff-dev:/opt/mss_tests/test_mss_trans_splicing.py
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_mss_trans_splicing.py::test_intron_features_emitted -v'
```
Expected: FAIL — no `intron` features in the output (introns are skipped as `_STRUCTURAL`).

- [ ] **Step 3: Add `build_intron_feature` + collect introns in `build_entry_features`**

In `ddbj_mss_tools/src/gff2mss/convert.py`, add the builder (near `build_noncoding_features`):
```python
def build_intron_feature(intron, seqlen) -> MssFeature:
    location = _location_attr(intron) or build_insdc_location(intron.spans, seqlen)
    quals = []
    if intron.gene:
        quals.append(MssQualifier("gene", intron.gene))
    lt = intron._first("locus_tag")
    if lt:
        quals.append(MssQualifier("locus_tag", lt))
    for n in intron.note:
        quals.append(MssQualifier("note", n))
    num = intron._first("number")
    if num:
        quals.append(MssQualifier("number", num))
    if intron.is_trans_spliced:
        quals.append(MssQualifier("trans_splicing", ""))
    return MssFeature("intron", location, quals)
```
In `build_entry_features`, inside the per-seqid loop, collect introns and add them to the `items` list so they sort by coordinate with the genes/parentless RNA. After the existing `genes`/`parentless` collection and before building `items`, add:
```python
        introns = [f for f in doc.features if f.type == "intron"
                   and any(s.seqid == seqid for s in f.spans)]
```
Change the `items` construction to include introns, and handle them in the emit loop:
```python
        items = ([(_span_start(g), g) for g in genes]
                 + [(_span_start(r), r) for r in parentless]
                 + [(_span_start(i), i) for i in introns])
        items.sort(key=lambda t: t[0])
        feats: list = []
        for _, feat in items:
            if feat.type == "gene":
                feats.extend(build_gene_features(feat, cfg.transcript_mode, assigner,
                                                 genome_seq, cfg, diagnostics))
            elif feat.type == "intron":
                feats.append(build_intron_feature(feat, len(genome_seq)))
            else:
                feats.append(build_rna_feature(feat, assigner.assign(feat),
                                                len(genome_seq), feat.id, feat.id, feat.gene))
```

- [ ] **Step 4: Sync + run to verify pass**

Run:
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_mss_trans_splicing.py -v'
```
Expected: PASS (both trans-splicing tests).

- [ ] **Step 5: Commit (gff2mss + test)**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add src/gff2mss/convert.py tests/test_mss_trans_splicing.py
git commit -m "feat(gff2mss): emit intron features (with /trans_splicing when trans-spliced)"
```

---

### Task 4: end-to-end .ann text assertion + full regression

**Files:**
- Test: `ddbj_mss_tools/tests/test_mss_trans_splicing.py` (extend with an `.ann`-text assertion)

**Interfaces:**
- Consumes: `gff2mss.convert.convert(doc, seqs, cfg, common_rows, *, strict=False) -> (MssDocument, diagnostics)` and `gff2mss.emit.emit_ann(MssDocument) -> str`. Reuses the `_load()` helper (parse + `normalize`) from Task 2.
- Produces: proof the whole chain yields, in the emitted `.ann` text, the trans CDS line + `/trans_splicing`, and both intron lines.

> **Architecture note (verified):** `gff2mss` never runs `normalize` — neither `convert` nor
> `assemble.build_ann_text` does; gff2mss consumes an **already-canonical** GFF (the
> heterosigma pipeline runs a separate `ddbj-gff` normalize step, producing
> `*.normalized.gff3`, before conversion). So `build_ann_text` on the RAW fixture would NOT
> have `location=` and would mis-emit the trans CDS. This task therefore normalizes in memory
> (via `_load()`), then `convert` + `emit_ann` — the same in-memory canonicalization Tasks 2/3
> use. (For real submission the marchantia GFF is normalized first, exactly like heterosigma.)

- [ ] **Step 1: Write the e2e test (append)**

Append to `ddbj_mss_tools/tests/test_mss_trans_splicing.py`:
```python
def test_ann_text_end_to_end():
    from gff2mss.convert import convert
    from gff2mss.emit import emit_ann
    doc, seqs, cfg = _load()                     # parse + normalize (location= built)
    mss_doc, _ = convert(doc, seqs, cfg, common_rows=[])
    text = emit_ann(mss_doc)
    assert "join(complement(1641..1754),93..324,829..854)" in text   # trans CDS
    assert "join(complement(855..1640),1..92)" in text               # trans intron 1
    assert "\ttrans_splicing\t" in text          # valueless qualifier row (empty value col)
    assert "\tintron\t" in text                  # intron feature key column
```

- [ ] **Step 2: Sync + run**

Run:
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/test_mss_trans_splicing.py ddbj-gff-dev:/opt/mss_tests/test_mss_trans_splicing.py
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_mss_trans_splicing.py -v'
```
Expected: PASS (3 tests). If the `\ttrans_splicing\t` assertion fails, inspect the emitted row — a valueless qualifier renders as `\t\t\ttrans_splicing\t` (empty key/location cols on a continuation row and an empty value col); adjust the substring to match the actual emitter output (confirm against `emit.feature_rows`), keeping the intent (a `/trans_splicing` row exists).

- [ ] **Step 3: Full gff2mss regression**

Run:
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest . -q'
```
Expected: PASS. If a pre-existing fixture that contains `intron` rows now emits new `intron` features and breaks an expected-output/snapshot test, update that expected output to include the introns (introns are now intentionally emitted) and note it in the report. heterosigma organelle (0 introns) is unaffected.

- [ ] **Step 4: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add tests/test_mss_trans_splicing.py
# include any updated fixture/expected-output files touched in Step 3 (name them explicitly; never git add -A)
git commit -m "test(gff2mss): end-to-end trans-splicing .ann (CDS + introns + /trans_splicing)"
```

---

## Self-Review

**Spec coverage:**
- A1 `pass_trans_splicing_location` → Task 1 (with `?`/`.`→forward, preserve existing `location=`, skip non-trans/single-span — all tested). A2 validate-clean → Task 1 Step 1 `test_full_normalize_and_validate_clean_on_fixture`. B1 `_location_attr` honor `location=` → Task 2 Step 4 (used by CDS) + Task 3 (used by intron). B2 trans CDS location+translation+`/trans_splicing`+internal-stop→misc_feature+remote-diagnostic → Task 2. B3 intron emission (trans + cis, `gene_biotype=other` covered by collecting introns doc-wide, not per-gene) → Task 3. Verified location strings + protein → Global Constraints, asserted in Tasks 1–4. End-to-end `.ann` → Task 4. Regression (ddbj-gff 154→158, gff2mss suite) → Task 1 Step 6, Task 4 Step 3.
- Non-goals (remote-part translation, AUGUSTUS intron stripping, part-row collapsing, frameshift) → not implemented; remote is diagnosed (Task 2).

**Placeholder scan:** No TBD/TODO. Every code step has complete code. Two contingency notes (Task 4 Step 2 normalize-wiring, Step 3 snapshot update) are tied to concrete files the implementer must confirm — they instruct exactly what to check and do, not vague "handle it".

**Type consistency:** `pass_trans_splicing_location(doc, ctx) -> list[Change]` matches sibling passes and `ALL_PASSES` usage. `_location_attr(feature) -> str|None` and `_trans_compound(parts)` defined in Task 2, reused in Task 3. `Span(seqid, start, end, strand, phase, part=...)` positional/keyword form matches existing tests and the parser. `MssFeature(key, location, quals)` / `MssQualifier(key, value)` match `model.py`. `build_entry_features(doc, seqs, cfg, diagnostics) -> dict[seqid, list[MssFeature]]` matches `convert.py`. Per-part strand mapping (`-1 if s.strand == "-" else 1`) is identical in normalize (Task 1) and gff2mss (`_trans_compound`, Task 2).
```
