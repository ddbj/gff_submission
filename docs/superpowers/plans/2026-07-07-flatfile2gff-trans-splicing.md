# flatfile2gff trans-splicing extension — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `flatfile2gff` (`ddbj-gff`) so a trans-spliced, intron-bearing DDBJ flatfile converts to correct canonical INSDC GFF, making the flatfile ⇄ GFF round-trip work for trans-splicing (verified on marchantia chloroplast rps12).

**Architecture:** All changes in `ddbj-gff`'s `src/ddbj_gff/flatfile/convert.py`. A `/trans_splicing` flatfile feature is marked canonically (`exception=trans-splicing` + `is_ordered=true` + per-span `part=`), so the existing `normalize.pass_trans_splicing_location` builds `location=join(...)`. Intron features are emitted (forward path). The trans-spliced gene keeps its segments (multi-part). The reverse chain (`gff2mss`→`mss2ff`) already handles trans-splicing + introns, so round-trip verification is in `ddbj_mss_tools`.

**Tech Stack:** Python 3.10+, BioPython 1.87, pytest 9. All tests run inside the `ddbj-gff-dev` Docker container.

## Global Constraints

- **Two repos, one-way dependency:** changes are in `ddbj-gff` (`gff_submission`), which MUST NOT import `gff2mss`/`common`. The round-trip test (Task 4) lives in `ddbj_mss_tools` and may import `ddbj_gff`.
- **Canonical trans-splicing marking:** a flatfile feature with a `/trans_splicing` qualifier (BioPython: `"trans_splicing" in feature.qualifiers`) → canonical GFF feature gets `exception=trans-splicing` + `is_ordered=true` + `part=1,2,…` on its spans **in BioPython `location.parts` order** (biological 5′→3′). The raw `trans_splicing` attribute is NOT emitted.
- **Do not build `location=` in flatfile2gff** — it is built by the existing `normalize.pass_trans_splicing_location` (which flatfile_to_gff already runs). flatfile2gff only sets `exception`/`is_ordered`/`part`.
- **Intron features** are emitted as children of their locus_tag's gene (`Parent=gene-<locus_tag>`), with `/gene /locus_tag /number /note`, plus trans marking when `/trans_splicing`.
- **Trans-spliced gene = segment-preserving:** its spans = the trans transcript's part-spans (`part=`/`is_ordered=true`), NOT a min..max single span. No `exception` on the gene (so no `location=` is built for it — matches the chloroplast.gff3 canonical form).
- **Verified values (rps12 fixture, pin in tests):** trans CDS `join(complement(1641..1754),93..324,829..854)` `/trans_splicing`, translation `MPTIQQLIRNKRQPIENRTKSPALKGCPQRRGVCTRVYTTTPKKPNSALRKIARVRLTSGFEITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVRYHIIRGTLDAVGVKDRQQGRSKYGVKKSK` (0 internal stops, M-start); trans intron `join(complement(855..1640),1..92)` `/trans_splicing` `/number=1`; cis intron `325..828` `/number=2`. transl_table 11.
- **Commit policy:** `gff_submission` (Tasks 1-3) = normal commits on branch `feat/flatfile2gff-trans-splicing`. `ddbj_mss_tools` (Task 4) = commit ONLY `tests/**`; NEVER `git add -A`; do NOT push.
- **Reference spec:** `docs/superpowers/specs/2026-07-07-flatfile2gff-trans-splicing-design.md`.
- **Fixture (already created, committed as branch prereq):** `tests/flatfile_fixtures/trans_splicing_rps12.gbk` — 1754 bp organelle excerpt; source (`/organelle=plastid:chloroplast`, `taxon:1480154`, circular) + trans CDS + trans intron1 + cis intron2, all `locus_tag=Mp_Cg00010`.

## Test Environment (`ddbj-gff-dev` container)

- `ddbj_gff` bind-mounted live at `/workspace/src`; tests at `/workspace/tests`; fixture at `/workspace/tests/flatfile_fixtures/trans_splicing_rps12.gbk` (Tasks 1-3 need NO sync).
- Task 4: re-sync BOTH `docker cp .../src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss` AND `.../src/common/. ...:/opt/mss_src/common`; copy the fixture + test to `/opt/mss_tests`.
- venv `/opt/ddbj-venv`.

**Run ddbj-gff tests (Tasks 1-3):**
```bash
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev \
  bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/<file> -v'
```

## File Structure

- `src/ddbj_gff/flatfile/convert.py` — add `_is_trans`, `_mark_trans_spliced`; extend `qualifiers_to_attrs`, `synthesize_features` (trans marking + intron emission + segment-preserving gene).
- Tests: `tests/test_flatfile_trans_splicing.py` (new); `ddbj_mss_tools/tests/test_flatfile_trans_roundtrip.py` (new, Task 4).

---

### Task 1: trans-splicing helpers + qualifier drop

**Files:**
- Modify: `src/ddbj_gff/flatfile/convert.py` (add `_is_trans`, `_mark_trans_spliced`; add `trans_splicing` to `_DROP_QUALS`)
- Test: `tests/test_flatfile_trans_splicing.py`

**Interfaces:**
- Produces: `_is_trans(feature) -> bool` (`"trans_splicing" in feature.qualifiers`); `_mark_trans_spliced(attrs: dict, spans: list) -> None` (sets `attrs["exception"]=["trans-splicing"]`, `attrs["is_ordered"]=["true"]`, and `spans[i].part = i+1`). `qualifiers_to_attrs` no longer emits a `trans_splicing` attribute.

- [ ] **Step 1: Write the failing test**

Create `tests/test_flatfile_trans_splicing.py`:
```python
from Bio import SeqIO
from ddbj_gff.model import Span
from ddbj_gff.flatfile.convert import _is_trans, _mark_trans_spliced, qualifiers_to_attrs

FIX = "tests/flatfile_fixtures/trans_splicing_rps12.gbk"


def _cds(rec):
    return next(f for f in rec.features if f.type == "CDS")


def test_is_trans_detects_qualifier():
    rec = SeqIO.read(FIX, "genbank")
    assert _is_trans(_cds(rec)) is True
    src = next(f for f in rec.features if f.type == "source")
    assert _is_trans(src) is False


def test_mark_trans_spliced_sets_exception_ordered_part():
    attrs = {"locus_tag": ["Mp_Cg00010"]}
    spans = [Span("c", 1641, 1754, "-"), Span("c", 93, 324, "+"), Span("c", 829, 854, "+")]
    _mark_trans_spliced(attrs, spans)
    assert attrs["exception"] == ["trans-splicing"]
    assert attrs["is_ordered"] == ["true"]
    assert [s.part for s in spans] == [1, 2, 3]


def test_qualifiers_to_attrs_drops_raw_trans_splicing():
    rec = SeqIO.read(FIX, "genbank")
    attrs = qualifiers_to_attrs(_cds(rec))
    assert "trans_splicing" not in attrs   # replaced by exception (set in synthesize_features)
    assert attrs["product"] == ["ribosomal protein S12"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_flatfile_trans_splicing.py -v'`
Expected: FAIL — `ImportError: cannot import name '_is_trans'`.

- [ ] **Step 3: Implement**

In `src/ddbj_gff/flatfile/convert.py`, change `_DROP_QUALS`:
```python
_DROP_QUALS = {"translation", "codon_start", "trans_splicing"}
```
Add the two helpers (near `_locus_tag`):
```python
def _is_trans(feature) -> bool:
    """True if the flatfile feature carries a /trans_splicing qualifier."""
    return "trans_splicing" in feature.qualifiers


def _mark_trans_spliced(attrs: dict, spans: list) -> None:
    """Mark a synthesized feature as trans-spliced in canonical form: set
    exception=trans-splicing + is_ordered=true and number the spans part=1,2,...
    in their given (biological 5'->3') order. normalize.pass_trans_splicing_location
    then builds location=join(...) from these."""
    attrs["exception"] = ["trans-splicing"]
    attrs["is_ordered"] = ["true"]
    for i, s in enumerate(spans, 1):
        s.part = i
```

- [ ] **Step 4: Run test to verify it passes**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_flatfile_trans_splicing.py -v'`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/gff_submission
git add src/ddbj_gff/flatfile/convert.py tests/test_flatfile_trans_splicing.py
git commit -m "feat(flatfile): trans-splicing helpers (_is_trans, _mark_trans_spliced); drop raw trans_splicing attr"
```

---

### Task 2: synthesize_features — trans marking + intron emission + segment-preserving gene

**Files:**
- Modify: `src/ddbj_gff/flatfile/convert.py` (`synthesize_features`)
- Test: `tests/test_flatfile_trans_splicing.py` (extend)

**Interfaces:**
- Consumes: `_is_trans`, `_mark_trans_spliced` (Task 1), `bio_location_to_spans`, `qualifiers_to_attrs`, `Feature`, `Span`.
- Produces: `synthesize_features` now (a) collects `intron` features and emits them as gene children (trans-marked if `/trans_splicing`); (b) marks a trans-spliced CDS and its synthesized mRNA; (c) emits a trans-spliced gene as segment-preserving multi-part.

- [ ] **Step 1: Write the failing test (extend)**

Append to `tests/test_flatfile_trans_splicing.py`:
```python
from collections import Counter
from ddbj_gff.flatfile.convert import synthesize_features


def test_synthesis_marks_trans_and_emits_introns():
    rec = SeqIO.read(FIX, "genbank")
    feats = synthesize_features(rec, "AP025455")
    counts = Counter(f.type for f in feats)
    assert counts["gene"] == 1 and counts["CDS"] == 1 and counts["intron"] == 2

    cds = next(f for f in feats if f.type == "CDS")
    assert cds.attributes.get("exception") == ["trans-splicing"]
    assert cds.attributes.get("is_ordered") == ["true"]
    assert [s.part for s in cds.spans] == [1, 2, 3]
    assert [s.strand for s in cds.spans] == ["-", "+", "+"]   # per-part strand preserved
    assert "trans_splicing" not in cds.attributes             # raw attr dropped

    gene = next(f for f in feats if f.type == "gene")
    assert len(gene.spans) == 3 and [s.part for s in gene.spans] == [1, 2, 3]   # segment-preserving

    introns = [f for f in feats if f.type == "intron"]
    trans_i = next(f for f in introns if len(f.spans) == 2)
    cis_i = next(f for f in introns if len(f.spans) == 1)
    assert trans_i.attributes.get("exception") == ["trans-splicing"]
    assert trans_i.attributes.get("number") == ["1"]
    assert cis_i.attributes.get("number") == ["2"] and "exception" not in cis_i.attributes
    assert all(i.parent_ids and i.parent_ids[0].startswith("gene-") for i in introns)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_flatfile_trans_splicing.py::test_synthesis_marks_trans_and_emits_introns -v'`
Expected: FAIL — no intron features and CDS lacks `exception` (current behavior emits `trans_splicing=`, no introns).

- [ ] **Step 3: Implement**

In `synthesize_features` (`src/ddbj_gff/flatfile/convert.py`):

(a) Include introns in the collection and group:
```python
    bio = [f for f in rec.features
           if f.type in _RNA_TYPES or f.type == "CDS" or f.type == "intron"]
```

(b) Inside the per-group loop, after `rnas = [...]`, add:
```python
        introns = [f for f in members if f.type == "intron"]
        group_trans = any(_is_trans(m) for m in members)
```
and fix biotype to ignore introns — replace the `biotype = _BIOTYPE.get(members[0].type, "other")` line with:
```python
        biotype_src = next((m.type for m in members if m.type != "intron"), members[0].type)
        biotype = _BIOTYPE.get(biotype_src, "other")
```

(c) Replace the single-span gene emission
```python
        out.append(Feature(gene_id, "DDBJ", "gene", [Span(seqid, g_lo, g_hi, g_strand)],
                           gene_attrs, []))
```
with segment-preserving-when-trans:
```python
        if group_trans:
            trans_src = next(m for m in members if _is_trans(m))
            gene_spans = bio_location_to_spans(trans_src.location, seqid,
                                               is_cds=(trans_src.type == "CDS"))
            gene_attrs["is_ordered"] = ["true"]
            for gi, s in enumerate(gene_spans, 1):
                s.part = gi
        else:
            gene_spans = [Span(seqid, g_lo, g_hi, g_strand)]
        out.append(Feature(gene_id, "DDBJ", "gene", gene_spans, gene_attrs, []))
```

(d) Mark a trans synth-mRNA — in the transcript loop's synth branch (`else:  # synth mRNA from the CDS`), after building `tx_spans`/`tx_attrs` and before setting `tx_attrs["ID"]`, add:
```python
                if _is_trans(member_cds[0]):
                    _mark_trans_spliced(tx_attrs, tx_spans)
```

(e) Mark a trans CDS — in the `for c in member_cds:` loop, after `c_attrs = qualifiers_to_attrs(c)`, add:
```python
                if _is_trans(c):
                    _mark_trans_spliced(c_attrs, cspans)
```

(f) Emit introns — at the END of the per-group loop (after the transcript `for` loop, still inside `for lt, members`), add:
```python
        for k, intr in enumerate(introns, 1):
            ispans = bio_location_to_spans(intr.location, seqid, is_cds=False)
            iattrs = qualifiers_to_attrs(intr)
            if _is_trans(intr):
                _mark_trans_spliced(iattrs, ispans)
            iid = f"intron-{lt}-{k}"
            iattrs["ID"] = [iid]
            iattrs["Parent"] = [gene_id]
            out.append(Feature(iid, "DDBJ", "intron", ispans, iattrs, [gene_id]))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_flatfile_trans_splicing.py -v'`
Expected: PASS (all). If the nuclear Citrus synthesis test (`tests/test_flatfile_synthesis.py`) breaks, check that `introns` collection did not disturb the nuclear path (Citrus has no introns, so `introns==[]` and `group_trans==False` → unchanged behavior).

- [ ] **Step 5: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/gff_submission
git add src/ddbj_gff/flatfile/convert.py tests/test_flatfile_trans_splicing.py
git commit -m "feat(flatfile): mark trans-spliced CDS/mRNA/gene + emit intron features"
```

---

### Task 3: end-to-end flatfile_to_gff (location= built, validate-clean) + regression

**Files:**
- Test: `tests/test_flatfile_trans_splicing.py` (extend)

**Interfaces:**
- Consumes: `flatfile_to_gff` (existing), `ddbj_gff.writer.write`, `ddbj_gff.validate.validate`. `flatfile_to_gff` already runs `normalize`, which now builds `location=` for the trans-marked features.

- [ ] **Step 1: Write the failing/e2e test (extend)**

Append to `tests/test_flatfile_trans_splicing.py`:
```python
from ddbj_gff.flatfile.convert import flatfile_to_gff
from ddbj_gff.writer import write
from ddbj_gff.validate import validate


def test_flatfile_to_gff_builds_location_and_validates():
    rec = SeqIO.read(FIX, "genbank")
    doc = flatfile_to_gff(rec)
    cds = next(f for f in doc.features if f.type == "CDS")
    intron_trans = next(f for f in doc.features if f.type == "intron" and len(f.spans) == 2)
    # normalize built the canonical location= from the part-ordered per-part-strand spans
    assert cds.attributes["location"] == ["join(complement(1641..1754),93..324,829..854)"]
    assert intron_trans.attributes["location"] == ["join(complement(855..1640),1..92)"]
    # no ERROR-level validate diagnostics (region -> feature-type-not-insdc WARNING is allowed)
    errors = [d for d in validate(doc) if getattr(d, "severity", None) and d.severity.name == "ERROR"]
    assert errors == [], f"unexpected validate errors: {[d.code for d in errors]}"
    # trans-splicing no longer flagged noncanonical (location= present)
    assert "noncanonical-special-case" not in {d.code for d in validate(doc)}
```

- [ ] **Step 2: Run test to verify it passes (after Tasks 1-2 it should)**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_flatfile_trans_splicing.py -v'`
Expected: PASS. If `cds.attributes["location"]` is missing, `pass_trans_splicing_location` did not fire — check the CDS has `exception=["trans-splicing"]` and ≥2 spans (Task 2). If validate reports an ERROR, read the code (e.g. `cds-invalid-phase`, `dangling-parent`) and fix the synthesis, not the test.

- [ ] **Step 3: Full ddbj-gff regression (not-slow)**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests -m "not slow" -q'`
Expected: PASS (prior 167 + the new trans-splicing test file). Zero failures. In particular `tests/test_flatfile_synthesis.py` and `tests/test_flatfile_convert.py` (nuclear Citrus) stay green.

- [ ] **Step 4: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/gff_submission
git add tests/test_flatfile_trans_splicing.py
git commit -m "test(flatfile): trans-splicing flatfile_to_gff builds location= and validates clean"
```

---

### Task 4: round-trip verification (ddbj_mss_tools)

**Files:**
- Create: `ddbj_mss_tools/tests/test_flatfile_trans_roundtrip.py`
- Create: `ddbj_mss_tools/tests/flatfile_fixtures/trans_splicing_rps12.gbk` (copy)

**Interfaces:**
- Consumes: `ddbj_gff.flatfile.flatfile_to_gff`, `ddbj_gff.writer.write`, `ddbj_gff.parse`, `gff2mss.convert.convert`, `gff2mss.emit.emit_ann`, `gff2mss.config.MssConfig`, `mss2ff.cli`.
- Produces: proof the trans-spliced CDS + introns survive `flatfile → GFF → gff2mss → mss2ff` (organelle, emit_mrna=false).

- [ ] **Step 1: Copy the fixture**

```bash
mkdir -p /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/flatfile_fixtures
cp /Users/tanizawa/projects/ddbj/gff_submission/tests/flatfile_fixtures/trans_splicing_rps12.gbk \
   /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/flatfile_fixtures/trans_splicing_rps12.gbk
```

- [ ] **Step 2: Write the round-trip test**

Create `ddbj_mss_tools/tests/test_flatfile_trans_roundtrip.py`:
```python
import os, tempfile
from Bio import SeqIO
from ddbj_gff import parse
from ddbj_gff.flatfile import flatfile_to_gff
from ddbj_gff.writer import write
from gff2mss.convert import convert
from gff2mss.emit import emit_ann
from gff2mss.config import MssConfig
import mss2ff.cli

FIX = os.path.join(os.path.dirname(__file__), "flatfile_fixtures", "trans_splicing_rps12.gbk")
PROT = ("MPTIQQLIRNKRQPIENRTKSPALKGCPQRRGVCTRVYTTTPKKPNSALRKIARVRLTSGF"
        "EITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVRYHIIRGTLDAVGVKDRQQGRSKYGVKKSK")


def test_trans_splicing_flatfile_roundtrip(tmp_path):
    rec = SeqIO.read(FIX, "genbank")
    # forward: flatfile -> canonical GFF -> gff2mss (.ann), organelle emit_mrna=false
    doc = parse(write(flatfile_to_gff(rec)))
    cfg = MssConfig(source={}, transl_table=11, product_default="hypothetical protein")
    cfg.emit_mrna = False
    mss_doc, _ = convert(doc, {rec.id: rec.seq}, cfg, common_rows=[])
    ann = emit_ann(mss_doc)

    # the .ann carries the trans CDS (join complement) + /trans_splicing + both introns
    assert "join(complement(1641..1754),93..324,829..854)" in ann
    assert "join(complement(855..1640),1..92)" in ann
    assert "325..828" in ann                     # cis intron
    assert "\ttrans_splicing\t" in ann
    assert ann.count("\tintron\t") == 2

    # reverse: .ann + fasta -> mss2ff -> flatfile'; CDS translation matches original
    ann_p = tmp_path / "o.ann"; ann_p.write_text(ann, encoding="utf-8")
    fa_p = tmp_path / "o.fasta"
    fa_p.write_text(f">{rec.id}\n{str(rec.seq)}\n", encoding="utf-8")
    ff_p = tmp_path / "o.ff"
    mss2ff.cli.main([str(ann_p), str(fa_p), "-o", str(ff_p)])
    rt = SeqIO.read(str(ff_p), "genbank")
    cds = [f for f in rt.features if f.type == "CDS"]
    assert len(cds) == 1
    prot = str(cds[0].extract(rt.seq).translate(table=11))
    prot = prot[:-1] if prot.endswith("*") else prot
    assert prot == PROT
    assert "trans_splicing" in cds[0].qualifiers
    assert sum(1 for f in rt.features if f.type == "intron") == 2
```

- [ ] **Step 3: Sync + run**

Run:
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/common/. ddbj-gff-dev:/opt/mss_src/common
docker exec ddbj-gff-dev bash -lc 'mkdir -p /opt/mss_tests/flatfile_fixtures'
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/flatfile_fixtures/trans_splicing_rps12.gbk ddbj-gff-dev:/opt/mss_tests/flatfile_fixtures/trans_splicing_rps12.gbk
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/test_flatfile_trans_roundtrip.py ddbj-gff-dev:/opt/mss_tests/test_flatfile_trans_roundtrip.py
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_flatfile_trans_roundtrip.py -v'
```
Expected: PASS. If `\ttrans_splicing\t` isn't found, inspect the emitted `.ann` row for the valueless qualifier form and adjust the substring (keep the intent). If the CDS translation mismatches, the phase/part-order in the trans marking is wrong (Task 1/2) — fix there, not the test.

- [ ] **Step 4: Full gff2mss regression**

Run:
```bash
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest . -q'
```
Expected: PASS (prior suite + new test; the pre-existing `@slow` skip). Zero failures.

- [ ] **Step 5: Commit (test + fixture only)**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add tests/test_flatfile_trans_roundtrip.py tests/flatfile_fixtures/trans_splicing_rps12.gbk
git commit -m "test(flatfile): round-trip trans-splicing flatfile -> GFF -> gff2mss -> mss2ff"
```

---

## Self-Review

**Spec coverage:** A (trans marking: exception+is_ordered+part, drop raw trans_splicing) → Tasks 1-2. B (intron emission) → Task 2 (f). C (segment-preserving trans gene) → Task 2 (c). D (organelle emit_mrna=false round-trip) → Task 4. `location=` built by existing normalize (not flatfile2gff) → Task 3 asserts it. Verified values (CDS/intron locations, protein) → Global Constraints, asserted in Tasks 2-4. Regression (nuclear Citrus unaffected) → Task 3 Step 3.

**Placeholder scan:** No TBD/TODO. Every code step is complete. The Task 3/4 contingency notes point at concrete diagnostics/values to check.

**Type consistency:** `_is_trans(feature) -> bool` and `_mark_trans_spliced(attrs, spans) -> None` (Task 1) used in Task 2. `Span.part` is a mutable dataclass field (model.py). `bio_location_to_spans(location, seqid, *, is_cds, codon_start=1)` and `qualifiers_to_attrs(feature)` unchanged from existing. `flatfile_to_gff(rec) -> GffDocument` unchanged. `MssConfig(source={}, transl_table=11, product_default=…)` + `cfg.emit_mrna=False` + `convert(doc, seqs, cfg, common_rows=[])` + `emit_ann` + `mss2ff.cli.main([ann, fasta, "-o", out])` match established usage.

**Non-goals (documented):** origin-spanning + trans combined; NCBI's separate gene_biotype=other intron-holder gene; multi-record WGS flatfiles.
```
