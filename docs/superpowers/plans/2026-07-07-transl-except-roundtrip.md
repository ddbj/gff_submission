# transl_except round-trip — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make gff2mss emit the `/transl_except` qualifier (from the recoded_codon/stop_codon children it already collects) so a CDS with a translational exception (e.g. Pyrrolysine read-through) round-trips flatfile ⇄ canonical GFF.

**Architecture:** One small change in `ddbj_mss_tools` `gff2mss/convert.py` `build_cds_feature`: append `MssQualifier("transl_except", spec)` for each already-computed `excepts` spec on the clean-translation CDS path. The trans-spliced + transl_except combination is out of scope (diagnostic only). The forward side (flatfile2gff→normalize→recoded_codon) and translation (`translate_cds_with_transl_except`) already work; only qualifier emission is missing.

**Tech Stack:** Python 3.10+, BioPython 1.87, pytest 9. All tests run inside the `ddbj-gff-dev` Docker container.

## Global Constraints

- **Repo:** all code changes are in `ddbj_mss_tools` (`gff2mss`). `ddbj-gff` is unchanged except the committed test fixture. Commit ONLY `src/gff2mss/**` and `tests/**` in ddbj_mss_tools; NEVER `git add -A` (~49 pre-existing uncommitted files); do NOT push.
- **The fix:** on `build_cds_feature`'s clean-translation path (the one that returns `MssFeature("CDS", location, quals)` at the end), append `MssQualifier("transl_except", spec)` for each `spec` in the already-computed `excepts` list. Do NOT touch the misc_feature (internal-stop) return path.
- **Scope:** non-trans-spliced CDS + transl_except is fully supported (qualifier + translation). A CDS that is BOTH trans-spliced AND has transl_except is out of scope: in `_build_trans_spliced_cds`, if `_collect_transl_excepts(cds_feat)` is non-empty, emit a `transl-except-trans-splicing` WARNING; do not otherwise change that path.
- **Verified fact:** `_collect_transl_excepts` reconstructs `(pos:746..748,aa:Pyl)` exactly from the `recoded_codon` child (aa_names round-trips pyrrolysine→Pyl), so the emitted qualifier matches the original.
- **Reference spec:** `docs/superpowers/specs/2026-07-07-transl-except-roundtrip-design.md`.
- **Fixture (already created, committed as branch prereq in ddbj-gff):** `tests/flatfile_fixtures/transl_except_p87.gbk` — 2400 bp LC757512 excerpt; source + CDS `20..2320` `/transl_except=(pos:746..748,aa:Pyl)` `/transl_table=1` `/gene=p87`. Naive translation has one internal stop at protein index 242.

## Test Environment (`ddbj-gff-dev` container)

- `ddbj_gff` live at `/workspace/src`; the fixture will be copied into ddbj_mss_tools tests.
- gff2mss: re-sync BOTH `docker cp .../src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss` AND `.../src/common/. ...:/opt/mss_src/common` (container `common` goes stale); copy fixture + test to `/opt/mss_tests`.
- venv `/opt/ddbj-venv`. Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest <file> -v'`.

## File Structure

- Modify: `ddbj_mss_tools/src/gff2mss/convert.py` (`build_cds_feature` emit; `_build_trans_spliced_cds` diagnostic).
- Create: `ddbj_mss_tools/tests/flatfile_fixtures/transl_except_p87.gbk` (copy of the ddbj-gff fixture).
- Create: `ddbj_mss_tools/tests/test_transl_except_roundtrip.py`.

---

### Task 1: gff2mss emits `/transl_except` (+ trans-spliced combo diagnostic)

**Files:**
- Modify: `ddbj_mss_tools/src/gff2mss/convert.py`
- Test: `ddbj_mss_tools/tests/test_transl_except_roundtrip.py` (unit part)
- Copy: `ddbj_mss_tools/tests/flatfile_fixtures/transl_except_p87.gbk`

**Interfaces:**
- Consumes: existing `_collect_transl_excepts`, `MssQualifier`, `Diagnostic`, `Severity`.
- Produces: `build_cds_feature` clean-CDS output now includes one `transl_except` qualifier per collected spec; `_build_trans_spliced_cds` emits a `transl-except-trans-splicing` WARNING when a trans-spliced CDS also has transl_except.

- [ ] **Step 1: Copy the fixture into ddbj_mss_tools**

```bash
mkdir -p /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/flatfile_fixtures
cp /Users/tanizawa/projects/ddbj/gff_submission/tests/flatfile_fixtures/transl_except_p87.gbk \
   /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/flatfile_fixtures/transl_except_p87.gbk
```

- [ ] **Step 2: Write the failing test**

Create `ddbj_mss_tools/tests/test_transl_except_roundtrip.py`:
```python
import os
from Bio import SeqIO
from ddbj_gff import parse
from ddbj_gff.flatfile import flatfile_to_gff
from ddbj_gff.writer import write
from gff2mss.convert import build_entry_features
from gff2mss.config import MssConfig

FIX = os.path.join(os.path.dirname(__file__), "flatfile_fixtures", "transl_except_p87.gbk")


def _cds_feats():
    rec = SeqIO.read(FIX, "genbank")
    doc = parse(write(flatfile_to_gff(rec)))               # flatfile -> canonical GFF (recoded_codon child)
    seqs = {rec.id: rec.seq}
    cfg = MssConfig(source={}, transl_table=1, product_default="hypothetical protein")
    cfg.emit_mrna = False                                   # virus 2-level
    per_entry = build_entry_features(doc, seqs, cfg, [])
    return [f for feats in per_entry.values() for f in feats]


def test_transl_except_qualifier_emitted_on_clean_cds():
    feats = _cds_feats()
    cds = [f for f in feats if f.key == "CDS"]
    assert len(cds) == 1                                    # translated cleanly via Pyl -> CDS, not misc_feature
    te = [q.value for q in cds[0].qualifiers if q.key == "transl_except"]
    assert te == ["(pos:746..748,aa:Pyl)"]                  # qualifier emitted, matches original
```

- [ ] **Step 3: Sync + run to verify failure**

```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/common/. ddbj-gff-dev:/opt/mss_src/common
docker exec ddbj-gff-dev bash -lc 'mkdir -p /opt/mss_tests/flatfile_fixtures'
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/flatfile_fixtures/transl_except_p87.gbk ddbj-gff-dev:/opt/mss_tests/flatfile_fixtures/transl_except_p87.gbk
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/test_transl_except_roundtrip.py ddbj-gff-dev:/opt/mss_tests/test_transl_except_roundtrip.py
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_transl_except_roundtrip.py::test_transl_except_qualifier_emitted_on_clean_cds -v'
```
Expected: FAIL — `te == []` (the CDS is emitted as a clean CDS, but with no `transl_except` qualifier).

- [ ] **Step 4: Implement the emit**

In `ddbj_mss_tools/src/gff2mss/convert.py`, `build_cds_feature`, the clean-translation path builds:
```python
    quals = [
        MssQualifier("locus_tag", locus_tag),
        MssQualifier("transl_table", str(table_id)),
        MssQualifier("codon_start", str(codon_start)),
        MssQualifier("product", _product(mrna, gene, cfg)),
    ]
    if gene.gene or mrna.gene:
        quals.append(MssQualifier("gene", gene.gene or mrna.gene))
    inference = mrna.attributes.get("inference")
    if inference:
        quals.append(MssQualifier("inference", inference[0]))
    quals.append(_submitter_note(gene, mrna))
    return MssFeature("CDS", location, quals)
```
Insert the transl_except emission just before `quals.append(_submitter_note(gene, mrna))` (so it groups with the CDS's own qualifiers, before the submitter note):
```python
    for spec in excepts:
        quals.append(MssQualifier("transl_except", spec))
    quals.append(_submitter_note(gene, mrna))
    return MssFeature("CDS", location, quals)
```
(`excepts` is already computed earlier in the function. Do NOT modify the misc_feature return path.)

In `_build_trans_spliced_cds`, add the out-of-scope-combo diagnostic. Right after `parts = cds_feat.ordered_spans()` near the top of that function, add:
```python
    if _collect_transl_excepts(cds_feat):
        diagnostics.append(Diagnostic(Severity.WARNING, None, "transl-except-trans-splicing",
                                      f"CDS {mrna.id!r} combines transl_except with trans-splicing; "
                                      f"the trans-spliced translation path does not apply transl_except "
                                      f"(qualifier round-trip only, translation may be incorrect)"))
```
(Do not change the rest of `_build_trans_spliced_cds`.)

- [ ] **Step 5: Sync + run to verify pass**

```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_transl_except_roundtrip.py -v'
```
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add src/gff2mss/convert.py tests/test_transl_except_roundtrip.py tests/flatfile_fixtures/transl_except_p87.gbk
git commit -m "feat(gff2mss): emit /transl_except from recoded_codon children (round-trip)"
```

---

### Task 2: end-to-end round-trip (mss2ff) + full regression

**Files:**
- Test: `ddbj_mss_tools/tests/test_transl_except_roundtrip.py` (extend)

**Interfaces:**
- Consumes: `ddbj_gff.flatfile.flatfile_to_gff`, `ddbj_gff.writer.write`, `ddbj_gff.parse`, `gff2mss.convert.convert`, `gff2mss.emit.emit_ann`, `gff2mss.config.MssConfig`, `mss2ff.cli`.
- Produces: proof the `/transl_except` qualifier + the `O`-containing translation survive flatfile → GFF → gff2mss → mss2ff.

- [ ] **Step 1: Write the round-trip test (extend)**

Append to `ddbj_mss_tools/tests/test_transl_except_roundtrip.py`:
```python
def test_transl_except_full_roundtrip(tmp_path):
    from gff2mss.convert import convert
    from gff2mss.emit import emit_ann
    import mss2ff.cli
    rec = SeqIO.read(FIX, "genbank")
    doc = parse(write(flatfile_to_gff(rec)))
    cfg = MssConfig(source={}, transl_table=1, product_default="hypothetical protein")
    cfg.emit_mrna = False
    mss_doc, _ = convert(doc, {rec.id: rec.seq}, cfg, common_rows=[])
    ann = emit_ann(mss_doc)
    assert "\ttransl_except\t(pos:746..748,aa:Pyl)" in ann     # qualifier in the .ann

    ann_p = tmp_path / "o.ann"; ann_p.write_text(ann, encoding="utf-8")
    fa_p = tmp_path / "o.fasta"; fa_p.write_text(f">{rec.id}\n{str(rec.seq)}\n", encoding="utf-8")
    ff_p = tmp_path / "o.ff"
    mss2ff.cli.main([str(ann_p), str(fa_p), "-o", str(ff_p)])
    rt = SeqIO.read(str(ff_p), "genbank")
    cds = [f for f in rt.features if f.type == "CDS"]
    assert len(cds) == 1
    assert cds[0].qualifiers.get("transl_except") == ["(pos:746..748,aa:Pyl)"]   # /transl_except preserved
    assert "O" in cds[0].qualifiers["translation"][0]           # Pyl -> O in the regenerated protein
```
> If the `.ann` substring assertion fails, inspect the emitted qualifier row (`emit.feature_rows` format) and adjust the exact whitespace while keeping the intent (a `transl_except` row with value `(pos:746..748,aa:Pyl)`). If `mss2ff` renders `/transl_except` differently, adjust the `.ff` assertion to the actual key it uses — but do NOT weaken it to a substring-only check.

- [ ] **Step 2: Sync + run**

```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/test_transl_except_roundtrip.py ddbj-gff-dev:/opt/mss_tests/test_transl_except_roundtrip.py
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_transl_except_roundtrip.py -v'
```
Expected: PASS (both tests).

- [ ] **Step 3: Full gff2mss regression**

```bash
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest . -q'
```
Expected: PASS (prior suite + new tests; the pre-existing `@slow` skip). Zero failures. In particular the existing trans-splicing and origin-spanning transl_except tests stay green (the new emit only adds a qualifier on the clean-CDS path).

- [ ] **Step 4: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add tests/test_transl_except_roundtrip.py
git commit -m "test(gff2mss): end-to-end transl_except round-trip (flatfile -> GFF -> gff2mss -> mss2ff)"
```

---

## Self-Review

**Spec coverage:** emit `/transl_except` from `excepts` on the clean-CDS path → Task 1 Step 4. misc_feature path untouched → Task 1 Step 4 (explicit). trans-spliced+transl_except diagnostic (out-of-scope combo) → Task 1 Step 4. Verified spec value `(pos:746..748,aa:Pyl)` → asserted Task 1 Step 2 + Task 2. Full round-trip through mss2ff (qualifier + `O` translation) → Task 2. Regression → Task 2 Step 3.

**Placeholder scan:** No TBD/TODO. Every code step is complete. The Task 2 contingency note points at the concrete emitter/mss2ff format to check, not a vague instruction.

**Type consistency:** `MssQualifier("transl_except", spec)` where `spec` is a str from `_collect_transl_excepts` (list[str]). `build_entry_features(doc, seqs, cfg, diagnostics) -> dict`; `MssFeature.key`/`.qualifiers`; `convert(doc, seqs, cfg, common_rows=[]) -> (MssDocument, diags)`; `emit_ann(MssDocument) -> str`; `mss2ff.cli.main([ann, fasta, "-o", out])` — all match established usage in prior plans. `cfg.emit_mrna=False` matches organelle/virus 2-level.
