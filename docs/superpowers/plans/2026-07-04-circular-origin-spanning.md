# Circular Origin-Spanning Features Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Correctly canonicalize and convert features that cross the origin of a circular molecule (INSDC `end > seqlen` single-row convention) so `ddbj-gff` validates them and `gff2mss` emits the proper `join()`/`complement(join())` MSS location and a correct translation.

**Architecture:** Two repos, one-way dependency (`gff2mss → ddbj-gff`). `ddbj-gff` (canonicalization) gains a shared `GffDocument.circular_seqids` property, a `pass_circular_origin` normalize pass that propagates `Is_circular` onto origin-spanning features (coordinates kept as `end>seqlen`), and a landmark-aware `rule_seqid_bounds`. `gff2mss` (conversion) gains a `_wrap_spans` helper that splits an `end>seqlen` span into two in-bounds pieces in biological 5′→3′ order, wired into `build_insdc_location` and `extract_seq`.

**Tech Stack:** Python 3.10+, BioPython 1.87 (`_insdc_location_string`, `CompoundLocation`), pydantic 2, pytest 9. All tests run inside the `ddbj-gff-dev` Docker container.

## Global Constraints

- **Two repos, one-way dependency:** `gff2mss` imports `ddbj_gff`; `ddbj_gff` MUST NOT import `gff2mss`/`common`. Never introduce a reverse or circular import.
- **Canonical representation:** origin-spanning features keep the INSDC single-row `end > seqlen` coordinate (`end = true_end + landmark_length`); do NOT rewrite coordinates into two spans in the canonical GFF. `Is_circular=true` is propagated onto the feature by the normalize pass.
- **Verified location strings (exact, do not alter):** for span `4447..5268` on a circular molecule of length `5125` →
  - minus: `complement(join(4447..5125,1..143))`
  - plus: `join(4447..5125,1..143)`
  - minus + 5′-partial: `complement(join(4447..5125,1..>143))`
  - minus + 3′-partial: `complement(join(<4447..5125,1..143))`
- **Biological part order (the rule that makes BioPython emit the above):** plus → `[head, tail]`, minus → `[tail, head]`, where `head = start..L`, `tail = 1..(end−L)`. `_ordered()` gives the WRONG order for wraps — bypass it on the wrap path.
- **Commit policy:**
  - `gff_submission` (Tasks 1–3): work on branch `feat/circular-origin-spanning` off `main` (@ `c689cd2`). Normal commits (`git add` the specific files).
  - `ddbj_mss_tools` (Tasks 4–5): work on branch `feat/circular-origin-spanning` off `main` (@ `093830e`). Commit ONLY `src/gff2mss/**` and `tests/**` files you create/modify. NEVER `git add -A` — there are ~49 pre-existing uncommitted files in `examples/`/`docs/`/`data/` that MUST stay untouched. Do NOT push.
- **Reference spec:** `docs/superpowers/specs/2026-07-04-circular-origin-spanning-design.md`.
- **Fixture:** `gff_submission/tests/normalize_fixtures/cp187952_origin.{gff3,fasta,README.md}` (real INSDC-derived; modA CDS `ACPZ3T_00005`, minus, `4447..5268`, L=5125, 274 codons, protein starts `MKL`).

## Test Environment (both repos, `ddbj-gff-dev` container)

The container is already provisioned:
- `ddbj_gff` source is bind-mounted live at `/workspace/src`; tests at `/workspace/tests` (Tasks 1–3 need NO sync).
- `gff2mss` + `common` live at `/opt/mss_src` (docker-cp'd — NOT live; re-sync after every edit). ddbj_mss_tools tests live at `/opt/mss_tests`.
- venv `/opt/ddbj-venv` (Bio 1.87, pydantic 2.13, pytest 9.1).

**Run ddbj-gff tests (Tasks 1–3):**
```bash
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev \
  bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/<file> -v'
```

**Sync + run gff2mss tests (Tasks 4–5):**
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/<file> ddbj-gff-dev:/opt/mss_tests/<file>
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest <file> -v'
```
> `docker cp .../gff2mss/.` (trailing `/.`) copies contents into the existing dir; without it the copy NESTS at `/opt/mss_src/gff2mss/gff2mss`.

## File Structure

- `gff_submission/src/ddbj_gff/model.py` — add `GffDocument.circular_seqids` property (shared by pass + rule).
- `gff_submission/src/ddbj_gff/normalize/passes.py` — add `pass_circular_origin`.
- `gff_submission/src/ddbj_gff/normalize/normalize.py` — register the pass in `ALL_PASSES`.
- `gff_submission/src/ddbj_gff/validate/rules.py` — make `rule_seqid_bounds` landmark-aware.
- `ddbj_mss_tools/src/gff2mss/convert.py` — add `_wrap_spans`; wire into `build_insdc_location`, `extract_seq`; add multi-span+wrap diagnostic in `build_cds_feature`.
- New tests: `gff_submission/tests/test_model_circular_seqids.py`, `.../tests/test_normalize_circular_origin.py`; extend `.../tests/test_validate_rules_header.py`; extend `ddbj_mss_tools/tests/test_mss_location.py`; new `ddbj_mss_tools/tests/test_gff2mss_origin_spanning.py`; copy fixture into `ddbj_mss_tools/tests/mss_fixtures/`.

---

### Task 1: `GffDocument.circular_seqids` (ddbj-gff, model)

**Files:**
- Modify: `gff_submission/src/ddbj_gff/model.py` (add property to `GffDocument`, near `sequence_regions`)
- Test: `gff_submission/tests/test_model_circular_seqids.py`

**Interfaces:**
- Consumes: `GffDocument.features` (list of `Feature`), `Feature.type`, `Feature.is_circular`, `Feature.spans` (list of `Span`, each has `.seqid`).
- Produces: `GffDocument.circular_seqids -> set[str]` — the set of seqids whose landmark feature (`type in {"region","source"}`) has `Is_circular=true`. Consumed by Tasks 2 and 3.

- [ ] **Step 1: Write the failing test**

Create `gff_submission/tests/test_model_circular_seqids.py`:
```python
from ddbj_gff.model import Feature, Span, GffDocument


def test_circular_seqids_from_region_landmark():
    region = Feature("r", "S", "region", [Span("CP", 1, 100, "+")], {"Is_circular": ["true"]}, [])
    gene = Feature("g", "S", "gene", [Span("CP", 90, 130, "+")], {}, [])
    linear = Feature("l", "S", "region", [Span("MT", 1, 50, "+")], {}, [])
    doc = GffDocument(features=[region, gene, linear])
    assert doc.circular_seqids == {"CP"}


def test_circular_seqids_empty_when_no_landmark_flag():
    gene = Feature("g", "S", "gene", [Span("CP", 1, 9, "+")], {"Is_circular": ["true"]}, [])
    # flag on a non-landmark feature type does not make the seqid circular
    doc = GffDocument(features=[gene])
    assert doc.circular_seqids == set()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_model_circular_seqids.py -v'`
Expected: FAIL — `AttributeError: 'GffDocument' object has no attribute 'circular_seqids'`.

- [ ] **Step 3: Add the property**

In `gff_submission/src/ddbj_gff/model.py`, in the `GffDocument` class immediately after the `sequence_regions` property, add:
```python
    @property
    def circular_seqids(self) -> set[str]:
        """Seqids whose landmark feature (region/source) is marked Is_circular=true."""
        out: set[str] = set()
        for f in self.features:
            if f.type in ("region", "source") and f.is_circular:
                for s in f.spans:
                    out.add(s.seqid)
        return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_model_circular_seqids.py -v'`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/gff_submission
git add src/ddbj_gff/model.py tests/test_model_circular_seqids.py
git commit -m "feat(model): GffDocument.circular_seqids from landmark Is_circular"
```

---

### Task 2: `pass_circular_origin` normalize pass (ddbj-gff)

**Files:**
- Modify: `gff_submission/src/ddbj_gff/normalize/passes.py` (add `pass_circular_origin`)
- Modify: `gff_submission/src/ddbj_gff/normalize/normalize.py` (import + register in `ALL_PASSES`)
- Test: `gff_submission/tests/test_normalize_circular_origin.py`

**Interfaces:**
- Consumes: `doc.circular_seqids` (Task 1), `doc.sequence_regions -> dict[seqid, (start, end)]`, `ctx.seq_lengths` (`dict|None`), `Feature.attributes` (`dict[str, list[str]]`), `Change(action, target, message)` from `.report`.
- Produces: `pass_circular_origin(doc, ctx) -> list[Change]`. Sets `f.attributes["Is_circular"] = ["true"]` on every feature with a span `end > seqlen` on a circular seqid. Registered in `ALL_PASSES` after `pass_wrap_cds_in_mrna`.

- [ ] **Step 1: Write the failing test**

Create `gff_submission/tests/test_normalize_circular_origin.py`:
```python
from ddbj_gff import parse
from ddbj_gff.model import Feature, Span, Directive, GffDocument
from ddbj_gff.normalize.passes import pass_circular_origin, NormalizeContext
from ddbj_gff.normalize.normalize import normalize


def _ctx():
    return NormalizeContext(vocab=None, seq_lengths=None, config=None)


def test_pass_flags_origin_spanning_feature_on_circular_seqid():
    region = Feature("r", "S", "region", [Span("CP", 1, 100, "+")], {"Is_circular": ["true"]}, [])
    wrap = Feature("w", "S", "gene", [Span("CP", 90, 130, "+")], {}, [])   # end 130 > 100
    inside = Feature("i", "S", "gene", [Span("CP", 10, 40, "+")], {}, [])  # within bounds
    seq_dir = Directive("x", "sequence-region", ("CP", 1, 100))
    doc = GffDocument(directives=[seq_dir], features=[region, wrap, inside])
    changes = pass_circular_origin(doc, _ctx())
    assert wrap.attributes.get("Is_circular") == ["true"]
    assert "Is_circular" not in inside.attributes
    assert len(changes) == 1


def test_pass_noop_when_seqid_not_circular():
    wrap = Feature("w", "S", "gene", [Span("MT", 90, 130, "+")], {}, [])
    seq_dir = Directive("x", "sequence-region", ("MT", 1, 100))
    doc = GffDocument(directives=[seq_dir], features=[wrap])
    assert pass_circular_origin(doc, _ctx()) == []
    assert "Is_circular" not in wrap.attributes


def test_full_normalize_flags_moda_on_cp187952():
    with open("tests/normalize_fixtures/cp187952_origin.gff3") as fh:
        doc = parse(fh.read())
    work, _report = normalize(doc)
    moda = [f for f in work.features
            if f.type in ("gene", "CDS") and f._first("locus_tag") == "ACPZ3T_00005"]
    assert moda, "modA gene/CDS not found after normalize"
    assert all(f.is_circular for f in moda)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_normalize_circular_origin.py -v'`
Expected: FAIL — `ImportError: cannot import name 'pass_circular_origin'`.

- [ ] **Step 3: Add the pass**

In `gff_submission/src/ddbj_gff/normalize/passes.py`, append:
```python
def pass_circular_origin(doc, ctx) -> list:
    """Propagate Is_circular=true onto origin-spanning features (a span with
    end>seqlen) on a circular landmark. Coordinates are left as-is (canonical keeps
    the INSDC end>seqlen convention); the flag lets validate/downstream treat the
    feature as circular."""
    changes: list = []
    circular = doc.circular_seqids
    if not circular:
        return changes
    regions = doc.sequence_regions
    for f in doc.features:
        for s in f.spans:
            if s.seqid not in circular:
                continue
            seqlen = regions.get(s.seqid, (None, None))[1]
            if seqlen is None and ctx.seq_lengths:
                seqlen = ctx.seq_lengths.get(s.seqid)
            if seqlen is not None and s.end > seqlen:
                if f.attributes.get("Is_circular") != ["true"]:
                    f.attributes["Is_circular"] = ["true"]
                    changes.append(Change("add-qualifier", f.id or "?",
                                          f"propagated Is_circular=true to origin-spanning "
                                          f"feature (span {s.start}..{s.end} > seqlen {seqlen})"))
                break
    return changes
```

- [ ] **Step 4: Register the pass**

In `gff_submission/src/ddbj_gff/normalize/normalize.py`, add `pass_circular_origin` to the import from `.passes` and insert it into `ALL_PASSES` immediately after `pass_wrap_cds_in_mrna`:
```python
from .passes import (NormalizeContext, pass_directives, pass_coerce_transcript_to_mrna,
                     pass_wrap_cds_in_mrna, pass_circular_origin, pass_so_terms,
                     pass_transl_except, pass_anticodon)

ALL_PASSES = [pass_directives, pass_coerce_transcript_to_mrna, pass_wrap_cds_in_mrna,
              pass_circular_origin, pass_so_terms, pass_transl_except, pass_anticodon]
```

- [ ] **Step 5: Run test to verify it passes**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_normalize_circular_origin.py -v'`
Expected: PASS (3 passed).

- [ ] **Step 6: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/gff_submission
git add src/ddbj_gff/normalize/passes.py src/ddbj_gff/normalize/normalize.py tests/test_normalize_circular_origin.py
git commit -m "feat(normalize): pass_circular_origin propagates Is_circular to end>seqlen features"
```

---

### Task 3: landmark-aware `rule_seqid_bounds` (ddbj-gff, validate)

**Files:**
- Modify: `gff_submission/src/ddbj_gff/validate/rules.py` (`rule_seqid_bounds`, currently ~lines 48–65)
- Test: `gff_submission/tests/test_validate_rules_header.py` (extend; imports `Feature, Span, Directive, GffDocument`, `rules`, helper `codes`, `V` already present)

**Interfaces:**
- Consumes: `doc.circular_seqids` (Task 1), `doc.sequence_regions`, `Feature.is_circular`, `make_diagnostic(code, message)`.
- Produces: `rule_seqid_bounds` no longer emits `feature-outside-region` for `end>hi` when the seqid is circular (landmark flag OR feature flag); still emits it for non-circular `end>hi` and for any `start<lo`.

- [ ] **Step 1: Write the failing tests**

Append to `gff_submission/tests/test_validate_rules_header.py`:
```python
def test_seqid_bounds_circular_landmark_allows_origin_spanning():
    region = Feature("r", "S", "region", [Span("c", 1, 100, "+")], {"Is_circular": ["true"]}, [])
    seq_dir = Directive("x", "sequence-region", ("c", 1, 100))
    cds = Feature("a", "S", "CDS", [Span("c", 90, 130, "+")], {}, [])  # flag on region, not cds
    doc = GffDocument(directives=[seq_dir], features=[region, cds])
    assert "feature-outside-region" not in codes(rules.rule_seqid_bounds(doc, V))


def test_seqid_bounds_noncircular_end_beyond_region_flagged():
    seq_dir = Directive("x", "sequence-region", ("c", 1, 100))
    cds = Feature("a", "S", "CDS", [Span("c", 90, 130, "+")], {}, [])
    doc = GffDocument(directives=[seq_dir], features=[cds])
    assert "feature-outside-region" in codes(rules.rule_seqid_bounds(doc, V))
```

- [ ] **Step 2: Run tests to verify the first fails**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_validate_rules_header.py -k "circular_landmark or noncircular_end_beyond" -v'`
Expected: `test_seqid_bounds_circular_landmark_allows_origin_spanning` FAILS (rule flags it because the flag is on the region, not the CDS); `test_seqid_bounds_noncircular_end_beyond_region_flagged` PASSES.

- [ ] **Step 3: Make the rule landmark-aware**

Replace `rule_seqid_bounds` in `gff_submission/src/ddbj_gff/validate/rules.py` with:
```python
def rule_seqid_bounds(doc, vocab) -> list:
    diags = []
    regions = doc.sequence_regions
    circular_seqids = doc.circular_seqids
    for f in doc.features:
        for s in f.spans:
            if s.seqid not in regions:
                diags.append(make_diagnostic("undefined-seqid",
                                             f"feature {f.id!r} references seqid {s.seqid!r} "
                                             f"with no ##sequence-region"))
                continue
            lo, hi = regions[s.seqid]
            circular = f.is_circular or s.seqid in circular_seqids
            if s.start < lo:
                diags.append(make_diagnostic("feature-outside-region",
                                             f"feature {f.id!r} span {s.start}..{s.end} is outside "
                                             f"sequence-region {s.seqid}:{lo}..{hi}"))
            elif s.end > hi and not circular:
                diags.append(make_diagnostic("feature-outside-region",
                                             f"feature {f.id!r} span {s.start}..{s.end} is outside "
                                             f"sequence-region {s.seqid}:{lo}..{hi}"))
    return diags
```

- [ ] **Step 4: Run tests to verify they pass (incl. the existing circular test)**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_validate_rules_header.py -v'`
Expected: PASS — including the pre-existing `test_seqid_bounds_circular_origin_spanning_allowed` (feature-level flag) and `test_start_gt_end`.

- [ ] **Step 5: Full ddbj-gff regression (not-slow)**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests -m "not slow" -q'`
Expected: PASS (149 passed — the prior 147 plus the 2 new files' tests). Zero failures.

- [ ] **Step 6: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/gff_submission
git add src/ddbj_gff/validate/rules.py tests/test_validate_rules_header.py
git commit -m "feat(validate): rule_seqid_bounds treats circular landmarks as origin-spanning-ok"
```

---

### Task 4: origin-spanning location + translation in gff2mss

**Files:**
- Modify: `ddbj_mss_tools/src/gff2mss/convert.py` (add `Span` import; add `_wrap_spans`; wire into `build_insdc_location` and `extract_seq`; add multi-span+wrap diagnostic in `build_cds_feature`)
- Test: `ddbj_mss_tools/tests/test_mss_location.py` (extend)

**Interfaces:**
- Consumes: `ddbj_gff.model.Span`, existing `_STRAND`, `_ordered`, `FeatureLocation`, `CompoundLocation`, `_insdc_location_string`, `Diagnostic`, `Severity`.
- Produces: `_wrap_spans(spans, seqlen) -> (list[Span], bool)`; unchanged signatures for `build_insdc_location(spans, seqlen, five_prime_partial=False, three_prime_partial=False) -> str` and `extract_seq(spans, genome_seq) -> Seq`, now wrap-aware.

- [ ] **Step 1: Write the failing tests**

Append to `ddbj_mss_tools/tests/test_mss_location.py`:
```python
def test_minus_strand_origin_spanning():
    assert build_insdc_location([Span("c", 4447, 5268, "-")], 5125) == "complement(join(4447..5125,1..143))"


def test_plus_strand_origin_spanning():
    assert build_insdc_location([Span("c", 4447, 5268, "+")], 5125) == "join(4447..5125,1..143)"


def test_origin_spanning_partials_minus():
    assert build_insdc_location([Span("c", 4447, 5268, "-")], 5125,
                                five_prime_partial=True) == "complement(join(4447..5125,1..>143))"
    assert build_insdc_location([Span("c", 4447, 5268, "-")], 5125,
                                three_prime_partial=True) == "complement(join(<4447..5125,1..143))"


def test_extract_origin_spanning_plus_translates():
    # 9 bp circular; plus CDS 7..12 wraps: head=7..9 "ATG", tail=1..3 "TAA" -> "ATGTAA" -> M*
    genome = Seq("TAACCCATG")
    ex = extract_seq([Span("c", 7, 12, "+")], genome)
    assert str(ex) == "ATGTAA"
    assert str(ex.translate(table=11)) == "M*"


def test_extract_origin_spanning_minus_translates():
    # minus CDS 7..12 on revcomp genome yields the same coding sequence
    genome = Seq("TAACCCATG").reverse_complement()  # so minus strand of 7..12 -> ATGTAA
    ex = extract_seq([Span("c", 7, 12, "-")], genome)
    assert str(ex.translate(table=11)) == "M*"
```

- [ ] **Step 2: Sync + run to verify failure**

Run:
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/test_mss_location.py ddbj-gff-dev:/opt/mss_tests/test_mss_location.py
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_mss_location.py -v'
```
Expected: the 5 new tests FAIL (`build_insdc_location` returns `complement(4447..5268)` / `extract_seq` truncates past the sequence end); existing tests still pass.

- [ ] **Step 3: Add `Span` import + `_wrap_spans`**

In `ddbj_mss_tools/src/gff2mss/convert.py`, add to the imports (with the other `ddbj_gff` imports):
```python
from ddbj_gff.model import Span
```
Add `_wrap_spans` immediately after the `_ordered` helper:
```python
def _wrap_spans(spans, seqlen):
    """Split an origin-spanning span (end>seqlen) into its two in-bounds pieces in
    biological 5'->3' order: plus [head, tail], minus [tail, head], where head=start..L
    and tail=1..(end-L). Non-wrapping spans are returned unchanged.
    Returns (spans, wrapped: bool)."""
    out, wrapped = [], False
    for s in spans:
        if s.end > seqlen:
            wrapped = True
            head = Span(s.seqid, s.start, seqlen, s.strand)
            tail = Span(s.seqid, 1, s.end - seqlen, s.strand)
            out += [tail, head] if s.strand == "-" else [head, tail]
        else:
            out.append(s)
    return out, wrapped
```

- [ ] **Step 4: Wire `_wrap_spans` into `build_insdc_location` and `extract_seq`**

In `build_insdc_location`, replace the `ordered = _ordered(spans)` line with:
```python
    wrapped_spans, wrapped = _wrap_spans(spans, seqlen)
    ordered = wrapped_spans if wrapped else _ordered(spans)
```
(The rest — `locs = [...]`, the partial branches, `CompoundLocation`, `_insdc_location_string` — is unchanged; on the wrap path `locs[0]` is the biological 5′ end and `locs[-1]` the 3′ end, so the existing partial logic applies as-is.)

In `extract_seq`, replace the `ordered = _ordered(spans)` line with:
```python
    wrapped_spans, wrapped = _wrap_spans(spans, len(genome_seq))
    ordered = wrapped_spans if wrapped else _ordered(spans)
```

- [ ] **Step 5: Sync + run to verify pass**

Run:
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_mss_location.py -v'
```
Expected: PASS (all, including the 5 new tests).

- [ ] **Step 6: Add multi-span + wrap diagnostic (non-goal safety)**

In `build_cds_feature` (`convert.py`), immediately after the `if not spans:` early-return block, add:
```python
    if len(spans) > 1 and any(s.end > len(genome_seq) for s in spans):
        diagnostics.append(Diagnostic(Severity.WARNING, None, "multi-exon-origin-spanning",
                                      f"CDS {mrna.id!r} is multi-exon and origin-spanning; "
                                      f"join order is best-effort (unsupported combination)"))
```
Add a test to `ddbj_mss_tools/tests/test_mss_location.py` (top of file needs `from gff2mss.convert import build_cds_feature` — add to the existing import line if not present; and `from ddbj_gff.model import Feature`):
```python
def test_multi_exon_wrap_warns():
    from ddbj_gff.model import Feature
    from gff2mss.convert import build_cds_feature
    cds = Feature("cds", "S", "CDS",
                  [Span("c", 1, 6, "+"), Span("c", 7, 14, "+")], {"transl_table": ["11"]}, [])  # 14 > 10
    mrna = Feature("m", "S", "mRNA", [Span("c", 1, 14, "+")], {}, [])
    mrna.children = [cds]
    diags = []
    build_cds_feature(mrna, Feature("g", "S", "gene", [Span("c", 1, 14, "+")], {}, []),
                      "T_0001", Seq("ATGAAATAACC"[:10] + "A"), _MinimalCfg(), diags)
    assert any(d.code == "multi-exon-origin-spanning" for d in diags)
```
> If `build_cds_feature`'s config/argument shape makes this call awkward, simplify the test to construct only what the diagnostic branch needs and assert the warning code; the diagnostic branch runs before translation. Confirm the exact `build_cds_feature` signature in `convert.py` and match it (do not invent parameters). If a minimal call is impractical, assert the branch via a direct unit check on `_wrap_spans` returning `wrapped=True` for the multi-span input plus a code-review note — but prefer exercising `build_cds_feature`.

- [ ] **Step 7: Sync + run to verify pass**

Run:
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/test_mss_location.py ddbj-gff-dev:/opt/mss_tests/test_mss_location.py
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_mss_location.py -v'
```
Expected: PASS.

- [ ] **Step 8: Commit (gff2mss files + test only)**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add src/gff2mss/convert.py tests/test_mss_location.py
git commit -m "feat(gff2mss): origin-spanning join()/complement + translation via _wrap_spans"
```

---

### Task 5: end-to-end integration on cp187952 + full regression

**Files:**
- Create: `ddbj_mss_tools/tests/mss_fixtures/cp187952_origin.gff3` (copy of the ddbj-gff fixture)
- Create: `ddbj_mss_tools/tests/mss_fixtures/cp187952_origin.fasta` (copy)
- Create: `ddbj_mss_tools/tests/test_gff2mss_origin_spanning.py`

**Interfaces:**
- Consumes: `ddbj_gff.parse`, `ddbj_gff.normalize.normalize.normalize`, `ddbj_gff.validate.validate`, `gff2mss.convert.build_entry_features`, `gff2mss.config.load_config`, `Bio.SeqIO`. `build_entry_features(doc, seqs, cfg, diagnostics) -> dict[seqid, list[MssFeature]]`; `MssFeature.key` (type) and `MssFeature.location`.
- Produces: proof that the full chain (parse → normalize → validate → gff2mss) yields modA as a `CDS` at `complement(join(4447..5125,1..143))` (a wrong translation would instead yield a `misc_feature`).

- [ ] **Step 1: Copy the fixture into ddbj_mss_tools**

```bash
cp /Users/tanizawa/projects/ddbj/gff_submission/tests/normalize_fixtures/cp187952_origin.gff3 \
   /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/mss_fixtures/cp187952_origin.gff3
cp /Users/tanizawa/projects/ddbj/gff_submission/tests/normalize_fixtures/cp187952_origin.fasta \
   /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/mss_fixtures/cp187952_origin.fasta
```

- [ ] **Step 2: Write the failing test**

Create `ddbj_mss_tools/tests/test_gff2mss_origin_spanning.py`:
```python
import os
from Bio import SeqIO
from ddbj_gff import parse
from ddbj_gff.normalize.normalize import normalize
from ddbj_gff.validate import validate
from gff2mss.convert import build_entry_features
from gff2mss.config import load_config

FIX = os.path.join(os.path.dirname(__file__), "mss_fixtures")


def _cfg(tmp_path):
    p = tmp_path / "cp.toml"
    p.write_text('[source]\norganism="Aliinostoc maniaoense"\nmol_type="genomic DNA"\n'
                 '[locus_tag]\nprefix="ACPZ3T"\n[cds]\ntransl_table=11\n'
                 '[transcript]\nemit_mrna=false\n', encoding="utf-8")
    cfg, _ = load_config(str(p))
    return cfg


def test_cp187952_origin_spanning_end_to_end(tmp_path):
    with open(os.path.join(FIX, "cp187952_origin.gff3")) as fh:
        doc = parse(fh.read())
    work, _ = normalize(doc)

    # canonicalization: origin-spanning modA no longer flagged out-of-region
    diags = validate(work)
    assert not any(d.code == "feature-outside-region" for d in diags)

    seqs = {rec.id: rec.seq for rec in
            SeqIO.parse(os.path.join(FIX, "cp187952_origin.fasta"), "fasta")}
    per_entry = build_entry_features(work, seqs, _cfg(tmp_path), [])
    feats = per_entry["CP187952.1"]

    # modA is emitted as a CDS (clean translation) at the wrapped location — not a misc_feature
    cds_locs = [f.location for f in feats if f.key == "CDS"]
    assert "complement(join(4447..5125,1..143))" in cds_locs
    assert not any(f.key == "misc_feature" for f in feats)
```

- [ ] **Step 3: Sync + run to verify failure/pass**

Run:
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/mss_fixtures/cp187952_origin.gff3 ddbj-gff-dev:/opt/mss_tests/mss_fixtures/cp187952_origin.gff3
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/mss_fixtures/cp187952_origin.fasta ddbj-gff-dev:/opt/mss_tests/mss_fixtures/cp187952_origin.fasta
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/test_gff2mss_origin_spanning.py ddbj-gff-dev:/opt/mss_tests/test_gff2mss_origin_spanning.py
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_gff2mss_origin_spanning.py -v'
```
Expected: PASS once Tasks 1–4 are in place (all four are prerequisites). If it fails, read the assertion: a `feature-outside-region` diag means Task 3 is incomplete; a `misc_feature` for modA means `extract_seq` (Task 4) is wrong; a missing CDS location means normalize/`build_entry_features` wiring is off.
> If `load_config` rejects an unknown `[transcript]`/`emit_mrna` key, drop that block (mRNA emission does not affect the CDS assertion); confirm the accepted config keys in `gff2mss/config.py`.

- [ ] **Step 4: Full gff2mss regression**

Run:
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest . -q'
```
Expected: PASS (prior 80 + the new location/e2e tests; 1 pre-existing skip). Zero failures.

- [ ] **Step 5: Commit (gff2mss test + fixtures only)**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add tests/test_gff2mss_origin_spanning.py tests/mss_fixtures/cp187952_origin.gff3 tests/mss_fixtures/cp187952_origin.fasta
git commit -m "test(gff2mss): end-to-end origin-spanning on cp187952 (normalize->convert)"
```

---

## Self-Review

**Spec coverage:**
- A1 `pass_circular_origin` → Task 2. A2 landmark-aware validate → Task 3. Shared circular-seqid detection → Task 1 (`circular_seqids`, DRY across A1/A2). B1 `_wrap_spans` → Task 4 Step 3. B2 `build_insdc_location` → Task 4 Step 4 (+ pinned strings in Step 1). B3 `extract_seq` → Task 4 Step 4 (+ translation tests Step 1). Multi-span+wrap non-goal diagnostic → Task 4 Step 6. End-to-end cp187952 verification → Task 5. Regression (ddbj-gff 147+, gff2mss 80) → Task 3 Step 5 and Task 5 Step 4. Fixture README correction → already applied during design (spec review). Heterosigma byte-identity: no origin-spanning present, so unaffected — the gff2mss regression suite (Task 5 Step 4) covers it.
- Optional ddbj-validator confirmation on the produced `.ann` is spec-optional; not a task (the e2e API-level assertion is the gate).

**Placeholder scan:** No TBD/TODO. Every code step shows complete code. The two `>`-noted fallbacks (Task 4 Step 6, Task 5 Step 3) are contingency instructions tied to concrete signatures the implementer must confirm, not placeholders for logic.

**Type consistency:** `circular_seqids -> set[str]` produced in Task 1, consumed identically in Tasks 2–3. `_wrap_spans(spans, seqlen) -> (list, bool)` defined and used consistently in Task 4. `Span(seqid, start, end, strand)` positional form matches the fixtures' usage in existing tests. `MssFeature.key`/`.location` used in Task 5 match `model.py`. `build_entry_features(doc, seqs, cfg, diagnostics)` signature matches `convert.py`.
