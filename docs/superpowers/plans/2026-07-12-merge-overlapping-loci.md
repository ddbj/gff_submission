# pass_merge_overlapping_loci — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in `ddbj-gff` normalize pass that merges same-strand gene loci whose mRNAs overlap (by a tunable percentage) into a single gene, excluding trans-spliced transcripts.

**Architecture:** All changes in `gff_submission` (`ddbj-gff`): two `NormalizeConfig` fields (flag + threshold) and a new `pass_merge_overlapping_loci` registered in `ALL_PASSES` after `pass_wrap_cds_in_mrna`. Default off → no behaviour change for existing pipelines (heterosigma etc.). The N. benthamiana pipeline enables it (its current input has zero overlaps, so it is a no-op there; this is a general/future-robustness rule).

**Tech Stack:** Python 3.10+, pytest 9. All tests run inside the `ddbj-gff-dev` Docker container.

## Global Constraints

- **Repo:** `gff_submission` only, branch `feat/merge-overlapping-loci` off `main`. Normal commits. No `ddbj_mss_tools` change.
- **Opt-in:** the pass is a no-op unless `NormalizeConfig.merge_overlapping_loci` is `True` (default `False`). Existing normalize behaviour must be **unchanged** when the flag is off.
- **Merge criterion:** same seqid **and** same strand; edge between two mRNAs when
  `overlap_bp / min(len_a, len_b) >= NormalizeConfig.merge_overlap_min_fraction` (float, default `0.0` = any ≥1 bp overlap). Connected components (transitive) form a locus.
- **Trans-spliced transcripts are EXEMPT:** an mRNA is skipped (never merged, never merges others) when `mrna.is_trans_spliced` OR any CDS child `is_trans_spliced`.
- **Representative gene** = the gene of the member mRNA with smallest `(lo, hi, mrna_id)`; its ID/attributes survive; other members' mRNAs are reparented to it; its span becomes the union extent; genes left with no children are removed (from `doc.features` and `doc.feature_index`).
- **Reference spec:** `docs/superpowers/specs/2026-07-12-merge-overlapping-loci-design.md`.

## Test Environment (`ddbj-gff-dev` container)

`ddbj_gff` bind-mounted live at `/workspace/src`, tests at `/workspace/tests` (no sync needed). Run:
```bash
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev \
  bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/<file> -v'
```

## File Structure

- `src/ddbj_gff/normalize/config.py` — add `merge_overlapping_loci` + `merge_overlap_min_fraction` (dataclass + loader).
- `src/ddbj_gff/normalize/passes.py` — add `pass_merge_overlapping_loci`.
- `src/ddbj_gff/normalize/normalize.py` — import + register in `ALL_PASSES`; add `"merge-loci"` to `_APPLIED`.
- Tests: `tests/test_normalize_merge_loci.py` (new).

---

### Task 1: NormalizeConfig fields (flag + threshold)

**Files:**
- Modify: `src/ddbj_gff/normalize/config.py`
- Test: `tests/test_normalize_merge_loci.py` (config part)

**Interfaces:**
- Produces: `NormalizeConfig.merge_overlapping_loci: bool = False`, `NormalizeConfig.merge_overlap_min_fraction: float = 0.0`; `load_normalize_config` reads them from the `[normalize]` table.

- [ ] **Step 1: Write the failing test**

Create `tests/test_normalize_merge_loci.py`:
```python
from ddbj_gff.normalize.config import NormalizeConfig


def test_merge_config_defaults_off():
    c = NormalizeConfig()
    assert c.merge_overlapping_loci is False
    assert c.merge_overlap_min_fraction == 0.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_normalize_merge_loci.py::test_merge_config_defaults_off -v'`
Expected: FAIL — `AttributeError: 'NormalizeConfig' object has no attribute 'merge_overlapping_loci'`.

- [ ] **Step 3: Add the fields**

In `src/ddbj_gff/normalize/config.py`, add the two fields to the dataclass (after `wrap_cds_in_mrna`):
```python
    wrap_cds_in_mrna: bool = True
    merge_overlapping_loci: bool = False
    merge_overlap_min_fraction: float = 0.0
```
and to `load_normalize_config`'s `NormalizeConfig(...)` call (after `wrap_cds_in_mrna=...`):
```python
        wrap_cds_in_mrna=n.get("wrap_cds_in_mrna", True),
        merge_overlapping_loci=n.get("merge_overlapping_loci", False),
        merge_overlap_min_fraction=n.get("merge_overlap_min_fraction", 0.0),
```

- [ ] **Step 4: Run test to verify it passes**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_normalize_merge_loci.py -v'`
Expected: PASS (1 passed).

- [ ] **Step 5: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/gff_submission
git add src/ddbj_gff/normalize/config.py tests/test_normalize_merge_loci.py
git commit -m "feat(normalize): NormalizeConfig merge_overlapping_loci + merge_overlap_min_fraction (default off)"
```

---

### Task 2: pass_merge_overlapping_loci + registration

**Files:**
- Modify: `src/ddbj_gff/normalize/passes.py` (add the pass)
- Modify: `src/ddbj_gff/normalize/normalize.py` (import + `ALL_PASSES` + `_APPLIED`)
- Test: `tests/test_normalize_merge_loci.py` (extend)

**Interfaces:**
- Consumes: `Feature`, `Span`, `Change`, `Feature.is_trans_spliced`, `doc.feature_index`, `doc.features`, mRNA `parent_ids`/`children`/`parents`/`attributes`.
- Produces: `pass_merge_overlapping_loci(doc, ctx) -> list[Change]`. Registered in `ALL_PASSES` after `pass_wrap_cds_in_mrna`; `"merge-loci"` added to `_APPLIED`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_normalize_merge_loci.py`:
```python
from ddbj_gff import parse
from ddbj_gff.normalize.normalize import normalize
from ddbj_gff.normalize.config import NormalizeConfig

HDR = "##gff-version 3\n##sequence-region c 1 100000\n"


def _gff(rows):
    return HDR + "".join(rows)


def _gene(gid, s, e, strand="+"):
    m = f"{gid}.m"
    return (f"c\tx\tgene\t{s}\t{e}\t.\t{strand}\t.\tID={gid}\n"
            f"c\tx\tmRNA\t{s}\t{e}\t.\t{strand}\t.\tID={m};Parent={gid}\n"
            f"c\tx\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID={gid}.c;Parent={m}\n")


def _norm(gff, **cfg):
    doc = parse(gff)
    work, _ = normalize(doc, config=NormalizeConfig(taxid=1, **cfg))
    return work


def _genes(doc):
    return [f for f in doc.features if f.type == "gene"]


def _mrnas(doc):
    return [f for f in doc.features if f.type == "mRNA"]


def test_merge_two_overlapping_same_strand():
    gff = _gff([_gene("gA", 100, 500), _gene("gB", 300, 700)])   # overlap 300..500
    doc = _norm(gff, merge_overlapping_loci=True)
    genes = _genes(doc)
    assert len(genes) == 1 and genes[0].id == "gA"               # rep = lowest start
    assert genes[0].spans[0].start == 100 and genes[0].spans[0].end == 700   # union
    mrnas = _mrnas(doc)
    assert len(mrnas) == 2 and all(m.parent_ids == ["gA"] for m in mrnas)     # both under gA
    assert "gB" not in doc.feature_index                          # merged-away gene removed
    codes = {d.code for d in validate_ok(doc)}
    assert "dangling-parent" not in codes


def validate_ok(doc):
    from ddbj_gff.validate import validate
    return validate(doc)


def test_flag_off_no_change():
    gff = _gff([_gene("gA", 100, 500), _gene("gB", 300, 700)])
    doc = _norm(gff)                                             # flag default off
    assert len(_genes(doc)) == 2 and len(_mrnas(doc)) == 2


def test_opposite_strand_not_merged():
    gff = _gff([_gene("gA", 100, 500, "+"), _gene("gB", 300, 700, "-")])
    doc = _norm(gff, merge_overlapping_loci=True)
    assert len(_genes(doc)) == 2                                 # antisense stays separate


def test_transitive_chain_merged():
    gff = _gff([_gene("gA", 100, 500), _gene("gB", 400, 800), _gene("gC", 700, 1000)])
    doc = _norm(gff, merge_overlapping_loci=True)                # A~B, B~C, A!~C
    genes = _genes(doc)
    assert len(genes) == 1 and genes[0].id == "gA"
    assert genes[0].spans[0].end == 1000
    assert len(_mrnas(doc)) == 3


def test_threshold_below_not_merged():
    # gA 100..500 (len401), gB 300..700 (len401); overlap 201; 201/401 = 0.50
    gff = _gff([_gene("gA", 100, 500), _gene("gB", 300, 700)])
    doc = _norm(gff, merge_overlapping_loci=True, merge_overlap_min_fraction=0.9)
    assert len(_genes(doc)) == 2                                 # 0.50 < 0.90 -> not merged


def test_trans_spliced_exempt():
    # a trans-spliced CDS (two parts) whose mRNA extent 100..900 overlaps a normal gene gN 200..600
    trans = ("c\tx\tgene\t100\t900\t.\t+\t.\tID=gT\n"
             "c\tx\tmRNA\t100\t900\t.\t+\t.\tID=gT.m;Parent=gT\n"
             "c\tx\tCDS\t100\t200\t.\t+\t0\tID=gT.c;Parent=gT.m;exception=trans-splicing;part=1\n"
             "c\tx\tCDS\t800\t900\t.\t+\t0\tID=gT.c;Parent=gT.m;exception=trans-splicing;part=2\n")
    gff = _gff([trans, _gene("gN", 200, 600)])
    doc = _norm(gff, merge_overlapping_loci=True)
    ids = {g.id for g in _genes(doc)}
    assert ids == {"gT", "gN"}                                  # trans-spliced gT exempt; gN untouched
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_normalize_merge_loci.py -v'`
Expected: the merge tests FAIL (no merge happens yet — the pass doesn't exist, so `test_merge_two_overlapping_same_strand` sees 2 genes). `test_flag_off_no_change` and `test_opposite_strand_not_merged` may already pass (no merge).

- [ ] **Step 3: Add the pass**

In `src/ddbj_gff/normalize/passes.py`, append:
```python
def _find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def pass_merge_overlapping_loci(doc, ctx) -> list:
    """Merge same-strand gene loci whose mRNAs overlap into one gene (opt-in).

    Edge between two mRNAs when overlap_bp/min(len) >= config.merge_overlap_min_fraction
    (default 0.0 = any overlap). Connected components form a locus; the gene of the
    lowest-coordinate mRNA is the representative, the others' mRNAs are reparented to it,
    and its span becomes the union. Trans-spliced transcripts (mRNA or a CDS child with
    exception=trans-splicing) are excluded entirely."""
    changes: list = []
    if not getattr(ctx.config, "merge_overlapping_loci", False):
        return changes
    frac = getattr(ctx.config, "merge_overlap_min_fraction", 0.0)

    def _is_trans(m):
        return m.is_trans_spliced or any(
            c.type == "CDS" and c.is_trans_spliced for c in m.children)

    def _gene_of(m):
        pid = m.parent_ids[0] if m.parent_ids else None
        return doc.feature_index.get(pid) if pid else None

    groups: dict = {}
    for f in doc.features:
        if f.type != "mRNA" or not f.spans or _is_trans(f):
            continue
        lo = min(s.start for s in f.spans)
        hi = max(s.end for s in f.spans)
        groups.setdefault((f.spans[0].seqid, f.spans[0].strand), []).append((lo, hi, f))

    touched: set = set()
    for (seqid, strand), items in groups.items():
        items.sort(key=lambda t: (t[0], t[1], t[2].id or ""))
        n = len(items)
        parent = list(range(n))
        for i in range(n):
            lo_i, hi_i, _ = items[i]
            for j in range(i + 1, n):
                lo_j, hi_j, _ = items[j]
                if lo_j > hi_i:
                    break
                ov = min(hi_i, hi_j) - max(lo_i, lo_j) + 1
                if ov > 0 and ov / min(hi_i - lo_i + 1, hi_j - lo_j + 1) >= frac:
                    ri, rj = _find(parent, i), _find(parent, j)
                    if ri != rj:
                        parent[rj] = ri
        comps: dict = {}
        for i in range(n):
            comps.setdefault(_find(parent, i), []).append(items[i][2])
        for members in comps.values():
            genes = {}
            for m in members:
                g = _gene_of(m)
                if g is not None:
                    genes[g.id] = g
            if len(genes) < 2:
                continue
            members.sort(key=lambda m: (min(s.start for s in m.spans),
                                        max(s.end for s in m.spans), m.id or ""))
            rep = _gene_of(members[0])
            u_lo = min(min(s.start for s in m.spans) for m in members)
            u_hi = max(max(s.end for s in m.spans) for m in members)
            for m in members:
                g = _gene_of(m)
                if g is None or g is rep:
                    continue
                g.children = [c for c in g.children if c is not m]
                m.parent_ids = [rep.id]
                m.attributes["Parent"] = [rep.id]
                m.parents = [rep]
                rep.children.append(m)
                touched.add(g.id)
            rep.spans = [Span(seqid, u_lo, u_hi, strand)]
            changes.append(Change("merge-loci", rep.id or "?",
                                  f"merged {len(genes)} loci into {rep.id!r} "
                                  f"({len(members)} mRNAs, {seqid}:{u_lo}..{u_hi})"))
    dead = {gid for gid in touched
            if doc.feature_index.get(gid) is not None and not doc.feature_index[gid].children}
    if dead:
        doc.features = [f for f in doc.features
                        if not (f.type == "gene" and f.id in dead)]
        for gid in dead:
            doc.feature_index.pop(gid, None)
    return changes
```

- [ ] **Step 4: Register the pass**

In `src/ddbj_gff/normalize/normalize.py`, add `pass_merge_overlapping_loci` to the `.passes` import, insert it into `ALL_PASSES` immediately after `pass_wrap_cds_in_mrna`, and add `"merge-loci"` to `_APPLIED`:
```python
from .passes import (NormalizeContext, pass_directives, pass_coerce_transcript_to_mrna,
                     pass_wrap_cds_in_mrna, pass_merge_overlapping_loci, pass_circular_origin,
                     pass_trans_splicing_location, pass_so_terms, pass_transl_except, pass_anticodon)

ALL_PASSES = [pass_directives, pass_coerce_transcript_to_mrna, pass_wrap_cds_in_mrna,
              pass_merge_overlapping_loci, pass_circular_origin, pass_trans_splicing_location,
              pass_so_terms, pass_transl_except, pass_anticodon]

_APPLIED = {"add-directive", "rename-type", "add-qualifier", "add-child-feature", "merge-loci"}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_normalize_merge_loci.py -v'`
Expected: PASS (all 6). If `test_merge_two_overlapping_same_strand` fails on `dangling-parent`, the reparent (mRNA `parent_ids`/`attributes["Parent"]`) or the merged-gene removal is off — fix the pass, not the test.

- [ ] **Step 6: Full ddbj-gff regression (not-slow)**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests -m "not slow" -q'`
Expected: PASS (prior suite + the 6 new tests). Zero failures — the flag defaults off, so all existing normalize/roundtrip/flatfile tests are unchanged.

- [ ] **Step 7: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/gff_submission
git add src/ddbj_gff/normalize/passes.py src/ddbj_gff/normalize/normalize.py tests/test_normalize_merge_loci.py
git commit -m "feat(normalize): pass_merge_overlapping_loci (same-strand, % threshold, trans-splicing exempt)"
```

---

## Self-Review

**Spec coverage:** opt-in flag + default off → Task 1 + Task 2 registration; same-strand grouping → pass step (group by `(seqid, strand)`); % threshold `overlap/min(len)` → edge condition; connected components (transitive) → union-find; representative = lowest-coord gene, reparent, union span, remove empty genes → merge block; trans-spliced exempt → `_is_trans` skip; graph maintenance (feature_index/children/parents, no dangling) → reparent + `dead` removal + the `dangling-parent` assertion. Verification cases (merge, flag-off, opposite-strand, transitive, below-threshold, trans-exempt) → the 6 tests. No-regression → Task 2 Step 6.

**Placeholder scan:** No TBD/TODO. Every code step complete. The one contingency note (Task 2 Step 5) names the concrete failure mode to check.

**Type consistency:** `pass_merge_overlapping_loci(doc, ctx) -> list[Change]` matches sibling passes and `ALL_PASSES`. `NormalizeConfig.merge_overlapping_loci: bool` / `merge_overlap_min_fraction: float` used identically in the pass (`getattr(ctx.config, ...)`). `Span(seqid, u_lo, u_hi, strand)` positional form matches `model.Span`. `Change(action, target, message)` matches `.report`. `Feature.is_trans_spliced` reads the `exception` attribute (present pre-`pass_trans_splicing_location`, which runs after this pass). `_find(parent, x)` is module-level, no loop-closure issue.
```
