# GFF repair / curation layer — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `ddbj_gff.repair` subpackage — modular, individually-invokable GFF→GFF operations (internal-stop CDS → `misc_feature`, UTR-absent → partial mRNA, missing start/stop codon → partial CDS) as two-phase **detect → apply** units in a registry, usable by both a human (CLI) and an AI agent.

**Architecture:** A registry maps a stable operation `name` to an `Operation` (a `detect` that returns non-destructive `Candidate`s and an `apply` that mutates the `GffDocument` and returns `Change`s). A CLI (`python -m ddbj_gff.repair`) exposes `--list` / `--detect` / `--apply`. Partiality is encoded with INSDC explicit attributes (`partial=true`, `start_range`, `end_range`); translation reuses a copy of `gff2mss/translate.py`.

**Tech Stack:** Python 3.11+, biopython (`Bio.Seq`, `Bio.SeqFeature`, `Bio.Data.CodonTable`), hatchling package, pytest. Existing modules reused: `ddbj_gff.parse`, `ddbj_gff.writer.write`, `ddbj_gff.io.open_text`, `ddbj_gff.model` (`Feature`/`Span`/`GffDocument`), `ddbj_gff.normalize.report.Change`, `ddbj_gff.aa_names`.

## Global Constraints

- Python `>=3.11`; dependency floor `biopython>=1.83` (no new runtime deps).
- Design source of truth: `docs/superpowers/specs/2026-07-23-gff-repair-curation-layer-design.md`.
- GFF I/O is ASCII: parse with `ddbj_gff.parse`, write with `ddbj_gff.writer.write`.
- Partiality encoding (INSDC GFF3 v0.5, verbatim from `docs/INSDC GFF3 Specification - v0.5.docx`): `partial` value is always `true`; `start_range` applies to column 4 (start), `end_range` to column 5 (end); each value is two integers (or `.` for unknown) comma-separated. Convention used here: start-partial → `start_range=.,<col4>`; end-partial → `end_range=<col5>,.`.
- 5′/3′ ↔ column mapping is strand-aware: on `+`/`.`/`?` strand 5′=col4, 3′=col5; on `-` strand 5′=col5, 3′=col4. (Mirror of `gff2mss.convert.mrna_partial_flags`.)
- `misc_feature` conversion retypes the **CDS only** and adds a `Note`; the enclosing `gene`/`mRNA` and all parent/child links are left intact.
- Follow existing patterns: pass/registry style like `ddbj_gff.normalize`, module-CLI like `ddbj_gff.normalize.cli` / `ddbj_gff.validate.cli`.
- Run tests with the `mss_tools` env: `/lustre9/open/home/yt/micromamba/envs/mss_tools/bin/pytest` (this plan writes `pytest` for brevity — use that binary). New tests live in `tests/` alongside existing ones. Default run excludes `slow` (pyproject `addopts = -m 'not slow'`).
- TDD, one behavior per test, frequent commits. Commit message trailer:
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>` and the `Claude-Session:` line as in repo history.

---

### Task 1: Scaffolding — context, report, registry

**Files:**
- Create: `src/ddbj_gff/repair/__init__.py`
- Create: `src/ddbj_gff/repair/context.py`
- Create: `src/ddbj_gff/repair/report.py`
- Create: `src/ddbj_gff/repair/registry.py`
- Test: `tests/test_repair_registry.py`

**Interfaces:**
- Produces:
  - `RepairContext(sequences: dict[str, Bio.Seq.Seq] | None = None, transl_table: int = 1)` (`context.py`)
  - `Candidate(operation: str, feature_id: str | None, seqid: str, detail: str, payload: dict = {})` and `candidates_to_json(list[Candidate]) -> str`, `render_candidates(list[Candidate]) -> str` (`report.py`)
  - `Operation(name, summary, requires_sequence, detect, apply)`, module-global `REGISTRY: dict[str, Operation]`, `register(op) -> Operation`, `get_operation(name) -> Operation`, `list_operations() -> list[Operation]` (`registry.py`)
  - `apply` uses `ddbj_gff.normalize.report.Change` for its return type.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_repair_registry.py
from ddbj_gff.repair.context import RepairContext
from ddbj_gff.repair.report import Candidate, candidates_to_json, render_candidates
from ddbj_gff.repair.registry import Operation, register, get_operation, list_operations, REGISTRY


def test_context_defaults():
    ctx = RepairContext()
    assert ctx.sequences is None
    assert ctx.transl_table == 1


def test_candidate_json_and_render():
    c = Candidate("op-x", "gene1", "chr1", "would do X", payload={"five": True})
    import json
    parsed = json.loads(candidates_to_json([c]))
    assert parsed == [{"operation": "op-x", "feature_id": "gene1",
                       "seqid": "chr1", "detail": "would do X",
                       "payload": {"five": True}}]
    text = render_candidates([c])
    assert "op-x" in text and "gene1" in text


def test_register_and_lookup():
    op = Operation("op-test", "test op", requires_sequence=False,
                   detect=lambda doc, ctx: [], apply=lambda doc, ctx, sel: [])
    register(op)
    assert get_operation("op-test") is op
    assert "op-test" in {o.name for o in list_operations()}
    del REGISTRY["op-test"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_repair_registry.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'ddbj_gff.repair'`

- [ ] **Step 3: Create the four modules**

```python
# src/ddbj_gff/repair/context.py
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class RepairContext:
    sequences: dict | None = None   # seqid -> Bio.Seq.Seq (nucleotide)
    transl_table: int = 1           # default table when a CDS omits transl_table
```

```python
# src/ddbj_gff/repair/report.py
from __future__ import annotations

import json
from dataclasses import dataclass, field


@dataclass
class Candidate:
    operation: str
    feature_id: str | None
    seqid: str
    detail: str
    payload: dict = field(default_factory=dict)


def candidates_to_json(cands: list[Candidate]) -> str:
    return json.dumps([
        {"operation": c.operation, "feature_id": c.feature_id, "seqid": c.seqid,
         "detail": c.detail, "payload": c.payload}
        for c in cands
    ])


def render_candidates(cands: list[Candidate]) -> str:
    lines = [f"repair: {len(cands)} candidate(s)"]
    for c in cands:
        lines.append(f"  [{c.operation}] {c.feature_id} ({c.seqid}): {c.detail}")
    return "\n".join(lines) + "\n"
```

```python
# src/ddbj_gff/repair/registry.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from ..model import GffDocument
from .context import RepairContext
from .report import Candidate


@dataclass
class Operation:
    name: str
    summary: str
    requires_sequence: bool
    detect: Callable[[GffDocument, RepairContext], list[Candidate]]
    apply: Callable[[GffDocument, RepairContext, "list[Candidate] | None"], list]


REGISTRY: dict[str, Operation] = {}


def register(op: Operation) -> Operation:
    REGISTRY[op.name] = op
    return op


def get_operation(name: str) -> Operation:
    try:
        return REGISTRY[name]
    except KeyError:
        raise KeyError(f"unknown repair operation {name!r}; "
                       f"available: {sorted(REGISTRY)}") from None


def list_operations() -> list[Operation]:
    return list(REGISTRY.values())
```

```python
# src/ddbj_gff/repair/__init__.py
"""GFF repair / curation layer.

Modular, individually-invokable GFF->GFF operations. Each Operation is a
two-phase detect(doc, ctx)->list[Candidate] (non-destructive) + apply(doc, ctx,
selection)->list[Change] (mutating) unit registered in REGISTRY. Add a new
operation by writing its detect/apply and calling register() in operations.py.
"""
from __future__ import annotations

from .context import RepairContext
from .registry import Operation, REGISTRY, register, get_operation, list_operations
from .report import Candidate, candidates_to_json, render_candidates

__all__ = ["RepairContext", "Operation", "REGISTRY", "register", "get_operation",
           "list_operations", "Candidate", "candidates_to_json", "render_candidates"]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_repair_registry.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/repair/ tests/test_repair_registry.py
git commit -m "feat(repair): scaffold context/report/registry for GFF curation layer"
```

---

### Task 2: Partiality attribute helper (`partial.py`)

**Files:**
- Create: `src/ddbj_gff/repair/partial.py`
- Test: `tests/test_repair_partial.py`

**Interfaces:**
- Produces:
  - `is_partial(feature) -> bool` — True when the feature already carries `partial=true`.
  - `partial_attrs(five_prime: bool, three_prime: bool, strand: str, start: int, end: int) -> dict[str, list[str]]` — INSDC attributes for the given partial ends (`start`/`end` are the feature's col4/col5).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_repair_partial.py
from ddbj_gff.repair.partial import is_partial, partial_attrs
from ddbj_gff.model import Feature


def test_is_partial():
    f = Feature("x", "src", "mRNA", [], {"partial": ["true"]})
    assert is_partial(f) is True
    assert is_partial(Feature("y", "src", "mRNA", [], {})) is False


def test_partial_attrs_plus_5prime():
    # + strand, 5' partial -> start_range on col4
    assert partial_attrs(True, False, "+", 100, 500) == {
        "partial": ["true"], "start_range": [".,100"]}


def test_partial_attrs_plus_3prime():
    assert partial_attrs(False, True, "+", 100, 500) == {
        "partial": ["true"], "end_range": ["500,."]}


def test_partial_attrs_minus_5prime_maps_to_end():
    # - strand: 5' is genomic right (col5) -> end_range
    assert partial_attrs(True, False, "-", 100, 500) == {
        "partial": ["true"], "end_range": ["500,."]}


def test_partial_attrs_both_ends():
    assert partial_attrs(True, True, "+", 100, 500) == {
        "partial": ["true"], "start_range": [".,100"], "end_range": ["500,."]}


def test_partial_attrs_none():
    assert partial_attrs(False, False, "+", 100, 500) == {}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_repair_partial.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'ddbj_gff.repair.partial'`

- [ ] **Step 3: Implement `partial.py`**

```python
# src/ddbj_gff/repair/partial.py
from __future__ import annotations


def is_partial(feature) -> bool:
    return feature.attributes.get("partial") == ["true"]


def partial_attrs(five_prime: bool, three_prime: bool, strand: str,
                  start: int, end: int) -> dict[str, list[str]]:
    """INSDC partial attributes for the given partial ends.

    5' maps to col4 (start) on +/./? strand and to col5 (end) on - strand;
    3' maps to the other. start_range applies to col4, end_range to col5.
    Value form: '.,<col4>' for start_range, '<col5>,.' for end_range.
    """
    if strand == "-":
        start_partial, end_partial = three_prime, five_prime
    else:
        start_partial, end_partial = five_prime, three_prime
    attrs: dict[str, list[str]] = {}
    if start_partial or end_partial:
        attrs["partial"] = ["true"]
    if start_partial:
        attrs["start_range"] = [f".,{start}"]
    if end_partial:
        attrs["end_range"] = [f"{end},."]
    return attrs
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_repair_partial.py -v`
Expected: PASS (6 tests)

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/repair/partial.py tests/test_repair_partial.py
git commit -m "feat(repair): strand-aware INSDC partial attribute helper"
```

---

### Task 3: Copy translation engine + CDS translation adapter

**Files:**
- Create: `src/ddbj_gff/repair/translate.py` (verbatim copy of `../ddbj_mss_tools/src/gff2mss/translate.py`)
- Create: `src/ddbj_gff/repair/cds.py` (adapter from model `Feature` → protein / coding sequence)
- Test: `tests/test_repair_cds.py`

**Interfaces:**
- Consumes: `translate.translate_cds_with_transl_except(sf, parent_seq)` (from the copied file); `ddbj_gff.aa_names`.
- Produces (`cds.py`):
  - `collect_transl_excepts(cds) -> list[str]` — gather `transl_except` specs from the CDS attribute and `recoded_codon`/`stop_codon` children.
  - `coding_sequence(cds, ctx) -> tuple[str, object]` — returns `(coding_after_codon_start_upper, codon_table)` using `ctx.sequences[seqid]` and `Bio.Data.CodonTable.ambiguous_generic_by_id[table]`.
  - `protein_of(cds, ctx) -> str` — translated protein (transl_except-aware), no forced trailing-stop stripping beyond biopython defaults.
  - `has_internal_stop(protein: str) -> bool` — `*` present anywhere except a trailing stop.

- [ ] **Step 1: Copy the translation engine (no test yet — done in Step 3's test)**

```bash
cp ../ddbj_mss_tools/src/gff2mss/translate.py src/ddbj_gff/repair/translate.py
```

The copied file already carries the NIG provenance header (lines 1–3). Leave it verbatim.

- [ ] **Step 2: Write the failing test**

```python
# tests/test_repair_cds.py
from Bio.Seq import Seq
from ddbj_gff.model import Feature, Span
from ddbj_gff.repair.context import RepairContext
from ddbj_gff.repair.cds import (collect_transl_excepts, protein_of,
                                  has_internal_stop, coding_sequence)


def _cds(seqid, start, end, strand="+", phase=0, **attrs):
    a = {"ID": ["c1"]}
    a.update({k: (v if isinstance(v, list) else [v]) for k, v in attrs.items()})
    return Feature("c1", "src", "CDS", [Span(seqid, start, end, strand, phase)], a)


def test_has_internal_stop():
    assert has_internal_stop("MKV*") is False        # trailing stop only
    assert has_internal_stop("MK*V") is True          # internal stop
    assert has_internal_stop("MKV") is False


def test_protein_of_clean_cds():
    # ATG AAA GTT TAA -> M K V (stop stripped by translate cds path)
    seq = Seq("ATGAAAGTTTAA")
    ctx = RepairContext(sequences={"s": seq}, transl_table=1)
    cds = _cds("s", 1, 12, "+", 0, transl_table="1")
    prot = protein_of(cds, ctx)
    assert prot.startswith("MKV")
    assert has_internal_stop(prot) is False


def test_protein_of_internal_stop():
    # ATG TAA AAA TAA -> M * K (internal stop at aa index 1)
    seq = Seq("ATGTAAAAATAA")
    ctx = RepairContext(sequences={"s": seq}, transl_table=1)
    cds = _cds("s", 1, 12, "+", 0, transl_table="1")
    prot = protein_of(cds, ctx)
    assert has_internal_stop(prot) is True


def test_coding_sequence_respects_codon_start():
    seq = Seq("GATGAAAGTTTAA")   # leading G, codon_start=2
    ctx = RepairContext(sequences={"s": seq})
    cds = _cds("s", 1, 13, "+", 1)   # phase 1 -> codon_start 2
    coding, table = coding_sequence(cds, ctx)
    assert coding.startswith("ATGAAAGTT")
```

- [ ] **Step 3: Implement `cds.py`**

```python
# src/ddbj_gff/repair/cds.py
from __future__ import annotations

from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature

from .. import aa_names
from . import translate as _translate


def collect_transl_excepts(cds) -> list[str]:
    """transl_except specs from the CDS attribute and recoded_codon/stop_codon children."""
    specs = list(cds.attributes.get("transl_except", []))
    for child in cds.children:
        if child.type not in ("recoded_codon", "stop_codon"):
            continue
        sp = child.spans[0]
        loc = f"{sp.start}..{sp.end}" if sp.start != sp.end else f"{sp.start}"
        if sp.strand == "-":
            loc = f"complement({loc})"
        if child.type == "stop_codon":
            aa = "Term"
        else:
            aa = aa_names.to_abbrev((child.attributes.get("codon_redefined") or [""])[0])
        specs.append(f"(pos:{loc},aa:{aa})")
    return specs


def _table_id(cds, ctx) -> int:
    return cds.transl_table or ctx.transl_table


def coding_sequence(cds, ctx):
    """Return (coding_sequence_after_codon_start_upper, codon_table)."""
    seqid = cds.spans[0].seqid
    parent = ctx.sequences[seqid]
    coding = str(cds.to_biopython_location().extract(parent)).upper()
    codon_start = cds.codon_start or 1
    table = CodonTable.ambiguous_generic_by_id[int(_table_id(cds, ctx))]
    return coding[codon_start - 1:], table


def protein_of(cds, ctx) -> str:
    """Translate a CDS (transl_except-aware) to protein, trailing stop stripped."""
    table_id = _table_id(cds, ctx)
    excepts = collect_transl_excepts(cds)
    if excepts:
        sf = SeqFeature(cds.to_biopython_location(), type="CDS",
                        qualifiers={"transl_table": [str(table_id)],
                                    "codon_start": [str(cds.codon_start or 1)],
                                    "transl_except": excepts})
        return str(_translate.translate_cds_with_transl_except(sf, ctx.sequences[cds.spans[0].seqid]))
    coding, _ = coding_sequence(cds, ctx)
    protein = str(Seq(coding).translate(table=int(table_id)))
    return protein[:-1] if protein.endswith("*") else protein


def has_internal_stop(protein: str) -> bool:
    body = protein[:-1] if protein.endswith("*") else protein
    return "*" in body
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_repair_cds.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/repair/translate.py src/ddbj_gff/repair/cds.py tests/test_repair_cds.py
git commit -m "feat(repair): copy translate engine + CDS protein adapter"
```

---

### Task 4: Operation `utr-absent-to-partial-mrna` (structural)

**Files:**
- Create: `src/ddbj_gff/repair/operations.py`
- Modify: `src/ddbj_gff/repair/__init__.py` (import `operations` so it registers)
- Test: `tests/test_repair_op_utr.py`

**Interfaces:**
- Consumes: `Candidate` (Task 1), `partial_attrs`/`is_partial` (Task 2), `register`/`Operation` (Task 1), `Change` from `ddbj_gff.normalize.report`.
- Produces: registered operation `"utr-absent-to-partial-mrna"` (`requires_sequence=False`) with `detect`/`apply`. Payload keys: `{"five": bool, "three": bool, "strand": str, "start": int, "end": int}` (start/end = mRNA col4/col5).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_repair_op_utr.py
from ddbj_gff import parse
from ddbj_gff.repair.context import RepairContext
from ddbj_gff.repair.registry import get_operation

HDR = "##gff-version 3\n##sequence-region s 1 100000\n"

# mRNA 1..1000; exon 1..1000; CDS 1..600  -> 3' UTR present (601..1000),
# 5' UTR absent (exon_lo==cds_lo==1) -> 5' partial on + strand.
GFF = HDR + "\n".join([
    "s\tsrc\tgene\t1\t1000\t.\t+\t.\tID=g1;locus_tag=X_0001",
    "s\tsrc\tmRNA\t1\t1000\t.\t+\t.\tID=m1;Parent=g1",
    "s\tsrc\texon\t1\t1000\t.\t+\t.\tID=e1;Parent=m1",
    "s\tsrc\tCDS\t1\t600\t.\t+\t0\tID=c1;Parent=m1;transl_table=1",
]) + "\n"


def test_detect_finds_five_prime_partial():
    doc = parse(GFF)
    op = get_operation("utr-absent-to-partial-mrna")
    cands = op.detect(doc, RepairContext())
    assert len(cands) == 1
    c = cands[0]
    assert c.feature_id == "m1"
    assert c.payload["five"] is True and c.payload["three"] is False


def test_apply_sets_partial_attrs_on_mrna():
    doc = parse(GFF)
    op = get_operation("utr-absent-to-partial-mrna")
    changes = op.apply(doc, RepairContext(), None)
    m1 = doc.feature_index["m1"]
    assert m1.attributes.get("partial") == ["true"]
    assert m1.attributes.get("start_range") == [".,1"]
    assert "end_range" not in m1.attributes
    assert len(changes) == 1


def test_apply_is_idempotent():
    doc = parse(GFF)
    op = get_operation("utr-absent-to-partial-mrna")
    op.apply(doc, RepairContext(), None)
    assert op.detect(doc, RepairContext()) == []   # already partial -> no candidate
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_repair_op_utr.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'ddbj_gff.repair.operations'` (import via get_operation raises KeyError once module missing; the import error surfaces first)

- [ ] **Step 3: Implement `operations.py` (first operation) and wire the import**

```python
# src/ddbj_gff/repair/operations.py
from __future__ import annotations

from ..normalize.report import Change
from .context import RepairContext
from .registry import Operation, register
from .report import Candidate
from .partial import is_partial, partial_attrs


def _child_spans(mrna, ftype):
    return [s for c in mrna.children if c.type == ftype for s in c.spans]


# --- utr-absent-to-partial-mrna -------------------------------------------

def _detect_utr(doc, ctx: RepairContext) -> list[Candidate]:
    out: list[Candidate] = []
    for f in doc.features:
        if f.type != "mRNA" or is_partial(f):
            continue
        exon = _child_spans(f, "exon")
        cds = _child_spans(f, "CDS")
        if not exon or not cds:
            continue
        exon_lo, exon_hi = min(s.start for s in exon), max(s.end for s in exon)
        cds_lo, cds_hi = min(s.start for s in cds), max(s.end for s in cds)
        strand = exon[0].strand
        left_partial = exon_lo == cds_lo
        right_partial = exon_hi == cds_hi
        five, three = ((right_partial, left_partial) if strand == "-"
                       else (left_partial, right_partial))
        if not (five or three):
            continue
        m_lo = min(s.start for s in f.spans) if f.spans else exon_lo
        m_hi = max(s.end for s in f.spans) if f.spans else exon_hi
        ends = ", ".join(e for e, on in (("5'", five), ("3'", three)) if on)
        out.append(Candidate(
            "utr-absent-to-partial-mrna", f.id, exon[0].seqid,
            detail=f"mRNA {f.id!r} missing UTR on {ends} end -> mark partial",
            payload={"five": five, "three": three, "strand": strand,
                     "start": m_lo, "end": m_hi}))
    return out


def _apply_utr(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_utr(doc, ctx)
    changes: list[Change] = []
    for c in cands:
        f = doc.feature_index.get(c.feature_id)
        if f is None or f.type != "mRNA" or is_partial(f):
            continue
        p = c.payload
        attrs = partial_attrs(p["five"], p["three"], p["strand"], p["start"], p["end"])
        f.attributes.update(attrs)
        changes.append(Change("mark-partial", f.id or "?",
                              f"mRNA marked partial ({', '.join(sorted(attrs))})"))
    return changes


register(Operation("utr-absent-to-partial-mrna",
                   "Mark an mRNA partial on ends where a UTR is absent (structural).",
                   requires_sequence=False, detect=_detect_utr, apply=_apply_utr))
```

Then wire the import so registration happens on `import ddbj_gff.repair`:

```python
# src/ddbj_gff/repair/__init__.py  — append at the very end of the file
from . import operations as _operations  # noqa: E402,F401  (populates REGISTRY)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_repair_op_utr.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/repair/operations.py src/ddbj_gff/repair/__init__.py tests/test_repair_op_utr.py
git commit -m "feat(repair): op utr-absent-to-partial-mrna"
```

---

### Task 5: Operation `missing-start-stop-to-partial-cds` (sequence-based)

**Files:**
- Modify: `src/ddbj_gff/repair/operations.py` (add the operation + register)
- Test: `tests/test_repair_op_startstop.py`

**Interfaces:**
- Consumes: `coding_sequence` (Task 3), `partial_attrs`/`is_partial` (Task 2).
- Produces: registered operation `"missing-start-stop-to-partial-cds"` (`requires_sequence=True`). Payload keys: `{"five": bool, "three": bool, "strand": str, "start": int, "end": int}` (CDS col4/col5). Skips trans-spliced CDS and CDS whose seqid is absent from `ctx.sequences`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_repair_op_startstop.py
from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.repair.context import RepairContext
from ddbj_gff.repair.registry import get_operation

HDR = "##gff-version 3\n##sequence-region s 1 100000\n"

# CDS 1..12 on + strand. Sequence ATG AAA GTT AAA -> starts with ATG (start codon),
# does NOT end in a stop codon -> 3' partial only.
GFF = HDR + "\n".join([
    "s\tsrc\tgene\t1\t12\t.\t+\t.\tID=g1;locus_tag=X_0001",
    "s\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=m1;Parent=g1",
    "s\tsrc\tCDS\t1\t12\t.\t+\t0\tID=c1;Parent=m1;transl_table=1",
]) + "\n"
SEQ = {"s": Seq("ATGAAAGTTAAA")}


def _ctx():
    return RepairContext(sequences=SEQ, transl_table=1)


def test_detect_three_prime_partial_missing_stop():
    doc = parse(GFF)
    op = get_operation("missing-start-stop-to-partial-cds")
    cands = op.detect(doc, _ctx())
    assert len(cands) == 1
    assert cands[0].payload["five"] is False    # ATG is a start codon
    assert cands[0].payload["three"] is True     # AAA is not a stop


def test_detect_five_prime_partial_missing_start():
    # sequence starts with CTT (not a start codon) and ends with a stop (TAA)
    doc = parse(GFF)
    op = get_operation("missing-start-stop-to-partial-cds")
    ctx = RepairContext(sequences={"s": Seq("CTTAAAGTTTAA")}, transl_table=1)
    cands = op.detect(doc, ctx)
    assert cands[0].payload["five"] is True      # CTT not a start codon
    assert cands[0].payload["three"] is False    # TAA is a stop


def test_apply_sets_end_range_on_cds():
    doc = parse(GFF)
    op = get_operation("missing-start-stop-to-partial-cds")
    op.apply(doc, _ctx(), None)
    c1 = doc.feature_index["c1"]
    assert c1.attributes.get("partial") == ["true"]
    assert c1.attributes.get("end_range") == ["12,."]
    assert "start_range" not in c1.attributes


def test_requires_sequence_flag():
    assert get_operation("missing-start-stop-to-partial-cds").requires_sequence is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_repair_op_startstop.py -v`
Expected: FAIL — `KeyError: "unknown repair operation 'missing-start-stop-to-partial-cds'"`

- [ ] **Step 3: Add the operation to `operations.py`**

```python
# src/ddbj_gff/repair/operations.py — add imports at top
from .cds import coding_sequence
```

```python
# src/ddbj_gff/repair/operations.py — add before the final register() calls block

# --- missing-start-stop-to-partial-cds ------------------------------------

def _detect_startstop(doc, ctx: RepairContext) -> list[Candidate]:
    out: list[Candidate] = []
    if ctx.sequences is None:
        return out
    for f in doc.features:
        if f.type != "CDS" or is_partial(f) or not f.spans:
            continue
        if f.is_trans_spliced:
            continue
        seqid = f.spans[0].seqid
        if seqid not in ctx.sequences:
            continue
        coding, table = coding_sequence(f, ctx)
        codon_start = f.codon_start or 1
        strand = f.spans[0].strand
        first = coding[:3]
        last = coding[-3:]
        five = codon_start > 1 or len(first) < 3 or first not in table.start_codons
        three = len(last) < 3 or last not in table.stop_codons
        if not (five or three):
            continue
        c_lo = min(s.start for s in f.spans)
        c_hi = max(s.end for s in f.spans)
        ends = ", ".join(e for e, on in (("5' (no start codon)", five),
                                         ("3' (no stop codon)", three)) if on)
        out.append(Candidate(
            "missing-start-stop-to-partial-cds", f.id, seqid,
            detail=f"CDS {f.id!r} partial: {ends}",
            payload={"five": five, "three": three, "strand": strand,
                     "start": c_lo, "end": c_hi}))
    return out


def _apply_startstop(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_startstop(doc, ctx)
    changes: list[Change] = []
    for c in cands:
        f = doc.feature_index.get(c.feature_id)
        if f is None or f.type != "CDS" or is_partial(f):
            continue
        p = c.payload
        attrs = partial_attrs(p["five"], p["three"], p["strand"], p["start"], p["end"])
        f.attributes.update(attrs)
        changes.append(Change("mark-partial", f.id or "?",
                              f"CDS marked partial ({', '.join(sorted(attrs))})"))
    return changes


register(Operation("missing-start-stop-to-partial-cds",
                   "Mark a CDS partial when its sequence lacks a start or stop codon.",
                   requires_sequence=True, detect=_detect_startstop, apply=_apply_startstop))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_repair_op_startstop.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/repair/operations.py tests/test_repair_op_startstop.py
git commit -m "feat(repair): op missing-start-stop-to-partial-cds"
```

---

### Task 6: Operation `internal-stop-to-misc` + validate accepts `misc_feature`

**Files:**
- Modify: `src/ddbj_gff/repair/operations.py` (add the operation + register)
- Modify: `src/ddbj_gff/validate/rules.py:84` (add `misc_feature` to `_INSDC_GFF3_SPECIAL`)
- Test: `tests/test_repair_op_internalstop.py`
- Test: `tests/test_validate_misc_feature.py`

**Interfaces:**
- Consumes: `protein_of`/`has_internal_stop` (Task 3).
- Produces: registered operation `"internal-stop-to-misc"` (`requires_sequence=True`). `apply` retypes the CDS to `misc_feature` and appends a `Note`; skips trans-spliced CDS and CDS whose seqid is absent from `ctx.sequences`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_repair_op_internalstop.py
from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.repair.context import RepairContext
from ddbj_gff.repair.registry import get_operation

HDR = "##gff-version 3\n##sequence-region s 1 100000\n"
# CDS 1..12: ATG TAA AAA TAA -> M * K  (internal stop at aa index 1)
GFF = HDR + "\n".join([
    "s\tsrc\tgene\t1\t12\t.\t+\t.\tID=g1;locus_tag=X_0001",
    "s\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=m1;Parent=g1",
    "s\tsrc\tCDS\t1\t12\t.\t+\t0\tID=c1;Parent=m1;transl_table=1",
]) + "\n"
SEQ = {"s": Seq("ATGTAAAAATAA")}


def _ctx():
    return RepairContext(sequences=SEQ, transl_table=1)


def test_detect_internal_stop():
    doc = parse(GFF)
    op = get_operation("internal-stop-to-misc")
    cands = op.detect(doc, _ctx())
    assert len(cands) == 1 and cands[0].feature_id == "c1"


def test_apply_retypes_cds_to_misc_feature_with_note():
    doc = parse(GFF)
    op = get_operation("internal-stop-to-misc")
    op.apply(doc, _ctx(), None)
    c1 = doc.feature_index["c1"]
    assert c1.type == "misc_feature"
    assert any("internal stop" in n for n in c1.note)
    # gene/mRNA and links intact
    g1 = doc.feature_index["g1"]
    m1 = doc.feature_index["m1"]
    assert g1.type == "gene" and m1.type == "mRNA"
    assert c1 in m1.children


def test_no_candidate_for_clean_cds():
    doc = parse(GFF)
    op = get_operation("internal-stop-to-misc")
    ctx = RepairContext(sequences={"s": Seq("ATGAAAGTTTAA")}, transl_table=1)  # M K V
    assert op.detect(doc, ctx) == []
```

```python
# tests/test_validate_misc_feature.py
from ddbj_gff import parse
from ddbj_gff.validate import validate

HDR = "##gff-version 3\n##sequence-region s 1 100000\n"
GFF = HDR + "s\tsrc\tmisc_feature\t1\t12\t.\t+\t.\tID=x;Note=nonfunctional\n"


def test_misc_feature_is_accepted_type():
    diags = validate(parse(GFF))
    assert not any(d.code == "feature-type-not-insdc" for d in diags)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_repair_op_internalstop.py tests/test_validate_misc_feature.py -v`
Expected: FAIL — internal-stop op unknown (`KeyError`); `test_misc_feature_is_accepted_type` fails because `misc_feature` currently triggers `feature-type-not-insdc`.

- [ ] **Step 3: Add the operation and update the validate special-case set**

```python
# src/ddbj_gff/repair/operations.py — add import at top
from .cds import protein_of, has_internal_stop
```

```python
# src/ddbj_gff/repair/operations.py — add before the register() call for this op

# --- internal-stop-to-misc -------------------------------------------------

def _detect_internal_stop(doc, ctx: RepairContext) -> list[Candidate]:
    out: list[Candidate] = []
    if ctx.sequences is None:
        return out
    for f in doc.features:
        if f.type != "CDS" or not f.spans or f.is_trans_spliced:
            continue
        if f.spans[0].seqid not in ctx.sequences:
            continue
        if has_internal_stop(protein_of(f, ctx)):
            out.append(Candidate(
                "internal-stop-to-misc", f.id, f.spans[0].seqid,
                detail=f"CDS {f.id!r} has an internal stop codon -> misc_feature",
                payload={}))
    return out


def _apply_internal_stop(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_internal_stop(doc, ctx)
    changes: list[Change] = []
    for c in cands:
        f = doc.feature_index.get(c.feature_id)
        if f is None or f.type != "CDS":
            continue
        f.type = "misc_feature"
        note = (f"nonfunctional CDS: internal stop codon(s) detected in {f.id}; "
                f"not translated")
        f.attributes.setdefault("Note", []).append(note)
        changes.append(Change("rename-type", f.id or "?", f"CDS -> misc_feature ({note})"))
    return changes


register(Operation("internal-stop-to-misc",
                   "Retype a CDS with an internal stop codon to misc_feature (+Note).",
                   requires_sequence=True, detect=_detect_internal_stop,
                   apply=_apply_internal_stop))
```

```python
# src/ddbj_gff/validate/rules.py:84  — add misc_feature
_INSDC_GFF3_SPECIAL = {"recoded_codon", "anticodon", "stop_codon", "start_codon",
                       "misc_feature"}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_repair_op_internalstop.py tests/test_validate_misc_feature.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/repair/operations.py src/ddbj_gff/validate/rules.py \
        tests/test_repair_op_internalstop.py tests/test_validate_misc_feature.py
git commit -m "feat(repair): op internal-stop-to-misc; validate accepts misc_feature"
```

---

### Task 7: Driver + CLI (`python -m ddbj_gff.repair`)

**Files:**
- Create: `src/ddbj_gff/repair/driver.py`
- Create: `src/ddbj_gff/repair/cli.py`
- Create: `src/ddbj_gff/repair/__main__.py`
- Modify: `src/ddbj_gff/repair/__init__.py` (export `run_detect`, `run_apply`, `DEFAULT_ORDER`)
- Test: `tests/test_repair_cli.py`
- Create (fixtures): `tests/fixtures/repair_internal_stop.gff3`, `tests/fixtures/repair_internal_stop.fasta`

**Interfaces:**
- Consumes: all registered operations; `ddbj_gff.parse`, `ddbj_gff.writer.write`, `ddbj_gff.io.open_text`, `ddbj_gff.validate.validate`.
- Produces (`driver.py`):
  - `DEFAULT_ORDER = ["internal-stop-to-misc", "utr-absent-to-partial-mrna", "missing-start-stop-to-partial-cds"]`
  - `run_detect(doc, ctx, names: list[str] | None = None) -> list[Candidate]`
  - `run_apply(doc, ctx, names: list[str]) -> list[Change]` — applies in the given order; the CLI resolves `all` to `DEFAULT_ORDER` (registry-intersected).
  - `load_sequences(path) -> dict[str, Seq]`
- CLI (`main(argv) -> int`): `--list`; `--gff IN [--fasta FA] --detect [--only a,b]`; `--gff IN [--fasta FA] --apply a,b|all --out OUT [--report R]`.

- [ ] **Step 1: Write the failing test + fixtures**

```
# tests/fixtures/repair_internal_stop.gff3
##gff-version 3
##sequence-region s 1 12
s	src	gene	1	12	.	+	.	ID=g1;locus_tag=X_0001
s	src	mRNA	1	12	.	+	.	ID=m1;Parent=g1
s	src	CDS	1	12	.	+	0	ID=c1;Parent=m1;transl_table=1
```

```
# tests/fixtures/repair_internal_stop.fasta
>s
ATGTAAAAATAA
```

```python
# tests/test_repair_cli.py
import os
from ddbj_gff.repair.cli import main
from ddbj_gff import parse
from ddbj_gff.validate import validate

FIX = os.path.join(os.path.dirname(__file__), "fixtures")
GFF = os.path.join(FIX, "repair_internal_stop.gff3")
FASTA = os.path.join(FIX, "repair_internal_stop.fasta")


def test_list(capsys):
    rc = main(["--list"])
    out = capsys.readouterr().out
    assert rc == 0
    for name in ("internal-stop-to-misc", "utr-absent-to-partial-mrna",
                 "missing-start-stop-to-partial-cds"):
        assert name in out


def test_detect_json_does_not_write(capsys, tmp_path):
    rc = main(["--gff", GFF, "--fasta", FASTA, "--detect", "--json"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "internal-stop-to-misc" in out and "c1" in out


def test_apply_writes_curated_gff_that_validates(tmp_path):
    out = tmp_path / "out.gff3"
    rc = main(["--gff", GFF, "--fasta", FASTA, "--apply", "all", "--out", str(out)])
    assert rc == 0
    before = {d.code for d in validate(parse(open(GFF).read()))
              if d.severity.name == "ERROR"}
    doc = parse(out.read_text())
    assert doc.feature_index["c1"].type == "misc_feature"
    after = {d.code for d in validate(doc) if d.severity.name == "ERROR"}
    assert after <= before   # repair introduced no new validation ERROR
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_repair_cli.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'ddbj_gff.repair.cli'`

- [ ] **Step 3: Implement `driver.py`, `cli.py`, `__main__.py`, and exports**

```python
# src/ddbj_gff/repair/driver.py
from __future__ import annotations

from Bio import SeqIO

from ..io import open_text
from .context import RepairContext
from .registry import get_operation, list_operations
from .report import Candidate

DEFAULT_ORDER = ["internal-stop-to-misc", "utr-absent-to-partial-mrna",
                 "missing-start-stop-to-partial-cds"]


def load_sequences(path: str) -> dict:
    with open_text(path) as fh:
        return {rec.id: rec.seq for rec in SeqIO.parse(fh, "fasta")}


def run_detect(doc, ctx: RepairContext, names=None) -> list[Candidate]:
    ops = [get_operation(n) for n in names] if names else list_operations()
    out: list[Candidate] = []
    for op in ops:
        out.extend(op.detect(doc, ctx))
    return out


def run_apply(doc, ctx: RepairContext, names) -> list:
    changes: list = []
    for n in names:
        changes.extend(get_operation(n).apply(doc, ctx, None))
    return changes
```

```python
# src/ddbj_gff/repair/cli.py
from __future__ import annotations

import argparse
import sys

from .. import parse
from ..writer import write
from ..normalize.report import NormalizationReport
from .context import RepairContext
from .registry import list_operations, REGISTRY
from .report import candidates_to_json, render_candidates
from .driver import run_detect, run_apply, load_sequences, DEFAULT_ORDER


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="ddbj_gff.repair",
                                 description="Apply modular GFF curation operations")
    ap.add_argument("--list", action="store_true", help="list available operations")
    ap.add_argument("--gff")
    ap.add_argument("--fasta")
    ap.add_argument("--transl-table", type=int, default=1, dest="transl_table")
    ap.add_argument("--detect", action="store_true", help="preview candidates; write nothing")
    ap.add_argument("--json", action="store_true", help="detect output as JSON")
    ap.add_argument("--only", help="comma-separated operation names (detect)")
    ap.add_argument("--apply", help="comma-separated operation names, or 'all'")
    ap.add_argument("--out")
    ap.add_argument("--report")
    args = ap.parse_args(argv)

    if args.list:
        for op in list_operations():
            seq = " (needs FASTA)" if op.requires_sequence else ""
            sys.stdout.write(f"{op.name}{seq}: {op.summary}\n")
        return 0

    if not args.gff:
        ap.error("--gff is required unless --list")

    with open(args.gff, encoding="ascii", errors="replace") as fh:
        doc = parse(fh.read())
    sequences = load_sequences(args.fasta) if args.fasta else None
    ctx = RepairContext(sequences=sequences, transl_table=args.transl_table)

    if args.detect:
        names = args.only.split(",") if args.only else None
        cands = run_detect(doc, ctx, names)
        sys.stdout.write(candidates_to_json(cands) + "\n" if args.json
                         else render_candidates(cands))
        return 0

    if args.apply:
        if args.apply == "all":
            names = [n for n in DEFAULT_ORDER if n in REGISTRY]
        else:
            names = args.apply.split(",")
        changes = run_apply(doc, ctx, names)
        out_text = write(doc)
        if args.out:
            with open(args.out, "w", encoding="ascii") as fh:
                fh.write(out_text)
        else:
            sys.stdout.write(out_text)
        report = NormalizationReport(applied=changes, unresolved=[]).render()
        if args.report:
            with open(args.report, "w", encoding="ascii") as fh:
                fh.write(report)
        else:
            sys.stderr.write(report)
        return 0

    ap.error("nothing to do: pass --list, --detect, or --apply")
    return 2
```

```python
# src/ddbj_gff/repair/__main__.py
import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
```

```python
# src/ddbj_gff/repair/__init__.py — extend the exports (add after the operations import)
from .driver import run_detect, run_apply, load_sequences, DEFAULT_ORDER  # noqa: E402

__all__ += ["run_detect", "run_apply", "load_sequences", "DEFAULT_ORDER"]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_repair_cli.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Run the whole suite to confirm no regressions**

Run: `pytest -q`
Expected: all pass (existing 193 + new repair tests), no failures.

- [ ] **Step 6: Commit**

```bash
git add src/ddbj_gff/repair/driver.py src/ddbj_gff/repair/cli.py \
        src/ddbj_gff/repair/__main__.py src/ddbj_gff/repair/__init__.py \
        tests/test_repair_cli.py tests/fixtures/repair_internal_stop.gff3 \
        tests/fixtures/repair_internal_stop.fasta
git commit -m "feat(repair): driver + module CLI (--list/--detect/--apply)"
```

---

### Task 8: Docs — README/usage note for the repair layer

**Files:**
- Modify: `README.md` (add a short "GFF repair / curation" usage section)
- Modify: `docs/development-setup.md` (mention the new `python -m ddbj_gff.repair` module CLI alongside normalize/validate)

**Interfaces:** none (documentation only).

- [ ] **Step 1: Add a usage section to `README.md`**

Add, near the existing normalize/validate usage, a section documenting:

```markdown
### GFF repair / curation (`python -m ddbj_gff.repair`)

Modular, individually-invokable curation operations. List them, preview
(`--detect`), then apply selected ones (`--apply`):

```bash
python -m ddbj_gff.repair --list
python -m ddbj_gff.repair --gff in.gff --fasta seq.fasta --detect --json
python -m ddbj_gff.repair --gff in.gff --fasta seq.fasta --apply all --out out.gff
```

Operations: `internal-stop-to-misc` (CDS with an internal stop → `misc_feature`),
`utr-absent-to-partial-mrna` (missing UTR → partial mRNA), and
`missing-start-stop-to-partial-cds` (missing start/stop codon → partial CDS).
Partiality is written as INSDC `partial=true` + `start_range`/`end_range`.
```

- [ ] **Step 2: Add the module CLI to `docs/development-setup.md`**

In the module-CLI list (the `python -m ddbj_gff.normalize` / `.validate` mention near the top),
add `python -m ddbj_gff.repair` with a one-line description.

- [ ] **Step 3: Verify docs render (no code to test)**

Run: `python -m ddbj_gff.repair --list`
Expected: prints the three operation names and summaries (matches the README).

- [ ] **Step 4: Commit**

```bash
git add README.md docs/development-setup.md
git commit -m "docs: document ddbj_gff.repair curation CLI"
```

---

## Self-Review

**1. Spec coverage:**
- Placement (GFF→GFF layer in `ddbj_gff`) → Task 1 package scaffold. ✓
- Registry + detect/apply two-phase → Tasks 1, 4–6. ✓
- Human + agent consumers (CLI list/detect/apply, JSON + human report) → Task 7. ✓
- Partiality via INSDC `partial`/`start_range`/`end_range` → Task 2 helper, applied in Tasks 4–5. ✓
- `misc_feature` = CDS-only retype + Note, links intact → Task 6 (asserted in test). ✓
- Translation copied from gff2mss (provenance kept) → Task 3. ✓
- Initial three operations → Tasks 4, 5, 6. ✓
- Extensibility (add detect/apply + register) → Task 1 `__init__` docstring + Task 4 pattern. ✓
- Validate rule for partial attrs / `misc_feature` (spec's "to pin during planning") → Task 6 adds `misc_feature` to `_INSDC_GFF3_SPECIAL`; `partial`/`start_range`/`end_range` are unrecognised attributes that `validate` already leaves untouched (no rule needed) and are verified by Task 7's round-trip validate assertion. ✓
- Testing (per-op detect/apply, idempotence, round-trip, integration) → idempotence in Task 4, round-trip validate in Tasks 6–7. ✓

**2. Placeholder scan:** No TBD/TODO; every code step has complete code; every command has expected output. ✓

**3. Type consistency:** `Candidate` fields (`operation`, `feature_id`, `seqid`, `detail`, `payload`) consistent across report/registry/operations/driver. `partial_attrs(five, three, strand, start, end)` signature matches all call sites (Tasks 4, 5). `Operation(name, summary, requires_sequence, detect, apply)` matches all `register(...)` calls. `run_detect`/`run_apply`/`DEFAULT_ORDER` in `driver.py` match CLI usage. `coding_sequence`/`protein_of`/`has_internal_stop`/`collect_transl_excepts` in `cds.py` match Task 5/6 call sites. ✓

## Notes for the implementer

- Run every `pytest` via `/lustre9/open/home/yt/micromamba/envs/mss_tools/bin/pytest` (the env with `ddbj-gff`, `ddbj-mss-tools`, biopython, pytest).
- The three operations already exist *implicitly* inside `gff2mss/convert.py` (`mrna_partial_flags`, the `if "*" in body` internal-stop→`misc_feature` branch at ~line 293, and the start-codon check). This plan lifts that logic onto the GFF as explicit, selectable operations. When in doubt about detection semantics, that file is the reference implementation.
- `Diagnostic.severity` is an enum with `.name`; Task 7's round-trip test filters `d.severity.name == "ERROR"`. Confirm the enum member name in `ddbj_gff/errors.py` / `validate/severities.py` if the test needs adjusting.
