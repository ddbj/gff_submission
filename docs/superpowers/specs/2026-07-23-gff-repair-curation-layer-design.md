# GFF repair / curation layer — design

**Date:** 2026-07-23
**Status:** design (approved, for spec review)
**Scope:** new `repair` subpackage in `ddbj-gff` (`gff_submission`) — a GFF→GFF curation
layer of modular, individually-invokable operations that apply *judgment-bearing* fixes
(internal-stop CDS → `misc_feature`, UTR-absent → partial mRNA, missing start/stop codon →
partial CDS). Each operation is a two-phase **detect → apply** unit registered in a registry
so that both a human (via CLI) and an AI agent can inspect a GFF, choose which fixes to run,
and apply them selectively.

## Problem / motivation

Preparing a canonical INSDC GFF3 for DDBJ (MSS) submission requires curation decisions that are
**not** deterministic canonicalization and are not always safe to apply blindly:

- A CDS with an **internal stop codon** cannot be translated to a valid protein; INSDC practice
  is to demote it to a generic `misc_feature` with an explanatory note.
- An mRNA lacking a **5′ or 3′ UTR** on a given end implies that end is **partial** (the model
  is truncated), which must be marked so downstream tools do not assert a complete transcript.
- A CDS whose sequence lacks a **start** or **stop** codon is **partial** on the corresponding end.

Today these are handled inconsistently and *implicitly*:

- `gff2mss` (sibling `ddbj_mss_tools`) already **re-infers** UTR-based partiality at MSS-emit
  time (`convert.py:mrna_partial_flags`, emitting `<`/`>` in the MSS location). It is buried in
  the emitter, invisible on the GFF, and not reviewable before submission.
- `gff2mss/translate.py` can translate a CDS (with `transl_except`) but **does not** detect
  internal stops or convert anything to `misc_feature` — an internal stop simply appears as `*`.
- `ddbj_gff` itself has **no** translation and no partiality handling; `normalize` receives only
  sequence *lengths*, not the sequences. The canonical model's `Span` has no partial flag, and
  the INSDC explicit attributes `partial` / `start_range` / `end_range` are not consumed
  (see `docs/spec-v0.5-gap-review.md` A-7).

The goal is a **curation layer that lifts these decisions up onto the canonical GFF**, makes
them explicit and reviewable, and exposes each as a discrete operation an agent (or a human) can
run on demand — distinct from `normalize`, which is always-safe structural canonicalization.

## Decisions (confirmed)

1. **Placement — GFF→GFF curation layer in `ddbj_gff`.** New subpackage `src/ddbj_gff/repair/`.
   Input: a canonical GFF (+ FASTA). Output: a curated GFF with partiality / `misc_feature`
   encoded explicitly, plus a report. `gff2mss` will (as a later, separate change) prefer these
   explicit attributes over its own inference.
2. **Interface — registry of discrete operations, two-phase detect → apply.** Each operation is
   individually addressable by a stable `name`; both a human (CLI) and an agent use the same
   `detect` (preview, non-destructive) → `apply` (mutate) workflow.
3. **Human *and* agent are first-class consumers.** The layer ships a CLI (list / detect / apply)
   and machine-readable (JSON) + human-readable reports.
4. **Partiality representation — INSDC explicit attributes.** `partial=true` plus
   `start_range` / `end_range` per INSDC GFF3 v0.5 (`docs/INSDC GFF3 Specification - v0.5.docx`).
   These are ordinary GFF attributes, so parser/writer round-trip them with **no model change**;
   this also matches gap-review A-7's "consume explicit partial attributes".
5. **`misc_feature` conversion — CDS only.** Retype the offending **CDS** to `misc_feature` and
   add a `Note`; leave the enclosing `gene` / `mRNA` and all parent/child links intact. Minimal,
   reversible.
6. **Translation home — copy now, unify later.** Copy `gff2mss/translate.py` into
   `repair/translate.py` (keeping the NIG provenance header). Switching `gff2mss` to import the
   `ddbj_gff` copy is deferred to a separate task (the dependency direction gff2mss → ddbj_gff
   makes this the correct eventual home).
7. **Initial operation set — the three above.** The framework must make adding a fourth operation
   trivial (write `detect`/`apply`, register in `REGISTRY`).

## Architecture

### Package layout

```
src/ddbj_gff/repair/
├── __init__.py       # public API: REGISTRY, run helpers; contract docstring
├── registry.py       # Operation dataclass + REGISTRY (name -> Operation)
├── context.py        # RepairContext (sequences, transl_table)
├── report.py         # Candidate (detect output); reuse normalize Change for apply
├── operations.py     # the three operations (detect + apply functions)
├── partial.py        # shared 5'/3' <-> start_range/end_range helpers (strand-aware)
├── translate.py      # copied from gff2mss (NIG provenance header)
├── cli.py            # argparse front-end: --list / --detect / --apply
└── __main__.py       # `python -m ddbj_gff.repair`
```

### Operation abstraction

```python
@dataclass
class Operation:
    name: str                  # stable id, e.g. "internal-stop-to-misc"
    summary: str               # one-line human/AI description
    requires_sequence: bool    # detect/apply need FASTA?
    detect: Callable[[GffDocument, RepairContext], list[Candidate]]
    apply:  Callable[[GffDocument, RepairContext, list[Candidate] | None], list[Change]]

REGISTRY: dict[str, Operation]   # extension = add an entry here
```

- **`detect(doc, ctx)`** — non-destructive. Returns a `Candidate` per affected feature:
  `operation`, `feature_id`, `seqid`, human/AI-readable `detail` (why + what would change), and a
  `payload` dict carrying enough info for `apply` to act (e.g. which side is partial).
- **`apply(doc, ctx, selection)`** — if `selection` is a candidate list, apply exactly those;
  if `None`, run `detect` and apply all. Mutates `doc` in place, returns `list[Change]` (the same
  `Change` type as `normalize/report.py`, for report consistency).
- Idempotence: `detect` judges the *current* doc state and `apply` re-checks the feature type, so
  re-running after an apply produces no further candidates/changes.

### Context and sequence access

```python
@dataclass
class RepairContext:
    sequences: dict[str, Seq] | None   # seqid -> nucleotide Seq (from FASTA)
    transl_table: int = 1              # default table when a CDS omits transl_table
```

The CLI loads the FASTA into full `Seq` objects (not just lengths — new relative to `normalize`),
reusing `io.open_text`. `requires_sequence=True` operations error clearly if `sequences is None`.

## Operations (initial set)

Recommended default apply order: **`internal-stop-to-misc` first** (after retyping, the feature is
no longer a CDS and is naturally excluded from the partial-CDS operation), then the partial ops.

### `internal-stop-to-misc` (sequence-based)

- **detect:** for each `CDS`, translate with `transl_except` applied
  (`translate_cds_with_transl_except`); if a stop symbol `*` appears at any position **other than
  the final residue**, emit a candidate (payload records the internal-stop aa position(s)).
- **apply:** set `feature.type = "misc_feature"`; append a `Note`
  (e.g. `nonfunctional CDS: internal stop codon at aa <pos>`). Leave `gene`/`mRNA` and all
  parent/child links unchanged. No-op if the feature is no longer a `CDS`.

### `utr-absent-to-partial-mrna` (structural, no sequence)

- **detect:** for each `mRNA`, determine whether the CDS extent reaches the mRNA/exon boundary on
  each genomic end (UTR absent on that end), mapping genomic left/right to 5′/3′ by strand — same
  logic as `gff2mss:mrna_partial_flags` (`left_partial = exon_lo == cds_lo`, etc.). Emit a
  candidate when either end is partial and the mRNA is not already marked partial.
- **apply:** set `partial=true` on the mRNA and the appropriate `start_range` / `end_range`
  attribute(s) for the partial end(s), via the strand-aware helper in `partial.py`.

### `missing-start-stop-to-partial-cds` (sequence-based)

- **detect:** for each `CDS`, inspect the coding sequence boundaries using the FASTA + `transl_table`:
  first codon (after `codon_start`) not in the table's start codons → **5′ partial**
  (`codon_start > 1` is treated as 5′ partial as well); sequence does not end in a stop codon →
  **3′ partial**. Emit a candidate with the partial end(s).
- **apply:** set `partial=true` + `start_range` / `end_range` on the CDS (strand-aware helper).

### Partiality encoding helper (`partial.py`)

Maps a (5′-partial, 3′-partial) pair on a strand to INSDC attributes:
`partial=true`, and `start_range` / `end_range` on the correct genomic coordinate
(on `+`: 5′→start_range, 3′→end_range; on `-`: 5′→end_range, 3′→start_range). Shared by both
partial operations.

> **To pin during planning:** the exact attribute *value* syntax for `start_range` / `end_range`
> (e.g. the `.,N` / `N,.` range form) must be extracted verbatim from
> `docs/INSDC GFF3 Specification - v0.5.docx` before implementation, and a `validate` rule to
> accept these attributes should be added (gap-review A-7 notes validate currently has no
> `partial`-family rule).

## CLI (`python -m ddbj_gff.repair`)

Follows the existing module-CLI convention (`ddbj_gff.normalize` / `ddbj_gff.validate`).

- `--list` — print each operation's `name`, `summary`, `requires_sequence` (discovery for humans
  and agents).
- `--gff IN [--fasta FA] --detect [--only a,b]` — print the candidate report in **JSON and
  human-readable** form; write **no** GFF (pure preview).
- `--gff IN [--fasta FA] --apply a,b|all --out OUT [--report R]` — apply the selected operations in
  order, write the curated GFF and a report of the `Change`s applied.

The detect JSON is stable and parseable so an agent can select candidates programmatically.

## Reporting

- `apply` reuses `normalize/report.py:Change` (action / target / message).
- `detect` introduces `Candidate` (`repair/report.py`); both render to human text and JSON.

## Testing

- **Per operation:** `detect` finds exactly the right candidates; `apply` writes the right
  attributes / retype + note; **re-`detect` after `apply` returns zero** (idempotence).
- **Fixtures:** minimal GFF + FASTA triples for (a) an internal-stop CDS, (b) a UTR-absent mRNA,
  (c) a CDS missing start and/or stop codon; include strand `+`/`-` variants for the partial
  coordinate mapping.
- **Round-trip:** the curated output GFF re-parses and passes `validate`.
- **Integration:** `detect → apply → write → re-parse` through the CLI.

## Non-goals / out of scope

- Switching `gff2mss` to consume the explicit `partial` attributes (separate follow-up task).
- Unifying the duplicated `translate.py` between `ddbj_gff` and `gff2mss` (separate follow-up).
- New operations beyond the initial three (the registry makes these additive).
- `pseudogene` handling and whole-locus `misc_feature` replacement (explicitly not chosen).

## Extension point

Adding an operation = write its `detect` and `apply` functions in `operations.py` and register an
`Operation` entry in `REGISTRY`. The contract is documented in `repair/__init__.py`.
