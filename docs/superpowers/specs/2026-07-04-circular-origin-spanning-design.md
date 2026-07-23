# Circular origin-spanning features — design

**Date:** 2026-07-04
**Status:** design (for review)
**Scope:** two repos — `gff_submission` (ddbj-gff canonicalization) + `ddbj_mss_tools` (gff2mss conversion)

## Problem

A feature on a **circular** molecule can cross the origin (e.g. an organelle CDS that
runs off the end of the sequence and continues from position 1). INSDC/NCBI GFF3
represents such a feature as a **single row with `end > seqlen`**, where
`end = true_end + landmark_length`. `Is_circular=true` sits on the **landmark**
(`region`/`source`) feature, not on the crossing feature.

This is the deferred "Phase 3-B-full ②" circular pass, exercised by the
`tests/normalize_fixtures/cp187952_origin.*` fixture (real INSDC record, re-coordinatized
to 5125 bp): CDS `modA` (`ACPZ3T_00005`, minus strand) is written `CP187952.1 … 4447 5268 … -`
with `end 5268 > seqlen 5125`, meaning `complement(join(4447..5125, 1..143))`.

**Not a blocker for the current Heterosigma submission:** `organelle.normalized.gff3`
(MT 39815 / CP 177054) contains **zero** origin-spanning features and no `Is_circular`
in the GFF (circularity there comes from `sequence_roles.tsv` on the gff2mss side).
This work is robustness/completeness, verified by the cp187952 fixture.

### Current behaviour (gaps)

- **ddbj-gff:** `feature.is_circular` reads the feature's own `Is_circular` attribute.
  On the cp187952 fixture the flag is on the `region` landmark, not the CDS, so
  `validate.rule_seqid_bounds` emits `feature-outside-region` for the modA gene+CDS
  (`end 5268 > hi 5125` while `not circular`). No normalize pass propagates circularity.
- **gff2mss:** `build_insdc_location` and `extract_seq` treat `end>seqlen` literally →
  `_insdc_location_string` yields an invalid MSS location (`complement(4447..5268)` on a
  5125 bp entry), and `extract_seq` reads past the sequence end → wrong translation.

## Decisions (confirmed)

1. **Scope = both repos**: canonicalization (ddbj-gff) + conversion (gff2mss), verified
   end-to-end on cp187952.
2. **Canonical GFF representation = keep `end>seqlen` + propagate `Is_circular`**
   (INSDC/NCBI-faithful, round-trippable). Coordinate→`join` splitting is a
   **conversion-side** (gff2mss) concern only. The canonical GFF body is not rewritten
   into two spans.

## Convention (single source of truth)

For a circular molecule of length `L`, an origin-spanning feature is a single span with
`start ≤ L < end`. Its two in-bounds genomic pieces are:

- `P_head = start .. L`
- `P_tail = 1 .. (end − L)`

INSDC location strings (`P_head` — the piece touching `L` — comes first for both strands;
`complement()` marks minus). This follows the INSDC identity
`complement(join(A,B)) = join(complement(B), complement(A))`, so the origin-crossing arc
(…→L | 1→…) stays **contiguous** only when `P_head` precedes `P_tail` inside the join:

- **plus strand:**  `join(start..L, 1..(end−L))`
- **minus strand:** `complement(join(start..L, 1..(end−L)))`

Fixture modA (minus, `4447..5268`, L=5125): **`complement(join(4447..5125,1..143))`**.
Empirically verified (BioPython in the ddbj-gff-dev container): this string, and only this
one, translates to the real protein (`MKL…*`, 274 codons, **0 internal stops**);
the reversed `complement(join(1..143,4447..5125))` translates to garbage (24 internal stops).

> **Ordering note (critical).** Feed BioPython the two pieces as `FeatureLocation`s in
> **biological 5′→3′ order** and `_insdc_location_string`/`CompoundLocation.extract` both
> produce the correct string *and* the correct coding sequence — no manual string building.
> Biological order is **plus → `[P_head, P_tail]`**, **minus → `[P_tail, P_head]`**
> (BioPython reverses parts under `complement(join(…))`). The existing `_ordered()`
> (coordinate sort: ascending for plus, descending for minus) yields the **opposite** order
> in the wrap case for *both* strands — so the wrap path must bypass `_ordered()` and use
> the explicit biological order. Verified strings: minus →
> `complement(join(4447..5125,1..143))`; plus → `join(4447..5125,1..143)`; minus 5′-partial →
> `complement(join(4447..5125,1..>143))`; minus 3′-partial → `complement(join(<4447..5125,1..143))`.

## Part A — ddbj-gff (canonicalization)

### A1. New normalize pass `pass_circular_origin`

**File:** `src/ddbj_gff/normalize/passes.py` (+ register in `normalize.py` `ALL_PASSES`).

Behaviour:
1. Build the set of **circular seqids**: a seqid is circular if some landmark feature on
   it (`type in {"region", "source"}`) has a truthy `Is_circular` attribute.
2. Resolve `seqlen` per seqid from `doc.sequence_regions` (`##sequence-region`), falling
   back to `ctx.seq_lengths` when available.
3. For every feature `F` with a span `s` where `s.seqid` is circular and `s.end > seqlen`:
   set `F.attributes["Is_circular"] = ["true"]` if not already present, and record a
   `Change("add-qualifier", F.id, "propagated Is_circular to origin-spanning feature …")`.
   **Coordinates are left unchanged** (canonical keeps `end>seqlen`).

Place the pass early (right after `pass_directives`) so later passes and `validate` see
the flag. Order among the biology passes is otherwise irrelevant.

**Out of scope:** `start>end`-style wrap encoding (only the `end>seqlen` convention);
features on a circular seqid that are within bounds (untouched).

### A2. `validate.rule_seqid_bounds` — landmark-aware circularity

**File:** `src/ddbj_gff/validate/rules.py`.

Today the rule suppresses `feature-outside-region` only when `f.is_circular` (the feature's
own attribute). Make it also suppress the `end > hi` case when the **seqid's landmark** is
circular, so validation is correct even when run on a not-yet-normalized doc (flag still on
the region). Compute circular seqids once (same rule as A1). Keep flagging `start < lo`
(a coordinate `< 1` is always invalid) and keep `end > hi` flagged for **non-circular**
seqids.

## Part B — gff2mss (conversion)

**File:** `ddbj_mss_tools/src/gff2mss/convert.py`.

### B1. Helper `_wrap_spans(spans, seqlen)`

Expand a span with `end > seqlen` into its two in-bounds pieces **in biological 5′→3′
order** (plus `[head, tail]`, minus `[tail, head]`), preserving `seqid`/`strand`:

```python
def _wrap_spans(spans, seqlen):
    """Split an origin-spanning span (end>seqlen) into its two in-bounds pieces,
    ordered biologically 5'->3' (plus [head,tail], minus [tail,head]).
    Returns (spans, wrapped: bool). Non-wrapping input is returned unchanged."""
    out, wrapped = [], False
    for s in spans:
        if s.end > seqlen:
            wrapped = True
            head = Span(s.seqid, s.start, seqlen, s.strand)      # start..L
            tail = Span(s.seqid, 1, s.end - seqlen, s.strand)    # 1..(end-L)
            out += ([tail, head] if s.strand == "-" else [head, tail])
        else:
            out.append(s)
    return out, wrapped
```

Assumption: the common case is a **single** wrapping span (organelle CDS/gene). A feature
that is both multi-exon and origin-spanning is out of scope; `build_cds_feature` logs a
`multi-exon-origin-spanning` WARNING when it sees `end>seqlen` combined with >1 input span,
then still produces **best-effort** wrapped output (there is no separate non-wrapped
fallback — `_wrap_spans` is applied unconditionally on the wrap path). Likewise, a CDS that
combines `transl_except` with an origin-spanning span is out of scope: the transl_except
translation path is not wrap-aware, so `build_cds_feature` emits a
`transl-except-origin-spanning` WARNING flagging that its protein may be incorrect.

### B2. `build_insdc_location(spans, seqlen, …)`

- If no span wraps: unchanged (existing `_ordered` path).
- If a span wraps: use the biologically-ordered spans from `_wrap_spans` **directly**
  (do NOT call `_ordered` — it inverts the wrap order). Build `FeatureLocation`s in that
  order and reuse the existing partial and `_insdc_location_string` logic — `locs[0]` is the
  biological 5′ end and `locs[-1]` the 3′ end, so the current partial branches apply as-is.

Pinned expectations (tests):
- `build_insdc_location([Span("c",4447,5268,"-")], 5125) == "complement(join(4447..5125,1..143))"`
- `build_insdc_location([Span("c",4447,5268,"+")], 5125) == "join(4447..5125,1..143)"`
- `build_insdc_location([Span("c",4447,5268,"-")], 5125, five_prime_partial=True) == "complement(join(4447..5125,1..>143))"`
- `build_insdc_location([Span("c",4447,5268,"-")], 5125, three_prime_partial=True) == "complement(join(<4447..5125,1..143))"`

### B3. `extract_seq(spans, genome_seq)`

Expand via `_wrap_spans` before building the `CompoundLocation` (same biological order),
so `.extract()` reads the pieces strand-aware and yields the correct 5′→3′ coding
sequence. This feeds translation, internal-stop→misc_feature, and partial detection.
**Verified** in-container: for the fixture modA (`4447..5268 -`, L=5125), the biologically-
ordered extraction translates to `MKL…*` (274 codons, 0 internal stops); the reversed order
gives 24 internal stops. The translation test below is the correctness gate.

Pinned expectation (test): the fixture modA CDS (`4447..5268 -`, 822 bp = 274 codons)
translates starting with `M`, with no internal stop before the terminal codon.

## Verification

- **Unit (gff2mss):** `test_mss_location.py` — minus + plus wrap strings above;
  `extract_seq` wrap translation (modA-derived or synthetic).
- **Normalize (ddbj-gff):** new test — `pass_circular_origin` on cp187952 sets
  `Is_circular` on the modA gene+CDS; other features untouched.
- **Validate (ddbj-gff):** `rule_seqid_bounds` on cp187952 (flag on region only) →
  no `feature-outside-region`.
- **End-to-end:** run cp187952 through parse→normalize→validate (clean) then gff2mss →
  `.ann` contains modA at `complement(join(4447..5125,1..143))`. Optionally confirm with
  `ghcr.io/ddbj/ddbj-validator:0.1.4-beta`.
- **Regression:** existing ddbj-gff (147 not-slow / 8 slow) and gff2mss (80) suites stay green;
  Heterosigma organelle `.ann` remains byte-identical (no origin-spanning present).

## Non-goals

- `start>end`-style origin encoding (only `end>seqlen`).
- Multi-exon features that also cross the origin (diagnostic + non-wrapped fallback).
- Changing how Heterosigma declares circularity (`sequence_roles.tsv`, unaffected).
