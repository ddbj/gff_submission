# pass_merge_overlapping_loci — design

**Date:** 2026-07-12
**Status:** design (for review)
**Scope:** new normalize pass in `ddbj-gff` (`gff_submission`) — merge gene loci whose mRNAs
overlap (same strand) into a single gene, so overlapping transcripts share one gene parent.

## Problem / motivation

Some upstream annotations (e.g. AUGUSTUS + evidence-based GMAP models merged into one file) can
produce **separate gene records whose mRNAs overlap the same genomic locus**. For DDBJ they should
be **one locus**: a single `gene` with the overlapping mRNAs as its children. The user wants this
as a general canonicalization rule (robustness / future inputs).

**Current N. benthamiana input is a no-op case:** `Nbe_v1.1.2.sorted.fixed.gff3` was verified to
have 1:1 gene:mRNA and **zero overlapping mRNAs** (the `.fixed` file already resolved them). So this
pass changes nothing there — it is built for future/other inputs and enabled opt-in.

## Decisions (confirmed)

1. **Opt-in, default off** — gated by `NormalizeConfig.merge_overlapping_loci` (bool, default
   `False`). Existing pipelines (heterosigma etc.) are unaffected; the N. benthamiana pipeline sets
   it `True`. This avoids wrongly merging legitimately-nested same-strand genes on curated inputs.
2. **Same strand only** — mRNAs on the same seqid AND same strand can merge. Opposite-strand
   (antisense) overlaps are left as distinct loci.
3. **Percentage overlap threshold (tunable)** — two same-strand mRNAs are "overlapping enough to
   merge" when their overlap fraction ≥ `NormalizeConfig.merge_overlap_min_fraction` (float
   0.0–1.0, **default 0.0** = any ≥1 bp overlap merges). Fraction is defined as
   `overlap_bp / min(len_a, len_b)` (fraction of the shorter mRNA that is covered). The threshold
   is the primary knob the user wants for future tuning; the definition is documented so it can be
   revised.
4. **mRNA-range based** (not CDS-only) — overlap is computed on the mRNA span extents
   (`min(start)`..`max(end)` over each mRNA's spans). (For N. benthamiana, mRNA==CDS extent, so
   equivalent.)
5. **Trans-spliced transcripts are exempt** — a trans-spliced mRNA is **excluded from the overlap
   graph entirely**: it is never merged with another locus, and it never causes others to merge
   (even if its extent overlaps them). A trans-spliced gene (e.g. rps12) spans far-apart segments,
   so its `min..max` extent can be huge and engulf unrelated loci in its gaps; merging on that
   basis would be wrong. A transcript is treated as trans-spliced when the mRNA itself
   `is_trans_spliced` **or any of its CDS children** `is_trans_spliced` (i.e. carries
   `exception=trans-splicing`).

## Behaviour

`pass_merge_overlapping_loci(doc, ctx) -> list[Change]`:

1. If `not ctx.config.merge_overlapping_loci`: return `[]` (no-op).
2. Collect mRNA features grouped by `(seqid, strand)`, **excluding trans-spliced transcripts**
   (an mRNA is skipped when `mrna.is_trans_spliced` or any CDS child `is_trans_spliced`). Each
   remaining mRNA has an extent `(lo, hi)` = `min(s.start)`..`max(s.end)` over its spans, and a
   current gene parent id. Trans-spliced transcripts and their genes are left entirely untouched.
3. Within each `(seqid, strand)` group, build an **overlap graph**: an edge between two mRNAs whose
   extents overlap with fraction ≥ threshold (`overlap/min(len) ≥ frac`). Take **connected
   components** (transitive: A~B, B~C ⇒ {A,B,C}).
4. For each component spanning **≥2 distinct gene parents**:
   - Choose the **representative gene** deterministically = the gene of the member mRNA with the
     smallest `(lo, hi, mrna_id)`.
   - Reassign every component mRNA's `Parent` (attribute + `parent_ids` + `parents`) to the
     representative gene; add them to the representative gene's `children`.
   - Recompute the representative gene's span = union extent `min(lo)`..`max(hi)` of the component's
     mRNAs (single `Span(seqid, lo, hi, strand)`).
   - **Remove** the other (merged-away) gene features from `doc.features` and `doc.feature_index`.
   - Record `Change("merge-loci", rep_gene_id, "merged N genes into locus … (M mRNAs)")`.
   - Components with a single gene parent (incl. existing alt-spliced genes) are untouched.
5. CDS/exon/other children stay under their mRNAs unchanged (only the mRNA→gene link and the gene
   layer change). No dangling parents result (every reassigned mRNA points at the surviving gene).

**Graph maintenance:** after merging, `doc.features` no longer contains merged-away genes,
`doc.feature_index` drops their ids, reassigned mRNAs' `parent_ids`/`attributes["Parent"]`/`parents`
point at the representative, and the representative's `children` includes all component mRNAs. This
keeps `validate.rule_parents` (which checks `feature_index`) clean and downstream passes correct.

**Placement:** register in `ALL_PASSES` **after `pass_wrap_cds_in_mrna`** (so every gene has its
mRNA layer before merging) and before the location-building passes (`pass_circular_origin`,
`pass_trans_splicing_location`).

## Config

Add to `NormalizeConfig`:
- `merge_overlapping_loci: bool = False`
- `merge_overlap_min_fraction: float = 0.0`

Both read from the `[normalize]` table (like existing `wrap_cds_in_mrna`), so the N. benthamiana
normalize invocation can enable them.

## Verification

Synthetic fixture (the real Nbe input has no overlaps): a small GFF with

- two same-strand genes whose mRNAs overlap ≥ threshold → merged into one gene (2 mRNAs, one
  parent), union span, other gene removed;
- a transitive chain A~B~C (A and C don't directly overlap) → all three merged;
- an opposite-strand overlapping gene → NOT merged (stays separate);
- a below-threshold partial overlap (with a raised `merge_overlap_min_fraction`) → NOT merged;
- a **trans-spliced** transcript whose extent overlaps a normal gene → **NOT merged**, and it does
  not pull that normal gene into a merge (both untouched);
- flag off → no change (byte-identical structure).

Assert gene/mRNA counts, parent reassignment, union spans, no dangling parents (`validate` clean),
and that the non-merge cases are untouched. Regression: existing ddbj-gff suite green (flag defaults
off → all existing normalize behaviour unchanged).

## Non-goals

- CDS-level or reciprocal-overlap criteria (single documented definition now; revisable via the
  threshold/definition later).
- Merging across strands (antisense stays separate).
- Choosing a "best" representative by model quality — deterministic lowest-coordinate rule only.
- Splitting (the inverse) — out of scope.
