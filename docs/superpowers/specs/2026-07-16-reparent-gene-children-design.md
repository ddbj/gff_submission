# pass_reparent_gene_children_to_mrna — design

**Date:** 2026-07-16
**Status:** design (approved, for spec review)
**Scope:** new normalize pass in `ddbj-gff` (`gff_submission`) — reparent a gene's direct
structural sub-features (CDS/exon/intron/UTR/codon) onto the gene's mRNA, so AUGUSTUS-dialect
GFFs where sub-features hang off the gene become the canonical `gene → mRNA → CDS/exon` hierarchy.

## Problem / motivation

Some upstream annotations (AUGUSTUS + GMAP, the N. benthamiana case) emit, for a subset of genes,
sub-features whose `Parent` is the **gene** ID instead of the **mRNA** ID:

```
gene   ID=Nbe.v1.1.s00360g00020
mRNA   ID=Nbe.v1.1.s00360g00020.1 ; Parent=Nbe.v1.1.s00360g00020
CDS    ID=...CDS1 ; Parent=Nbe.v1.1.s00360g00020     <- points at the GENE, not the mRNA
exon   ID=...exon1; Parent=Nbe.v1.1.s00360g00020     <- points at the GENE, not the mRNA
```

The mRNA then has **no CDS/exon children** (it is "empty"), while the CDS/exon are direct children
of the gene, siblings of the mRNA. Downstream `gff2mss` iterates the gene's mRNA child, finds it
empty, hits its `no-exon` path (`convert.py:469`), and **silently skips** the whole gene. In the
N. benthamiana genome this dropped **467 of 84,570 genes** from the `.ann` (verified: 1125 CDS +
1125 exon + 558 intron + 427 start_codon + 435 stop_codon mis-parented under 467 genes; the other
311,622 CDS were parented normally to their mRNA).

The canonical INSDC GFF contract is a strict `gene → mRNA → CDS/exon` hierarchy; making non-canonical
input conform is exactly `normalize`'s job. `pass_wrap_cds_in_mrna` already handles the case where a
gene has CDS/exon and **no** mRNA (it synthesizes an mRNA wrapper). This pass is its **complement**:
the gene has an mRNA **and** gene-level CDS/exon; the fix is to reparent, not to synthesize.

## Decisions (confirmed)

1. **Default on** — gated by `NormalizeConfig.reparent_gene_children` (bool, default `True`),
   mirroring `wrap_cds_in_mrna`. On a well-formed 3-level input the trigger condition is never met,
   so the pass is a **no-op** (structure unchanged) — existing pipelines (heterosigma etc.) are
   unaffected. It only mutates the malformed dialect case.
2. **Reparent, never guess** — the pass acts only when the target mRNA is unambiguous:
   - Gene has **exactly one** transcript child **and it is an `mRNA`**, and the gene has ≥1 direct
     structural child → reparent those structural children onto that mRNA.
   - Gene has **≥2 mRNA** children + gene-level structural children → target ambiguous → **leave
     untouched**, record a `needs-manual` attention diagnostic.
   - Gene's sole transcript is **non-mRNA** (tRNA/rRNA/ncRNA/…) + gene-level structural children →
     **leave untouched**, record `needs-manual` (this is not our case; a tRNA does not take CDS/exon).
   - Gene has **no transcript** at all (only gene-level CDS/exon) → **no-op here**; this is
     `pass_wrap_cds_in_mrna`'s job (which runs immediately after).
3. **Structural sub-feature set** — the direct gene children eligible for reparenting are exactly
   `{CDS, exon, intron, three_prime_UTR, five_prime_UTR, start_codon, stop_codon}`. Transcript-level
   children (`mRNA`, `transcript`, `tRNA`, `rRNA`, `ncRNA`, …) are never reparented.
4. **mRNA span: extend-only** — after reparenting, the mRNA's span grows to cover any reparented
   child that falls outside its current extent (`lo = min(mRNA.lo, children.lo)`,
   `hi = max(mRNA.hi, children.hi)`, single `Span` on the mRNA's seqid/strand). It never shrinks the
   mRNA (a legitimate UTR-bearing mRNA may extend beyond its CDS). For N. benthamiana the mRNA already
   spans its CDS/exon exactly, so this is a no-op there.
5. **Trans-spliced / multi-part mRNA is protected** — if the target mRNA is trans-spliced
   (`mrna.is_trans_spliced`) or already carries a multi-part span structure (`len(mrna.spans) > 1`),
   the pass **skips the gene entirely** (silent `continue`, no change). Such an mRNA carries a
   carefully built per-part span structure — e.g. `flatfile2gff` synthesizes a trans-spliced gene as
   an mRNA with multi-part spans plus `intron` features **intentionally parented to the gene** — and
   the single-`Span` recompute of Decision 4 would collapse it, breaking
   `pass_trans_splicing_location`'s `location=join(...)`. This is not the empty-mRNA dialect bug the
   pass targets, so the gene-level siblings are left as-is. (Discovered during implementation; the
   narrower `is_trans_spliced or len(spans) > 1` guard also lets the pass still reparent stray
   gene-level structural features onto a populated **single-span** mRNA, rather than skipping any
   mRNA that merely has children.)

## Behaviour

`pass_reparent_gene_children_to_mrna(doc, ctx) -> list[Change]`:

1. If `not ctx.config.reparent_gene_children`: return `[]` (no-op).
2. For each `gene` feature in `doc.features`:
   - Partition its `children` into transcripts (`mRNA`/`transcript`/`tRNA`/`rRNA`/`ncRNA`/…) and
     structural children (types in the set from Decision 3).
   - If there are no structural children → skip (nothing to do).
   - If there is **not exactly one transcript that is an `mRNA`**: record a `needs-manual` attention
     (`≥2 mRNA` or `sole transcript non-mRNA`) and skip; if there are **zero** transcripts, skip
     silently (leave for `pass_wrap_cds_in_mrna`).
   - Otherwise let `mrna` be that single mRNA. For each structural child `c`:
     - `c.parent_ids = [mrna.id]`; `c.parents = [mrna]`;
       `c.attributes["Parent"] = [mrna.id]` (only if `Parent` was present).
     - append `c` to `mrna.children`.
   - Remove the reparented children from `gene.children` (gene keeps only its transcript children).
   - Extend `mrna.spans` to a single `Span(seqid, lo, hi, strand)` per Decision 4
     (seqid/strand taken from the mRNA's existing span).
   - Record `Change("reparent-to-mrna", mrna.id, "reparented N gene-level CDS/exon/… of gene G to mRNA M")`.
3. No features are created or removed. `doc.features`, `doc.feature_index`, and `doc.roots` are
   unchanged (only `parent_ids`/`attributes["Parent"]`/`parents`/`children` links and the mRNA span
   are mutated). This keeps `validate.rule_parents` (which checks `feature_index`) clean.

## Config

Add to `NormalizeConfig`:
- `reparent_gene_children: bool = True`

Read from the `[normalize]` table (like `wrap_cds_in_mrna`), so a pipeline can disable it if ever
needed. Default on.

## Placement

Register in `ALL_PASSES` **after `pass_coerce_transcript_to_mrna`** (so `transcript` children have
already become `mRNA` and are recognized) and **before `pass_wrap_cds_in_mrna`** (so this pass fills
existing empty mRNAs first, and `wrap` then synthesizes mRNAs only for genes that still have
gene-level CDS/exon and no transcript). Add `"reparent-to-mrna"` to the `_APPLIED` set in
`normalize.py`.

Resulting order: `pass_directives, pass_coerce_transcript_to_mrna,
pass_reparent_gene_children_to_mrna, pass_wrap_cds_in_mrna, pass_merge_overlapping_loci, …`.

## Verification

Synthetic fixtures (small GFFs) asserting gene/mRNA/CDS counts, parent reassignment, mRNA span, and
`validate` clean:

1. **Empty mRNA + gene-level CDS/exon** → reparented; the mRNA gains the CDS/exon as children; gene
   keeps only the mRNA; `validate` clean; `Change("reparent-to-mrna", …)` recorded.
2. **Well-formed 3-level gene** (CDS/exon already under the mRNA) → **no-op** (structure identical,
   no changes recorded).
3. **Gene with 2 mRNAs + gene-level CDS** → **not** reparented; `needs-manual` attention recorded;
   structure unchanged.
4. **Gene with gene-level CDS and no transcript** → this pass no-op; then `pass_wrap_cds_in_mrna`
   synthesizes the mRNA (integration test running both passes in order).
5. **Gene with a single tRNA + gene-level exon** → **not** reparented (sole transcript non-mRNA);
   `needs-manual` recorded.
6. **Multi-exon mis-parented gene** → all CDS/exon reparented; mRNA span covers them
   (`min(start)..max(end)`).
7. **Flag off** → no change.
8. Regression: existing `ddbj-gff` suite green (flag defaults on → well-formed inputs are no-ops).

## Non-goals

- Choosing among multiple mRNAs (ambiguous case is reported, not resolved).
- Reparenting across genes, or moving transcripts between genes.
- Splitting one mis-parented CDS set across several transcripts.
- Fixing `gff2mss`'s silent-skip behaviour (the diagnostic surfacing is a separate concern; this pass
  removes the *cause* for the common dialect, so the drop no longer occurs).
