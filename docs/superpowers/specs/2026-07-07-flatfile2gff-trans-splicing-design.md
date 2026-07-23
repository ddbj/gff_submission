# flatfile2gff trans-splicing extension — design

**Date:** 2026-07-07
**Status:** design (for review)
**Scope:** extend `flatfile2gff` (`ddbj-gff`) to convert trans-spliced (and intron-bearing) DDBJ
flatfile features into canonical INSDC GFF, so the **flatfile ⇄ GFF round-trip works for
trans-splicing**. Verified on marchantia chloroplast rps12.

## Problem

`flatfile2gff` (nuclear 3-level) does not yet handle trans-spliced features or introns. Running
it on the real trans-spliced input (AP025455 rps12) empirically produces a **broken** canonical
GFF:
- the `/trans_splicing` qualifier maps to a bare `trans_splicing=` attribute — **not** the
  canonical `exception=trans-splicing`; so downstream `is_trans_spliced` is False;
- the CDS spans get no `part=` / `is_ordered`, so `normalize.pass_trans_splicing_location` never
  builds `location=join(...)`; the mixed-strand CDS is emitted as 3 independent-strand rows;
- **intron features are dropped** (intron is not collected by `synthesize_features`);
- the gene spans nearly the whole molecule (min..max of the two far-apart trans segments).

Consequently `gff2mss` does not recognize the CDS as trans-spliced and mis-converts it — the
round-trip fails. (The **reverse** chain `gff2mss → mss2ff` already handles trans-splicing +
introns; only the forward converter is missing this.)

## Decisions (confirmed)

1. **Test data = a small rps12 flatfile fixture** (`tests/flatfile_fixtures/trans_splicing_rps12.gbk`,
   already built): 1754 bp re-coordinatized excerpt of AP025455 (same coordinate system as the
   existing `trans_splicing_rps12.gff3` GFF fixture). Contains source (`/organelle=plastid:chloroplast`,
   `taxon:1480154`), a trans-spliced CDS `join(complement(1641..1754),93..324,829..854)`
   (`/trans_splicing`, transl_table 11, → protein `MPTIQQ…`, 124 codons, 0 internal stops), a trans
   intron `join(complement(855..1640),1..92)` (`/trans_splicing`, `/number=1`), and a cis intron
   `325..828` (`/number=2`). All under `locus_tag=Mp_Cg00010`.
2. **Trans-spliced gene = segment-preserving** (multi-part with `part=`/`is_ordered`), not a
   min..max single span. (Cosmetic for round-trip since gff2mss does not emit gene, but matches the
   canonical chloroplast.gff3 form and avoids a whole-molecule gene.)

## Approach

All changes are in `ddbj-gff`'s `src/ddbj_gff/flatfile/convert.py`; the reverse chain is untouched.

### A. Mark trans-splicing canonically
A flatfile feature carrying a `/trans_splicing` qualifier (BioPython: `"trans_splicing"` in
`feature.qualifiers`) is, in canonical GFF:
- given `exception=trans-splicing` (INSDC canonical marker; drives `Feature.is_trans_spliced`);
- given `is_ordered=true`;
- each of its spans assigned `part=1,2,…` **in BioPython `location.parts` order** (which is the
  biological 5′→3′ order as written in the flatfile join), preserving each part's own strand.

`qualifiers_to_attrs` stops emitting the raw `trans_splicing` attribute (it is replaced by
`exception`). Then the existing `normalize.pass_trans_splicing_location` builds
`location=join(...)` from the part-ordered, per-part-strand spans (verified: for rps12 →
`join(complement(1641..1754),93..324,829..854)`).

This marking applies to the trans-spliced **CDS**, the **synthesized mRNA** for it, the **gene**,
and any trans-spliced **intron**.

### B. Emit intron features (forward path)
`synthesize_features` collects `intron` features (currently dropped) and emits them as children of
their locus_tag's gene (`Parent=gene-<locus_tag>`), carrying `/gene`, `/locus_tag`, `/number`,
`/note`, and — when `/trans_splicing` is present — the canonical trans marking from (A). Introns
are NOT transcripts and get no exon/CDS children. (gff2mss already emits intron features doc-wide,
so they round-trip.)

### C. Trans-spliced gene span
When a locus_tag group is trans-spliced, the synthesized `gene` keeps the child segments as a
multi-part feature (`part=`/`is_ordered=true`) rather than a single min..max span.

### D. Organelle hierarchy (unchanged mechanism)
`flatfile2gff` still emits canonical 3-level (gene→mRNA→exon/CDS). For the organelle round-trip,
`gff2mss` runs with `emit_mrna=false`, collapsing to gene→CDS + introns in the `.ann` — matching
the organelle flatfile (which has no gene/mRNA). `transl_table` comes from the CDS (11 here);
`detect_molecule` already reports `compartment=organelle`, `topology=circular`.

## Verification (round-trip)

Forward: `flatfile_to_gff(rps12.gbk)` → canonical GFF with the trans CDS/intron marked
`exception=trans-splicing` + `location=join(...)`, gene segment-preserving, 2 introns emitted;
passes `ddbj_gff.validate` (no ERROR).
Full loop: GFF → `gff2mss` (emit_mrna=false, transl_table=11) → `.ann` → `mss2ff` → `flatfile'`.
Assert against the original fixture:
- CDS at `join(complement(1641..1754),93..324,829..854)` with `/trans_splicing`, translation
  `MPTIQQ…` (124 codons, 0 internal stop) — matches;
- trans intron at `join(complement(855..1640),1..92)` with `/trans_splicing` `/number=1`;
- cis intron at `325..828` `/number=2`;
- feature-type parity for the biological features (CDS + 2 introns); gene/mRNA/exon are
  intermediate-only (absent from both flatfile endpoints, as for the nuclear case).

Regression: existing flatfile2gff (nuclear Citrus) + ddbj-gff + gff2mss suites stay green.

## Non-goals
- Origin-spanning trans-splicing combined with the circular wrap (rps12 does not cross the origin).
- Reconstructing the separate `gene_biotype=other` intron-holder gene that NCBI annotwriter adds
  for cis introns (our intron is emitted as a child of the main locus_tag gene).
- Multi-record WGS flatfiles (single-record `SeqIO.read`, per the existing follow-up).
