# flatfile2gff — design (DDBJ flatfile → INSDC canonical GFF)

**Date:** 2026-07-05
**Status:** design (for review)
**Scope:** new converter in `gff_submission` (`ddbj-gff`); nuclear eukaryote (3-level) first.

## Goal

Convert a **DDBJ flatfile** (INSDC feature table) into the project's **canonical INSDC GFF3**,
so that the pair round-trips at the **biological-feature level**:

```
DDBJ flatfile ──[flatfile2gff (NEW)]──▶ canonical INSDC GFF ──[gff2mss]──▶ .ann ──[mss2ff]──▶ DDBJ flatfile'
```

The reverse chain (`gff2mss` → `mss2ff`) already works end-to-end (verified, incl. join/
complement/trans-splicing and `/translation` regeneration). The missing forward converter is
`flatfile2gff`. Success = the biological features (CDS / mRNA / tRNA / rRNA …) survive the full
loop with matching coordinates and translations, and the generated GFF passes `ddbj_gff.validate`.

## Scope decisions (confirmed)

- **Nuclear eukaryote, 3-level (gene → mRNA → exon/CDS) first.** Example: `BDQV01000200.1`
  (*Citrus unshiu*, WGS scaffold, 387,150 bp, linear, division PLN), fetched via DDBJ getentry
  into `examples/citrus_unshiu/`. Organelle (2-level, trans-splicing and other special rules) is
  a **later extension**, deliberately deferred — it has more special cases.
- **Output is the project's canonical INSDC GFF** (feature-mapping.tsv conformant, directives,
  `ddbj_gff.validate`-clean). **NCBI's annotwriter GFF is a logic reference only, NOT the target
  format** — the INSDC-GFF spec is still being drafted and differs from NCBI (see Evidence).
- **COMMON is ignored** (submitter, references, BioProject/BioSample/SRA accessions) — not needed.
  But molecule-type metadata IS extracted: **organism / taxid / division** (for codon table and
  hierarchy) and **topology (linear/circular)** and the **source feature's chromosome-vs-organelle**
  signal (drives processing).
- **Entry / locus_tag IDs are provisional**; they are only finalized when the record is actually
  loaded into the DB. `mss2ff` cannot reproduce DB-assigned IDs — accepted.
- **Success criterion (primary):** ① round-trip biological-feature match (coordinates +
  translation) via `flatfile2gff → gff2mss → mss2ff`; ② generated GFF passes `ddbj_gff.validate`.
  NCBI-GFF comparison is a sanity reference only; divergences are allowed.
- **CDS↔mRNA pairing under alternative splicing:** by **coordinate containment** (recommended
  default; see Feature Mapping).

## Evidence (real data, DDBJ vs NCBI for the same accession)

Fetched `BDQV01000200.1` from both DDBJ (getentry) and NCBI (annotwriter GFF3):

| | DDBJ flatfile | NCBI annotwriter GFF3 |
|---|---|---|
| features | mRNA 37, CDS 37, assembly_gap 36, source 1 | gene 33, mRNA 37, exon 193, CDS 187 (rows), region 1 |
| gene / exon | **none** | added by NCBI |
| distinct locus_tags | **30** | (NCBI genes: 33) |

Key facts this establishes:
- **The DDBJ flatfile has mRNA + CDS but no gene and no exon** — `flatfile2gff` must synthesize
  gene and exon.
- **NCBI's gene count (33) ≠ DDBJ locus_tag count (30)** — NCBI splits genes on a different
  granularity. We must **not** follow NCBI; we group by **`/locus_tag`** (the DDBJ/INSDC key) → 30
  genes.
- **Alternative splicing is present**: 5 locus_tags carry >1 transcript (e.g. `CUMW_191580` has 4);
  37 mRNA over 30 genes.
- **mRNA and CDS are not adjacently ordered** in the flatfile, and there is **no explicit
  transcript↔CDS link** (mRNA has `/note "Transcript ID: …"`, the CDS does not reference it).
  Pairing must be geometric.

Also (organelle reference, `AP025455.1`): DDBJ has no gene; NCBI adds one gene per feature with
`order(...)` locations and a separate gene for the cis-intron holder — informative for the later
organelle extension, not used now.

## Placement & dependency

- Lives in **`gff_submission` (`ddbj-gff`)**: it produces canonical GFF, so it belongs with the
  canonicalizer (`parser` / `model` / `normalize` / `validate`). New module e.g.
  `src/ddbj_gff/flatfile/` (reader + gene/transcript synthesis + GFF emit via the existing `writer`).
- Imports **BioPython** (already a dependency; `SeqIO.read(path, "genbank")` reads DDBJ flatfiles).
- **One-way dependency preserved**: `ddbj-gff` still does not import `gff2mss`/`common`.

## Input & molecule-type determination

Parse with BioPython genbank parser. From the record derive:
- **taxid** — from source `/db_xref="taxon:NNN"` → `##species …?id=NNN` directive + region Dbxref.
- **organism** — source `/organism`.
- **division** — `rec.annotations["data_file_division"]` (e.g. `PLN`).
- **topology** — `rec.annotations["topology"]` (`linear`/`circular`) → region `Is_circular` when circular.
- **compartment** — source `/organelle` present ⇒ organelle (deferred); absent ⇒ **nuclear /
  chromosome ⇒ 3-level** (this spec's path).
- **transl_table** — primary from the CDS `/transl_table` (each CDS carries it); the molecule type
  provides the default when absent.

Molecule type + division select the **hierarchy**: nuclear eukaryote ⇒ gene → mRNA → exon/CDS
(3-level). (Organelle ⇒ 2-level, later.)

## Feature mapping (flatfile → canonical GFF)

Emit canonical GFF via the existing `ddbj_gff.writer`. Per feature:

1. **source → `region`** landmark feature (`##sequence-region` + a region row): carry
   `mol_type`, `organism`, `Dbxref=taxon:…`, `submitter_seqid`/`chromosome`, `Is_circular` (if
   circular). Submission-only source qualifiers (collection_date, geo_loc_name, cultivar, etc.)
   are carried as attributes where harmless but are **not** round-trip-critical (gff2mss rebuilds
   source from its own config).
2. **gene synthesis** — group mRNA/CDS/RNA features by **`/locus_tag`**; one `gene` per group.
   - `gene` location = the union of the group's child segments (single span when contiguous;
     multi-part with `is_ordered=true` when the children are segmented, mirroring the reference
     canonical form — NOT a min..max span that would swallow large gaps).
   - `gene_biotype` from the child type (`protein_coding` for CDS-bearing, `tRNA`/`rRNA`/…).
   - qualifiers: `gene` (from `/gene`), `locus_tag`, `gene_synonym` (if present). No `/product`,
     no `/translation` on the gene.
   - Deterministic synthetic `ID=gene-<locus_tag>`; children get `Parent=` this ID.
3. **mRNA** — one canonical `mRNA` per flatfile mRNA feature (`Parent`=gene). Keep `/product`,
   `/note` (incl. the transcript-ID note), `/locus_tag`.
4. **CDS↔mRNA pairing (alt-splicing)** — within a locus_tag group, assign each CDS to the mRNA
   whose exon structure **contains** it (every CDS segment lies within an mRNA exon and shares
   splice boundaries). If exactly one mRNA qualifies → pair it. If several qualify → pick the
   best-matching (most shared boundaries) and emit a diagnostic. If none qualify → emit a
   diagnostic and fall back to a synthesized mRNA equal to the CDS location. Single-transcript
   genes are the trivial case.
5. **exon synthesis** — for each mRNA, split its location into one `exon` per segment
   (`Parent`=mRNA). **If a gene/transcript has no mRNA feature**, synthesize an `mRNA` equal to the
   CDS location first, then derive exons from it (per the confirmed rule).
6. **CDS** — canonical `CDS` (`Parent`=mRNA), multi-span preserved from the flatfile join (each
   segment a span with `phase` from `codon_start`: `codon_start=1→phase 0`, `2→1`, `3→2`). Keep
   `/product`, `/protein_id`, `/transl_table`, `/note`, `/db_xref`→`Dbxref`. Drop `/translation`
   (derivable; regenerated by mss2ff on the way back).
7. **tRNA / rRNA / other RNA** — mapped via feature-mapping.tsv; gene synthesized as above.
8. **assembly_gap** — **dropped** (round-trip regenerates gaps from N-runs via gff2mss's gap
   annotator; carrying them as features is redundant).
9. **Qualifier → attribute** mapping follows INSDC qualifier ↔ GFF attribute conventions already
   encoded in the project (Name/gene/locus_tag/product/Note/Dbxref/protein_id/gene_synonym/
   transl_table/codon_start→phase).

Directives emitted: `##gff-version 3`, `#!insdc-gff-version <v>`, `##sequence-region <id> 1 <len>`,
`##species …?id=<taxid>`, `#!transl_table primary:<n>` — via the existing normalize/writer path.

## Round-trip boundary (what does / doesn't survive the loop)

- **Survives (biological, the target):** sequence; mRNA + CDS features; coordinates incl.
  join/complement; `/product`, `/locus_tag`, `/transl_table`, `/codon_start`, `/protein_id`,
  `/gene`; the CDS translation (regenerated by mss2ff, must match).
- **Intermediate-only (in the GFF, absent from both flatfile endpoints):** `gene`, `exon` — the
  DDBJ flatfile has neither, and gff2mss does not emit them to `.ann`, so `flatfile'` also lacks
  them. Their correctness is validated by `ddbj_gff.validate` (and, informatively, by comparison
  to the NCBI GFF), NOT by the round-trip.
- **Side-channel (not in GFF):** COMMON (references/DBLINK/comment/dates) and gff2mss's source
  construction (from config/common.json). Ignored per scope.
- **Regenerated, not carried:** `/translation` (recomputed) and `assembly_gap` (from N-runs).

## Verification

- **Fixture:** `examples/citrus_unshiu/BDQV01000200.ddbj.gbk` (real DDBJ flatfile). A small
  re-coordinatized excerpt (a handful of genes incl. one alt-spliced locus like `CUMW_191340`)
  will be the committed test fixture, to keep it light while exercising multi-exon + alt-splicing.
- **① round-trip:** `flatfile2gff(flatfile)` → canonical GFF; then `gff2mss` → `.ann` → `mss2ff`
  → `flatfile'`. Assert the mRNA and CDS feature sets match between the original flatfile and
  `flatfile'` by location and by regenerated translation (BioPython feature-by-feature compare).
- **② conformance:** the generated canonical GFF passes `ddbj_gff.validate` with no errors.
- **Sanity (informative, non-gating):** compare gene grouping / feature coverage against the NCBI
  annotwriter GFF3; note (do not fail on) divergences (e.g. NCBI's 33 vs our 30 genes).

## Out of scope (this spec)

- Organelle (2-level) and trans-splicing / `order()` special rules (later extension; the reverse
  side already handles trans-splicing).
- COMMON reconstruction and DB-assigned ID finalization.
- Frameshift / non-standard `/transl_except` beyond what the existing translate path covers.

## Risks / open points

- **CDS↔mRNA containment ambiguity** for heavily overlapping alt-transcripts — mitigated by the
  best-boundary-match + diagnostic rule; rare.
- **gene location for spatially separated child segments** — use multi-part/`is_ordered`, not
  min..max, to avoid a gene spanning unrelated regions.
- **`/translation` mismatch** would indicate a coordinate/phase/codon-table error — this is the
  point of the round-trip check.
