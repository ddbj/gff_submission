# Trans-splicing support — design

**Date:** 2026-07-05
**Status:** design (for review)
**Scope:** two repos — `gff_submission` (ddbj-gff canonicalization) + `ddbj_mss_tools` (gff2mss conversion)

## Problem

A **trans-spliced** feature is assembled from segments that lie at different genomic
locations (and possibly different strands, or even different sequences). INSDC/DDBJ
represents it as a single feature whose location is a `join(...)` with a `/trans_splicing`
qualifier — the GFF column-4/5 coordinates of any single row cannot express it. In canonical
INSDC GFF3 the true location is carried in a **`location=` attribute** (spec v0.5; gap review
A-3), and the multi-segment feature appears as **multiple GFF rows sharing one ID**, each row
a `part=N` segment (the parser already merges same-ID rows into one multi-span `Feature`,
preserving `part`). `exception=trans-splicing` marks the feature; `is_ordered=true` means the
parts are in biological order and must not be reordered.

**Worked example (real, verified end-to-end):** *Marchantia polymorpha* chloroplast
`AP025455.1` (120,306 bp, circular), gene `rps12` (`locus_tag=Mp_Cg00010`), fetched from
DDBJ getentry. Its CDS is trans-spliced across three segments on mixed strands, and it has
two introns (one trans, one cis):

```
CDS     join(complement(66689..66802),93..324,829..854)   /trans_splicing  /transl_table=11
        -> MPTIQQLIRNKRQPIENRTKSPALKGCPQRRGVCTRVYTTTPKKPNSALRKIARVRLTSGFEITAYIPGIGH
           NLQEHSVVLVRGGRVKDLPGVRYHIIRGTLDAVGVKDRQQGRSKYGVKKSK   (124 codons incl. stop)
intron  join(complement(65903..66688),1..92)              /trans_splicing  /number=1
intron  325..828                                          /number=2   (cis; no trans_splicing)
```

The CDS location string and the protein were **confirmed** by extracting the location from
the real `AP025455.1` sequence with BioPython and translating with table 11 (M-start,
0 internal stops, matches the record's `/translation` exactly).

### Current behaviour (gaps)

- **ddbj-gff:** parser merges the multi-part rows and `model.is_trans_spliced` /
  `is_ordered` / `ordered_spans()` (part-sorted) all work. `validate.rule_special_case`
  emits `noncanonical-special-case` when a trans-spliced feature has **no** `location=`
  attribute. But **no normalize pass builds** the `location=` attribute — so the canonical
  form is never produced from part rows.
- **gff2mss:** never reads `location=`, `is_trans_spliced`, or `part`. `build_insdc_location`
  and `extract_seq` assume a **single strand** for all spans (`spans[0].strand`), so a
  mixed-strand trans-spliced feature gets a wrong location and a wrong translation. It emits
  no `/trans_splicing` qualifier. Introns are treated as `_STRUCTURAL` and **skipped**
  entirely (never emitted as features).

## Decisions (confirmed)

1. **When trans-splicing is present, `normalize` canonicalizes it into a `location=join(...)`
   attribute** (on the ddbj-gff side). This is the boundary contract gff2mss consumes.
2. **Translation is in scope**: a trans-spliced CDS must translate correctly (per-part
   strand), verified against the real rps12 protein.
3. **Intron emission**: the canonical GFF keeps intron features when the input has them, and
   gff2mss emits them. The AUGUSTUS/BRAKER "one intron row per intron" explosion (heterosigma
   nuclear = 101,627) is handled **upstream** by stripping introns during that tool's
   conversion — NOT by a gff2mss filter. heterosigma is out of scope here. (See the
   intron-feature-policy memory.)
4. **Canonical GFF keeps the part rows as-is** (coordinates unchanged); the `location=`
   attribute is the added canonical artifact. If a `location=` attribute already exists on
   input it is **authoritative and preserved** (it may contain remote seqids that col-4/5
   cannot express — gap review A-3).

## Convention (single source of truth)

For a trans-spliced feature with ordered parts `p1, p2, …` (by `part`, honoring
`is_ordered`), build a BioPython `CompoundLocation` of per-part `FeatureLocation`s **in part
order**, each with its own strand:

- GFF strand `-` → `strand=-1` (emitted as a per-part `complement(...)`)
- GFF strand `+`, `?`, `.` → `strand=+1` (no complement)

Then `Bio.SeqIO.InsdcIO._insdc_location_string(compound, seqlen)` yields the INSDC string,
and `compound.extract(genome_seq)` yields the 5′→3′ coding sequence. **Verified** (BioPython,
in-container):

- part order is preserved and per-part `complement` is applied for mixed strands:
  `[−(66689..66802), +(93..324), +(829..854)]` → `join(complement(66689..66802),93..324,829..854)`
- the two-part intron `[−(65903..66688), +(1..92)]` → `join(complement(65903..66688),1..92)`
- extraction + `translate(table=11)` of the CDS compound → the exact rps12 protein.

> **Same-strand note:** if *all* parts share one strand, `_insdc_location_string` applies its
> canonical sort (e.g. all-minus → `complement(join(ascending…))`). That is the standard
> INSDC form for a same-strand join and is accepted (confirmed acceptable by the user). The
> mixed-strand case — the one that matters for trans-splicing — preserves part order as above.

## Part A — ddbj-gff (canonicalization)

### A1. New normalize pass `pass_trans_splicing_location`

**File:** `src/ddbj_gff/normalize/passes.py` (+ register in `normalize.py` `ALL_PASSES`).

For every feature `F` with `F.is_trans_spliced` **and more than one span** **and no existing
`location` attribute**:
1. Order parts via `F.ordered_spans()` (part-sorted).
2. Build a `CompoundLocation` of per-part `FeatureLocation`s using the strand mapping above.
3. Set `F.attributes["location"] = [ _insdc_location_string(compound, seqlen) ]`.
4. Record `Change("add-qualifier", F.id, "built location=join(...) for trans-spliced feature")`.

`seqlen` per seqid from `doc.sequence_regions` (fallback `ctx.seq_lengths`); if unknown, the
`_insdc_location_string` call still works (seqlen only affects origin wrap, not trans-splice).
Features that already have `location=` are left untouched (authoritative). Single-span
trans-spliced features (unusual) are left untouched (col-4/5 suffices).

Registration order: after the structural passes and after `pass_wrap_cds_in_mrna` (so it sees
the merged feature set); before/independent of the biology passes otherwise. Place it
immediately after `pass_circular_origin`.

### A2. validate — no change required

`rule_special_case` already flags trans-spliced-without-`location=`. After A1 the canonical
doc carries `location=`, so the warning disappears. Add a regression test asserting this.

## Part B — gff2mss (conversion)

**File:** `ddbj_mss_tools/src/gff2mss/convert.py`.

### B1. Honor `location=` for the MSS location string

Add a helper used by the CDS/RNA/intron builders: if a feature carries a `location=`
attribute, its MSS location string is that value **verbatim** (this is the A-3 fix and covers
trans-splicing and remote locations uniformly); otherwise fall back to
`build_insdc_location(spans, seqlen, …)` as today.

### B2. Trans-spliced CDS: location + translation + `/trans_splicing`

In `build_cds_feature` (the CDS feature is `cds_feat = next(c for c in mrna.children if
c.type == "CDS")`), when `cds_feat.is_trans_spliced`:
- **Location:** use `cds_feat.attributes["location"][0]` verbatim (guaranteed by A1).
- **Translation (per-part strand):** build a `CompoundLocation` from the CDS spans **ordered
  by `part`**, each `FeatureLocation` with its own strand (mapping above), and
  `.extract(genome_seq)`. `codon_start` from the part-ordered first span's phase. Run the
  existing translation checks (internal stop → `misc_feature` fallback; M-start; multiple-of-3),
  driven by this per-part extraction — NOT `extract_seq`/`_ordered` (single-strand, wrong here).
- **Qualifiers:** the usual `locus_tag / transl_table / codon_start / product / gene /
  inference / submitter-note`, **plus `MssQualifier("trans_splicing", "")`** (a valueless
  qualifier — confirmed: `emit.feature_rows` renders `q.key`/`q.value` directly, so an empty
  value produces a MSS row with an empty qualifier-value column = a bare `/trans_splicing`;
  no emitter change needed).
- **Remote parts out of scope:** if any part references a remote seqid (span seqid ≠ entry
  seqid, or `location=` contains `:`), skip the translation extraction, emit the CDS with the
  verbatim location + `/trans_splicing`, and log a `trans-splicing-remote` WARNING (cannot
  fetch remote sequence to translate). marchantia is all-local and translates.

### B3. Emit intron features

Introns are currently in `_STRUCTURAL` and skipped. Add emission: collect `intron` features
across the document per seqid (they may be children of a trans-spliced gene, a normal gene,
or a `gene_biotype=other` gene — so collect by type, not by walking one gene's children), and
for each emit an `MssFeature("intron", location, quals)`:
- **Location:** via B1 — `location=` verbatim if present (trans intron), else
  `build_insdc_location(spans, seqlen)` (cis intron: single span or same-strand join).
- **Qualifiers** from the intron's own attributes: `gene` (→ `/gene`), `locus_tag`
  (→ `/locus_tag`), `Note` (→ `/note`), `number` (→ `/number`), and `/trans_splicing` when
  `is_trans_spliced`.
- Do **not** re-introduce introns into `build_gene_features`; keep them a separate emission so
  `gene_biotype=other` intron-only genes (rps12 cis intron 2) are covered without the "no-rna"
  skip. Emit in a stable position (e.g. interleaved by start coordinate with the other
  features of the seqid, consistent with existing ordering).

## Verification

Fixture: a minimal, self-contained excerpt of **AP025455.1** preserving the rps12 windows
(positions `1..854` and `65903..66802`), built from the DDBJ-getentry sequence (BioPython
`genbank` parser). Re-coordinatize to keep the fixture small (cp187952 precedent) **or** use a
`1..N` prefix that preserves rps12 coordinates verbatim — the plan chooses; translation must
reproduce the known protein either way. Fixture GFF = the rps12 rows from
`examples/marchantia/chloroplast.gff3` (gene ×2, CDS ×3 parts, intron ×3 rows) adjusted to the
fixture coordinates.

- **normalize (ddbj-gff):** `pass_trans_splicing_location` sets
  `location=join(complement(…),…)` on the trans CDS and the trans intron; leaves the cis
  intron and single-span features untouched; `validate` then emits no `noncanonical-special-case`.
- **gff2mss:** the produced `.ann` contains
  - a `CDS` at `join(complement(66689..66802),93..324,829..854)` with `/trans_splicing`,
    `/transl_table=11`, `/codon_start=1`, `/product`, and a clean translation (a wrong
    translation would instead emit a `misc_feature`);
  - an `intron` at `join(complement(65903..66688),1..92)` with `/trans_splicing` `/number=1`;
  - an `intron` at the cis-intron location with `/number=2` and no `/trans_splicing`.
- **Unit:** location-string building for mixed-strand parts (plus the two verified strings);
  per-part translation of the rps12 CDS → known protein; `?`/`.`→no-complement.
- **Regression:** ddbj-gff (154 not-slow) and gff2mss (154) suites stay green; any existing
  gff2mss fixture that contains intron rows may now emit intron features — update expected
  outputs as needed. heterosigma organelle `.ann` (0 introns, 0 trans) stays byte-identical.

## Non-goals

- Translating trans-spliced CDS whose parts include a **remote seqid** (sequence unavailable
  → diagnostic + verbatim location, no translation).
- The AUGUSTUS/BRAKER intron-stripping preprocessing profile (separate upstream work).
- Rewriting the canonical GFF part rows into a single row (parts are preserved; only
  `location=` is added).
- Frameshift/ribosomal-slippage `exception` handling (separate from trans-splicing).
