# cp187952_origin — circular origin-spanning test fixture

Minimal, self-contained test data for the **circular / origin-spanning** case
(deferred Phase 3-B-full ②), extracted from a real INSDC record.

## Source

- Accession: **CP187952.1** — *Aliinostoc maniaoense* strain TIOX60 chromosome,
  complete genome. **circular**, real length **6,707,124 bp**, taxid 3416947.
- Annotation fetched as NCBI `annotwriter` GFF3 (sviewer `report=gff3`);
  sequence windows fetched via NCBI efetch (`rettype=fasta`, `seq_start`/`seq_stop`).

## The origin-spanning feature

Gene/CDS **ACPZ3T_00005** (`modA`, molybdate ABC transporter substrate-binding
protein, `protein_id=YFD46256.1`, `transl_table=11`) is a minus-strand CDS that
wraps the origin. In the real record it is `complement(join(6706446..6707124, 1..143))`
(the `P_head` piece touching the origin comes first inside the join, so the arc
…6707124 | 1… stays contiguous), which NCBI GFF3 writes as a **single row with
`end > seqlen`**:

```
CP187952.1  Genbank  gene  6706446  6707267  .  -  .  ...locus_tag=ACPZ3T_00005
CP187952.1  Protein Homology  CDS  6706446  6707267  .  -  0  ...
```

`6707267 = 6707124 (seqlen) + 143` — the SO/INSDC end-beyond-length convention.
`Is_circular=true` sits on the `region` landmark feature, **not** on the CDS.

## What was extracted (re-coordinatized)

To keep the fixture small **and** genuinely origin-spanning, the genome was reduced
to a 5,125 bp circular molecule that preserves the real origin:

```
small genome = real[1..2000]  ++  real[6704000..6707124]
             (start window)        (end window)
small length = 5125 ;  small origin (5125|1) == real origin (6707124|1)
```

- start window is identity (`small = real`); end-window mapping is `small = real - 6701999`.
- A fake junction exists at small `2000|2001` (real `2000|6704000`); no kept feature crosses it.
- Sequence bases are the real bases (so the modA CDS still translates correctly).

Re-based features (real → small):

| locus_tag | product | real | small | strand |
|---|---|---|---|---|
| ACPZ3T_00010 | TOBE domain-containing protein | 562..771 | 562..771 | + |
| ACPZ3T_29370 | molybdate ABC transporter substrate-binding protein | 6704448..6705371 | 2449..3372 | + |
| ACPZ3T_29375 | PEP-CTERM sorting domain-containing protein | 6705466..6706278 | 3467..4279 | + |
| **ACPZ3T_00005 (modA)** | molybdate ABC transporter substrate-binding protein | 6706446..6707267 | **4447..5268** | − |

The origin-spanning modA CDS in the fixture is `complement(join(4447..5125, 1..143))`,
written as `4447..5268` with `end 5268 > seqlen 5125` (822 bp = 274 codons). Verified:
this order translates to `MKL…*` (274 codons, 0 internal stops); the reversed
`complement(join(1..143,4447..5125))` gives 24 internal stops.

## Files

- `cp187952_origin.gff3` — NCBI-style GFF3 (5 features: region + 4 gene/CDS pairs).
- `cp187952_origin.fasta` — 5,125 bp (real sequence, re-coordinatized).

## Current behaviour (baseline, before the circular pass)

`validate(parse(cp187952_origin.gff3))` reports `feature-outside-region` for the
gene+CDS of ACPZ3T_00005 (end 5268 > seqlen 5125, and `is_circular` is on the
region, not the feature). Phase 3-B-full ②'s circular pass should recognize the
origin-spanning feature on a circular landmark and mark it `is_circular=true`
(and/or normalize the wrapped coordinates) so the validator no longer flags it.

## Regeneration

`scratchpad/build_cp_fixture.py` (from the fetched windows `cp_start.fa` = real
1..2000 and `cp_endwin.fa` = real 6704000..6707124) rebuilds these files.
