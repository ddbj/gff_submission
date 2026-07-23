# transl_except round-trip — design

**Date:** 2026-07-07
**Status:** design (for review)
**Scope:** `ddbj_mss_tools` `gff2mss` — emit the `/transl_except` qualifier so a CDS carrying a
translational exception (e.g. Pyrrolysine read-through) round-trips flatfile ⇄ canonical GFF.

## Problem

For a CDS with `/transl_except=(pos:746..748,aa:Pyl)` (LC757512, Jasmine virus H p87 — an internal
UAG read as Pyl), the round-trip flatfile → GFF (flatfile2gff+normalize) → gff2mss → mss2ff is
**almost** complete:
- flatfile2gff keeps the `transl_except` attribute; `normalize.pass_transl_except` converts it to
  the canonical `recoded_codon` child (`codon_redefined=pyrrolysine`) — **works**;
- gff2mss's `_collect_transl_excepts` reconstructs `(pos:746..748,aa:Pyl)` from that child and
  `translate_cds_with_transl_except` uses it, so the CDS translates correctly (Pyl→`O`, no internal
  stop, emitted as a real CDS not misc_feature) — **works**;
- BUT `build_cds_feature` uses `excepts` **only for translation** and never adds a
  `/transl_except` qualifier to the CDS — so the `.ann` (and the mss2ff `.ff`) **lose the
  qualifier**. The regenerated protein has the `O`, but the feature no longer documents why.

Verified: `_collect_transl_excepts` on the flatfile2gff output reproduces `(pos:746..748,aa:Pyl)`
exactly (aa_names round-trips pyrrolysine→Pyl), so emitting it will match the original.

## Decision

Emit `/transl_except` from the already-collected `excepts` in gff2mss `build_cds_feature`. One
small change; the collection and translation logic already exist.

- **In scope (full round-trip, tested):** a non-trans-spliced CDS with `transl_except`
  (the LC757512 p87 case) — qualifier + correct translation.
- **Out of scope (deferred):** a CDS that is **both trans-spliced and has transl_except**.
  `_build_trans_spliced_cds` translates via `_trans_compound`+`Seq.translate` (no transl_except),
  so a stop-read-through would falsely internal-stop → misc_feature; emitting the qualifier alone
  would be inconsistent. Correct trans+transl_except translation needs a transl_except-aware trans
  path and there is no test data. When the combination is detected, emit a
  `transl-except-trans-splicing` WARNING and do NOT claim the qualifier is honored.

## Approach

In `ddbj_mss_tools/src/gff2mss/convert.py`, `build_cds_feature`: after the CDS `quals` list is
built (the clean-translation path, before returning the `CDS` MssFeature), append one
`MssQualifier("transl_except", spec)` per spec in the already-computed `excepts`:
```python
    for spec in excepts:
        quals.append(MssQualifier("transl_except", spec))
```
No change to the misc_feature (internal-stop) path — a feature that failed translation should not
advertise a transl_except. In `_build_trans_spliced_cds`, if `_collect_transl_excepts(cds_feat)` is
non-empty, emit the `transl-except-trans-splicing` WARNING (deferred combo) and otherwise leave that
path unchanged.

`emit.feature_rows` already renders `MssQualifier` verbatim, and `mss2ff` already parses/writes
`/transl_except` (it round-trips other CDS qualifiers), so no emitter/mss2ff change is needed.

## Test data

`tests/flatfile_fixtures/transl_except_p87.gbk` (already built): a 2400 bp excerpt of LC757512
(real sequence 1..2400) with source + one CDS `20..2320` `/transl_except=(pos:746..748,aa:Pyl)`
`/transl_table=1` `/gene=p87`. The naive (no-exception) translation has exactly one internal stop
at the Pyl position (protein index 242), confirming the exception is load-bearing.

## Verification (round-trip)

flatfile → `flatfile_to_gff` (→ recoded_codon child) → `gff2mss` (emit_mrna=false, transl_table=1)
→ `.ann` → `mss2ff` → `.ff`. Assert against the original:
- the `.ann` CDS carries `transl_except` with value `(pos:746..748,aa:Pyl)`;
- the regenerated `.ff` CDS has `/transl_except=(pos:746..748,aa:Pyl)` and a `/translation`
  containing `O` (Pyl) at the read-through position, and is a `CDS` (not misc_feature);
- feature-type parity for the biological feature (CDS).
Regression: gff2mss + ddbj-gff suites stay green; the marchantia rps12 (transl_except-free)
trans-splicing round-trip and the origin-spanning `transl-except-origin-spanning` diagnostic are
unaffected.

## Non-goals
- Correct translation for a CDS that is simultaneously trans-spliced and transl_except (diagnostic
  only).
- Multi-record / other transl_except aa types beyond what aa_names already supports (Sec, Pyl,
  Term/stop are handled; the mechanism is aa-agnostic via aa_names).
