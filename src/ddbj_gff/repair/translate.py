# Vendored verbatim from https://github.com/nigyta/translate_with_exception
#   commit d3c382242f1372afb2b49b47a245ba8dcf548cf4, file translate_with_transl_except.py
# NIG (National Institute of Genetics) own code, reused with authorization.
# Public API used here: translate_cds_with_transl_except(feature, parent_seq, stop_symbol="*") -> Seq

#!/usr/bin/env python3
"""
translate_with_transl_except.py

Translate CDS SeqFeature objects that carry transl_except qualifiers,
using Bio.Seq.translate as the core translation engine and applying
transl_except overrides as a post-processing step.

Supported transl_except forms:
  (pos:START..END,aa:TERM)            - partial terminal stop codon (+)
  (pos:POS,aa:TERM)                   - single-base terminal stop codon (+)
  (pos:complement(START..END),aa:Sec) - internal special codon on (-) strand
  (pos:START..END,aa:OTHER)           - non-standard amino acid → X
  (pos:START..END,aa:Met)             - non-AUG alternative start codon → M

Note on region-fetched GenBank files
--------------------------------------
When an NCBI GenBank file is obtained via efetch with seq_start / seq_stop
parameters (e.g. a sub-region of NC_000001.11), BioPython adjusts feature
coordinates to be relative to the region start, but transl_except positions
remain as absolute chromosome coordinates.  The ACCESSION line of such files
reads e.g. "NC_000001 REGION: 25799000..25816000".

_get_region_offset() detects this and returns the 0-based offset to subtract
from transl_except positions before coordinate mapping.
"""

import re
import sys
import warnings

from Bio.Data import CodonTable, IUPACData
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

# 3-letter (and special keyword) → 1-letter amino acid mapping
_AA_3TO1: dict[str, str] = {**IUPACData.protein_letters_3to1_extended}
_AA_3TO1.update({
    "Term":  "*",   # stop codon (partial or otherwise)
    "TERM":  "*",
    "Other": "X",   # non-standard / unspecified amino acid
    "OTHER": "X"
})


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_codon_table(table_id: int | str) -> CodonTable.CodonTable:
    """Return a CodonTable given an NCBI integer ID or name string."""
    try:
        return CodonTable.ambiguous_generic_by_id[int(table_id)]
    except (ValueError, TypeError):
        return CodonTable.ambiguous_generic_by_name[str(table_id)]


def _get_region_offset(parent: Seq | SeqRecord) -> int:
    """
    Return the 0-based offset to subtract from absolute transl_except positions
    so they become relative to the parent sequence passed to feature.extract().

    When NCBI efetch returns a sub-region, BioPython stores the region info in
    record.annotations['accessions'] as e.g.
        ['NC_000001', 'REGION:', '25799000..25816000']
    where 25799000 is the 1-based start of the fetched region on the chromosome.

    Returns 0 for full-sequence (non-region) records or non-SeqRecord inputs.
    """
    if not isinstance(parent, SeqRecord):
        return 0
    accessions = parent.annotations.get("accessions", [])
    try:
        region_idx = accessions.index("REGION:")
        region_str = accessions[region_idx + 1]        # e.g. '25799000..25816000'
        start_1based = int(region_str.split("..")[0])
        return start_1based - 1                        # convert to 0-based offset
    except (ValueError, IndexError, AttributeError):
        return 0


def _parent_pos_to_extracted_pos(parent_pos_0based: int, feature: SeqFeature) -> int:
    """
    Convert a 0-based parent-sequence position to the corresponding 0-based
    position in the sequence returned by feature.extract().

    Works for both SimpleLocation and CompoundLocation (join/order).
    For CompoundLocation the parts are traversed in the order that extract()
    concatenates them (i.e. the order stored in location.parts).
    """

    loc = feature.location

    if isinstance(loc, SimpleLocation):
        if loc.strand == -1:
            return int(loc.end) - 1 - parent_pos_0based
        else:
            return parent_pos_0based - int(loc.start)

    if isinstance(loc, CompoundLocation):
        cumulative = 0
        for part in loc.parts:
            p_start = int(part.start)
            p_end   = int(part.end)
            length  = p_end - p_start
            if part.strand == -1:
                if p_start <= parent_pos_0based < p_end:
                    return cumulative + (p_end - 1 - parent_pos_0based)
            else:
                if p_start <= parent_pos_0based < p_end:
                    return cumulative + (parent_pos_0based - p_start)
            cumulative += length
        raise ValueError(
            f"Parent position {parent_pos_0based} is not within any part of "
            f"feature location {loc}"
        )

    raise TypeError(f"Unsupported location type: {type(loc)}")


def _parse_transl_except(
    qualifier_value: str,
    feature: SeqFeature,
    start_offset: int,
    region_offset: int = 0,
) -> tuple[int, str]:
    """
    Parse one transl_except qualifier value and return (codon_index, one_letter_aa).

    codon_index is 0-based in the protein sequence (after start_offset is applied).

    Args:
        qualifier_value: Raw qualifier string, e.g. '(pos:4263..4264,aa:TERM)'.
        feature:         The parent SeqFeature (used for location mapping).
        start_offset:    codon_start - 1 (0-based frame offset).
        region_offset:   0-based offset to subtract from absolute transl_except
                         positions when the GenBank file was fetched as a
                         sub-region (see _get_region_offset).
    """
    m = re.fullmatch(
        r"\(pos:(complement\()?(\d+)(?:\.\.(\d+))?(?(1)\)),aa:(\w+)\)",
        qualifier_value.strip(),
    )
    if not m:
        raise ValueError(f"Cannot parse transl_except qualifier: {qualifier_value!r}")

    is_complement = bool(m.group(1))
    pos1 = int(m.group(2))          # 1-based, absolute (chromosome) coordinates
    pos2 = int(m.group(3)) if m.group(3) else pos1
    aa_name = m.group(4)

    aa = _AA_3TO1.get(aa_name) or _AA_3TO1.get(aa_name.capitalize())
    if aa is None:
        raise ValueError(f"Unknown amino acid in transl_except: {aa_name!r}")

    # Choose which chromosome position is the "first base of the codon" in the
    # extracted sequence, then convert to 0-based and apply region offset so
    # the result is relative to the parent sequence passed to feature.extract().
    #   (+) strand / no complement → leftmost position (pos1)
    #   (-) strand / complement    → rightmost position (pos2), because the
    #                                extracted (-) sequence is reverse-complemented
    if is_complement:
        parent_pos_0based = (pos2 - 1) - region_offset
    else:
        parent_pos_0based = (pos1 - 1) - region_offset

    pos_in_feat = _parent_pos_to_extracted_pos(parent_pos_0based, feature)
    codon_idx   = (pos_in_feat - start_offset) // 3
    return codon_idx, aa


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def translate_cds_with_transl_except(
    feature: SeqFeature,
    parent_seq: Seq | SeqRecord,
    stop_symbol: str = "*",
) -> Seq:
    """
    Translate a CDS SeqFeature, correctly handling transl_except qualifiers.

    Compared with SeqFeature.translate(cds=True) this function additionally:
      - Handles partial terminal stop codons (1 or 2 genomic bases, aa:TERM)
      - Handles internal special codons such as selenocysteine (aa:Sec → 'U')
      - Handles non-standard amino acids (aa:OTHER → 'X')
      - Handles non-AUG alternative start codons (aa:Met)
      - Handles CompoundLocation (join/order) features
      - Handles region-fetched GenBank files where transl_except positions are
        absolute chromosome coordinates but feature coordinates are region-relative

    The result mimics cds=True behaviour:
      - The initiator codon is always translated as M.
      - The trailing stop symbol is stripped from the returned sequence.

    Args:
        feature:     A SeqFeature of type CDS.
        parent_seq:  The parent DNA Seq or SeqRecord.  Pass the SeqRecord
                     (not just .seq) when working with region-fetched files so
                     that the region offset can be detected automatically.
        stop_symbol: Single character used for in-frame stop codons (default '*').

    Returns:
        Seq: The translated protein sequence (no trailing stop symbol).
    """
    # Detect region offset before stripping the SeqRecord wrapper.
    region_offset = _get_region_offset(parent_seq)
    if isinstance(parent_seq, SeqRecord):
        parent_seq = parent_seq.seq

    # --- codon_start / start_offset ---
    start_offset = int(feature.qualifiers.get("codon_start", [1])[0]) - 1

    # --- extract CDS subsequence ---
    feat_seq: Seq = feature.extract(parent_seq)[start_offset:]

    # Pad to a multiple of 3 so that partial terminal codons (covered by a
    # transl_except with aa:TERM) do not cause a warning from Seq.translate.
    remainder = len(feat_seq) % 3
    if remainder:
        feat_seq = feat_seq + Seq("N" * (3 - remainder))

    # --- codon table ---
    codon_table_id = feature.qualifiers.get("transl_table", ["Standard"])[0]
    codon_table = _get_codon_table(codon_table_id)

    # --- parse transl_except qualifiers ---
    exceptions: dict[int, str] = {}
    for te in feature.qualifiers.get("transl_except", []):
        idx, aa = _parse_transl_except(te, feature, start_offset, region_offset)
        exceptions[idx] = aa

    # --- translate (no CDS validation; we handle that ourselves) ---
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        protein = list(
            str(
                feat_seq.translate(
                    table=codon_table_id,
                    stop_symbol=stop_symbol,
                    to_stop=False,
                    cds=False,
                )
            )
        )

    # --- apply transl_except overrides ---
    for idx, aa in exceptions.items():
        if 0 <= idx < len(protein):
            protein[idx] = aa
        else:
            warnings.warn(
                f"transl_except codon index {idx} is out of range "
                f"(protein length {len(protein)}); qualifier ignored.",
                stacklevel=2,
            )

    # --- mimic cds=True: initiator codon → M ---
    # Only force M for non-ATG alternative start codons (ATG already maps to M).
    first_codon = str(feat_seq[:3]).upper()
    if first_codon in codon_table.start_codons and 0 not in exceptions:
        protein[0] = "M"

    # --- mimic cds=True: strip trailing stop symbol ---
    if protein and protein[-1] == stop_symbol:
        protein.pop()

    return Seq("".join(protein))


# ---------------------------------------------------------------------------
# Command-line demo / test
# ---------------------------------------------------------------------------

def _run_tests() -> bool:
    from Bio import SeqIO

    test_files = [
        ("test_data/GQ332765.gb",          "genbank"),
        ("test_data/CP028126.gb",          "genbank"),
        ("test_data/NC_004353.gb",         "genbank"),
        ("test_data/NC_000001_SELENON.gb",          "genbank"),
        ("test_data/NC_000001_SELENON_relative.gb", "genbank"),
        ("test_data/NC_000913.3.gb",       "genbank"),
        ("test_data/CP000099.1.gb",        "genbank"),
        ("test_data/CP009512.1.gb",        "genbank"),
        ("test_data/NC_004750.gb",         "genbank"),
        ("test_data/U00096.3.gb",          "genbank"),
        ("test_data/AP010953.1.ddbj",      "genbank"),
    ]

    total = passed = 0

    for path, fmt in test_files:
        print(f"\n{'='*60}")
        print(f"File: {path}")
        print("=" * 60)
        try:
            record = SeqIO.read(path, fmt)
        except FileNotFoundError:
            print(f"  [SKIP] file not found: {path}")
            continue

        region_offset = _get_region_offset(record)
        if region_offset:
            print(f"  region offset: {region_offset} (0-based)")

        for feature in record.features:
            if feature.type != "CDS" or "transl_except" not in feature.qualifiers:
                continue

            locus_tag = feature.qualifiers.get("locus_tag", ["?"])[0]
            gene      = feature.qualifiers.get("gene", [""])[0]
            label     = locus_tag if locus_tag != "?" else gene or str(feature.location)

            print(f"\n  CDS: {label}  location={feature.location}")

            try:
                our_prot = str(translate_cds_with_transl_except(feature, record))
            except Exception as exc:
                for te in feature.qualifiers["transl_except"]:
                    print(f"    transl_except: {te}")
                print(f"    [ERROR] {exc}")
                total += 1
                continue

            ref_prot = feature.qualifiers.get("translation", [None])[0]

            # Per-transl_except detail: show qualifier, codon index, and the
            # amino acid at that position in both translations.
            start_offset = int(feature.qualifiers.get("codon_start", [1])[0]) - 1
            for te in feature.qualifiers["transl_except"]:
                try:
                    idx, declared_aa = _parse_transl_except(te, feature, start_offset, region_offset)
                except Exception as exc:
                    print(f"    transl_except: {te}  [parse error: {exc}]")
                    continue
                our_aa = our_prot[idx] if 0 <= idx < len(our_prot) else "?"
                ref_aa = ref_prot[idx] if ref_prot and 0 <= idx < len(ref_prot) else "?"
                print(
                    f"    transl_except: {te}"
                    f"  →  aa_index={idx}"
                    f"  declared={declared_aa!r}"
                    f"  ref={ref_aa!r}"
                    f"  ours={our_aa!r}"
                )

            total += 1
            if ref_prot is not None:
                if our_prot == ref_prot:
                    print(f"    [PASS] translation matches /translation qualifier")
                    passed += 1
                else:
                    print(f"    [FAIL] mismatch")
                    for i, (a, b) in enumerate(zip(our_prot, ref_prot)):
                        if a != b:
                            print(f"           first diff at aa index {i}: ours={a!r}, ref={b!r}")
                            break
                    if len(our_prot) != len(ref_prot):
                        print(f"           length: ours={len(our_prot)}, ref={len(ref_prot)}")
            else:
                print(f"    [INFO] no /translation qualifier to compare against")
                print(f"           result={our_prot[:30]}{'...' if len(our_prot)>30 else ''}")
                passed += 1

    print(f"\n{'='*60}")
    print(f"Results: {passed}/{total} passed")

    # ------------------------------------------------------------------
    # Comprehensive CDS test: translate every CDS in U00096.3.gb and
    # compare against the /translation qualifier.
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print("Comprehensive CDS test: all CDS in test_data/U00096.3.gb")
    print("=" * 60)

    all_path = "test_data/U00096.3.gb"
    try:
        all_record = SeqIO.read(all_path, "genbank")
    except FileNotFoundError:
        print(f"  [SKIP] file not found: {all_path}")
        return passed == total

    all_total = all_passed = all_failed = 0
    failures: list[str] = []

    for feature in all_record.features:
        if feature.type != "CDS":
            continue
        ref_prot = feature.qualifiers.get("translation", [None])[0]
        if ref_prot is None:
            continue

        all_total += 1
        locus_tag = feature.qualifiers.get("locus_tag", ["?"])[0]

        try:
            our_prot = str(translate_cds_with_transl_except(feature, all_record))
        except Exception as exc:
            all_failed += 1
            failures.append(f"  [ERROR] {locus_tag}: {exc}")
            continue

        if our_prot == ref_prot:
            all_passed += 1
        else:
            all_failed += 1
            diff_idx = next(
                (i for i, (a, b) in enumerate(zip(our_prot, ref_prot)) if a != b),
                min(len(our_prot), len(ref_prot)),
            )
            msg = f"  [FAIL]  {locus_tag}: first diff at aa_index={diff_idx}"
            if len(our_prot) != len(ref_prot):
                msg += f" (len ours={len(our_prot)}, ref={len(ref_prot)})"
            failures.append(msg)

    for msg in failures:
        print(msg)

    print(f"\n  CDS with /translation: {all_total}")
    print(f"  Passed: {all_passed}  Failed: {all_failed}")

    all_ok = all_failed == 0
    overall_ok = (passed == total) and all_ok
    print(f"\n{'='*60}")
    print(f"Overall: transl_except tests {passed}/{total}, "
          f"comprehensive CDS test {all_passed}/{all_total}")
    return overall_ok


if __name__ == "__main__":
    ok = _run_tests()
    sys.exit(0 if ok else 1)
