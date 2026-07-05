from Bio import SeqIO
from ddbj_gff.flatfile.convert import bio_location_to_spans, qualifiers_to_attrs

FIX = "tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk"


def _feat(rec, ftype, lt, seg):
    for f in rec.features:
        if f.type == ftype and f.qualifiers.get("locus_tag", [None])[0] == lt \
                and len(f.location.parts) == seg:
            return f
    raise AssertionError("feature not found")


def test_cds_location_to_spans_1based_phase():
    rec = SeqIO.read(FIX, "genbank")
    cds = _feat(rec, "CDS", "CUMW_191330", 2)          # 1229..2241, 2 segments, + strand
    spans = bio_location_to_spans(cds.location, "S", is_cds=True, codon_start=1)
    assert len(spans) == 2
    assert (spans[0].start, spans[0].strand, spans[0].phase) == (1229, "+", 0)  # codon_start 1 -> phase 0
    assert all(s.strand == "+" for s in spans)


def test_mrna_location_no_phase():
    rec = SeqIO.read(FIX, "genbank")
    mrna = _feat(rec, "mRNA", "CUMW_191340", 3)        # 3 exons
    spans = bio_location_to_spans(mrna.location, "S", is_cds=False)
    assert len(spans) == 3
    assert all(s.phase is None for s in spans)


def test_qualifiers_to_attrs_maps_and_drops():
    rec = SeqIO.read(FIX, "genbank")
    cds = _feat(rec, "CDS", "CUMW_191330", 2)
    attrs = qualifiers_to_attrs(cds)
    assert attrs["locus_tag"] == ["CUMW_191330"]
    assert attrs["product"] == ["hypothetical protein"]
    assert "protein_id" in attrs and "transl_table" in attrs
    assert "translation" not in attrs and "codon_start" not in attrs   # dropped
    assert "Note" in attrs or "note" not in cds.qualifiers             # /note -> Note
