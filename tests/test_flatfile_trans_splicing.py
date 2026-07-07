from Bio import SeqIO
from ddbj_gff.model import Span
from ddbj_gff.flatfile.convert import _is_trans, _mark_trans_spliced, qualifiers_to_attrs

FIX = "tests/flatfile_fixtures/trans_splicing_rps12.gbk"


def _cds(rec):
    return next(f for f in rec.features if f.type == "CDS")


def test_is_trans_detects_qualifier():
    rec = SeqIO.read(FIX, "genbank")
    assert _is_trans(_cds(rec)) is True
    src = next(f for f in rec.features if f.type == "source")
    assert _is_trans(src) is False


def test_mark_trans_spliced_sets_exception_ordered_part():
    attrs = {"locus_tag": ["Mp_Cg00010"]}
    spans = [Span("c", 1641, 1754, "-"), Span("c", 93, 324, "+"), Span("c", 829, 854, "+")]
    _mark_trans_spliced(attrs, spans)
    assert attrs["exception"] == ["trans-splicing"]
    assert attrs["is_ordered"] == ["true"]
    assert [s.part for s in spans] == [1, 2, 3]


def test_qualifiers_to_attrs_drops_raw_trans_splicing():
    rec = SeqIO.read(FIX, "genbank")
    attrs = qualifiers_to_attrs(_cds(rec))
    assert "trans_splicing" not in attrs   # replaced by exception (set in synthesize_features)
    assert attrs["product"] == ["ribosomal protein S12"]


from collections import Counter
from ddbj_gff.flatfile.convert import synthesize_features


def test_synthesis_marks_trans_and_emits_introns():
    rec = SeqIO.read(FIX, "genbank")
    feats = synthesize_features(rec, "AP025455")
    counts = Counter(f.type for f in feats)
    assert counts["gene"] == 1 and counts["CDS"] == 1 and counts["intron"] == 2

    cds = next(f for f in feats if f.type == "CDS")
    assert cds.attributes.get("exception") == ["trans-splicing"]
    assert cds.attributes.get("is_ordered") == ["true"]
    assert [s.part for s in cds.spans] == [1, 2, 3]
    assert [s.strand for s in cds.spans] == ["-", "+", "+"]   # per-part strand preserved
    assert "trans_splicing" not in cds.attributes             # raw attr dropped

    gene = next(f for f in feats if f.type == "gene")
    assert len(gene.spans) == 3 and [s.part for s in gene.spans] == [1, 2, 3]   # segment-preserving

    introns = [f for f in feats if f.type == "intron"]
    trans_i = next(f for f in introns if len(f.spans) == 2)
    cis_i = next(f for f in introns if len(f.spans) == 1)
    assert trans_i.attributes.get("exception") == ["trans-splicing"]
    assert trans_i.attributes.get("number") == ["1"]
    assert cis_i.attributes.get("number") == ["2"] and "exception" not in cis_i.attributes
    assert all(i.parent_ids and i.parent_ids[0].startswith("gene-") for i in introns)


def test_trans_gene_spans_have_no_phase():
    rec = SeqIO.read(FIX, "genbank")
    feats = synthesize_features(rec, "AP025455")
    gene = next(f for f in feats if f.type == "gene")
    assert all(s.phase is None for s in gene.spans)   # gene rows must not carry CDS phase


from ddbj_gff.flatfile.convert import flatfile_to_gff
from ddbj_gff.writer import write
from ddbj_gff.validate import validate


def test_flatfile_to_gff_builds_location_and_validates():
    rec = SeqIO.read(FIX, "genbank")
    doc = flatfile_to_gff(rec)
    cds = next(f for f in doc.features if f.type == "CDS")
    intron_trans = next(f for f in doc.features if f.type == "intron" and len(f.spans) == 2)
    # normalize built the canonical location= from the part-ordered per-part-strand spans
    assert cds.attributes["location"] == ["join(complement(1641..1754),93..324,829..854)"]
    assert intron_trans.attributes["location"] == ["join(complement(855..1640),1..92)"]
    # no ERROR-level validate diagnostics (region -> feature-type-not-insdc WARNING is allowed)
    errors = [d for d in validate(doc) if getattr(d, "severity", None) and d.severity.name == "ERROR"]
    assert errors == [], f"unexpected validate errors: {[d.code for d in errors]}"
    # trans-splicing no longer flagged noncanonical (location= present)
    assert "noncanonical-special-case" not in {d.code for d in validate(doc)}
