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
