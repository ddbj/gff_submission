from ddbj_gff.mss.model import MssQualifier, MssFeature, MssEntry, MssDocument
from ddbj_gff.mss.emit import emit_ann, emit_fasta, feature_rows


def test_emit_ann_column_rules():
    doc = MssDocument(
        common_rows=["COMMON\tDBLINK\t\tproject\tPRJDB1"],
        entries=[MssEntry("chr1", [
            MssFeature("source", "1..100", [MssQualifier("organism", "Foo bar"),
                                            MssQualifier("mol_type", "genomic DNA")]),
            MssFeature("CDS", "join(1..10,20..30)", [MssQualifier("locus_tag", "P_000010"),
                                                     MssQualifier("codon_start", "1")]),
        ])],
    )
    out = emit_ann(doc)
    lines = out.splitlines()
    assert lines[0] == "COMMON\tDBLINK\t\tproject\tPRJDB1"
    # entry name only on the entry's first row; feature/location only on feature's first row
    assert lines[1] == "chr1\tsource\t1..100\torganism\tFoo bar"
    assert lines[2] == "\t\t\tmol_type\tgenomic DNA"
    assert lines[3] == "\tCDS\tjoin(1..10,20..30)\tlocus_tag\tP_000010"
    assert lines[4] == "\t\t\tcodon_start\t1"


def test_emit_fasta_wrap_and_terminator():
    seqs = {"chr1": "A" * 130}
    out = emit_fasta(seqs)
    lines = out.splitlines()
    assert lines[0] == ">chr1"
    assert lines[1] == "A" * 60
    assert lines[2] == "A" * 60
    assert lines[3] == "A" * 10
    assert lines[4] == "//"


def test_emit_ann_entry_name_resets_per_entry():
    doc = MssDocument(
        common_rows=[],
        entries=[
            MssEntry("chr1", [MssFeature("source", "1..10", [MssQualifier("organism", "X")])]),
            MssEntry("chr2", [MssFeature("source", "1..20", [MssQualifier("organism", "X")])]),
        ],
    )
    lines = emit_ann(doc).splitlines()
    assert lines[0] == "chr1\tsource\t1..10\torganism\tX"
    assert lines[1] == "chr2\tsource\t1..20\torganism\tX"  # entry name on chr2's first row too


def test_emit_fasta_multi_entry():
    out = emit_fasta({"chr1": "AAA", "chr2": "CCC"})
    assert out.splitlines() == [">chr1", "AAA", "//", ">chr2", "CCC", "//"]


def test_feature_rows_first_row_carries_key_and_location():
    feat = MssFeature("CDS", "1..9", [MssQualifier("locus_tag", "L_1"),
                                       MssQualifier("product", "x")])
    rows = feature_rows(feat)
    assert rows[0] == ["", "CDS", "1..9", "locus_tag", "L_1"]
    assert rows[1] == ["", "", "", "product", "x"]
