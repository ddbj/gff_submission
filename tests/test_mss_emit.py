from ddbj_gff.mss.model import MssQualifier, MssFeature, MssEntry, MssDocument
from ddbj_gff.mss.emit import emit_ann, emit_fasta


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
