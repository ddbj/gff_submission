import pytest

from ddbj_gff.errors import GffParseError, Severity
from ddbj_gff.parser import parse

CANONICAL = """\
##gff-version 3
##sequence-region chr1 1 100
chr1\tS\tgene\t1\t99\t.\t+\t.\tID=g1;Name=x
chr1\tS\tmRNA\t1\t99\t.\t+\t.\tID=m1;Parent=g1
chr1\tS\tCDS\t1\t10\t.\t+\t0\tID=c1;Parent=m1;part=1
chr1\tS\tCDS\t50\t99\t.\t+\t2\tID=c1;Parent=m1;part=2
"""


def test_directives_and_features_parsed():
    doc = parse(CANONICAL)
    assert doc.gff_version == "3"
    assert doc.sequence_regions == {"chr1": (1, 100)}
    assert {f.id for f in doc.features} == {"g1", "m1", "c1"}


def test_same_id_rows_aggregate_into_one_feature_with_spans():
    doc = parse(CANONICAL)
    cds = doc.get("c1")
    assert len(cds.spans) == 2
    assert [s.part for s in cds.ordered_spans()] == [1, 2]


def test_parent_child_graph_resolved():
    doc = parse(CANONICAL)
    g1 = doc.get("g1")
    assert [c.id for c in g1.children] == ["m1"]
    assert doc.get("m1").parents[0].id == "g1"
    assert [r.id for r in doc.roots] == ["g1"]


def test_forward_parent_reference_resolves():
    text = (
        "chr1\tS\tmRNA\t1\t99\t.\t+\t.\tID=m1;Parent=g1\n"
        "chr1\tS\tgene\t1\t99\t.\t+\t.\tID=g1\n"
    )
    doc = parse(text)
    assert doc.get("g1").children[0].id == "m1"


def test_dangling_parent_records_warning():
    doc = parse("chr1\tS\tmRNA\t1\t9\t.\t+\t.\tID=m1;Parent=ghost\n")
    assert any(d.code == "dangling-parent" and d.severity == Severity.WARNING
               for d in doc.diagnostics)


def test_no_id_rows_become_standalone_features():
    doc = parse("chr1\tS\tregion\t1\t9\t.\t+\t.\tNote=x\n")
    assert len(doc.features) == 1
    assert doc.features[0].id is None


def test_peptide_fasta_loaded():
    text = CANONICAL + "##FASTA\n>c1\nMAAA\n"
    doc = parse(text)
    assert str(doc.fasta["c1"]) == "MAAA"


def test_strict_raises_on_error_diagnostic():
    with pytest.raises(GffParseError):
        parse("too\tfew\tcols\n", strict=True)
