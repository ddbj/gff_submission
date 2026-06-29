from ddbj_gff.model import Directive, Feature, GffDocument, Span
from ddbj_gff.writer import write


def test_write_one_line_per_span_and_part_injected():
    f = Feature(
        "c1", "S", "CDS",
        [Span("chr1", 1, 10, "+", 0, part=1), Span("chr1", 50, 99, "+", 2, part=2)],
        {"ID": ["c1"], "Parent": ["m1"]},
        ["m1"],
    )
    doc = GffDocument(directives=[Directive("##gff-version 3", "gff-version", "3")], features=[f])
    text = write(doc)
    rows = [l for l in text.splitlines() if not l.startswith("#")]
    assert len(rows) == 2
    assert rows[0].split("\t")[3:5] == ["1", "10"]
    assert "part=1" in rows[0]
    assert "part=2" in rows[1]
    assert rows[1].split("\t")[7] == "2"  # phase of span 2


def test_header_directives_emitted_fasta_moved_to_end():
    f = Feature("g1", "S", "gene", [Span("c", 1, 9, "+")], {"ID": ["g1"]}, [])
    doc = GffDocument(
        directives=[
            Directive("##gff-version 3", "gff-version", "3"),
            Directive("##FASTA", "FASTA", None),
        ],
        features=[f],
        fasta={"g1": "MAAA"},
    )
    text = write(doc)
    lines = text.splitlines()
    assert lines[0] == "##gff-version 3"
    assert "##FASTA" in lines
    assert lines.index("##gff-version 3") < lines.index("##FASTA")
    assert lines.index("##FASTA") > lines.index([l for l in lines if l.startswith("c\t")][0])
    assert ">g1" in lines


def test_score_and_phase_dot_when_none():
    f = Feature("e", "S", "exon", [Span("c", 1, 9, "-")], {"ID": ["e"]}, [])
    doc = GffDocument(features=[f])
    row = [l for l in write(doc).splitlines() if l.startswith("c\t")][0]
    cols = row.split("\t")
    assert cols[5] == "." and cols[7] == "."
