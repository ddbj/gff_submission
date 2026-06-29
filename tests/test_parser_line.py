from ddbj_gff.errors import Severity
from ddbj_gff.parser import parse_feature_line


def test_parse_basic_cds_line():
    line = "chr1\tDDBJ\tCDS\t10\t20\t.\t+\t0\tID=cds1;Parent=g1;part=1;product=PsbA"
    diags = []
    row = parse_feature_line(line, 5, diags)
    assert row.id == "cds1"
    assert row.type == "CDS"
    assert row.parent_ids == ["g1"]
    assert row.span.start == 10 and row.span.end == 20
    assert row.span.strand == "+" and row.span.phase == 0
    assert row.span.part == 1
    assert "part" not in row.attributes  # part lifted onto the span
    assert row.attributes["product"] == ["PsbA"]
    assert diags == []


def test_wrong_column_count_records_error_and_returns_none():
    diags = []
    row = parse_feature_line("a\tb\tc", 3, diags)
    assert row is None
    assert diags[0].severity == Severity.ERROR
    assert diags[0].code == "col-count"


def test_non_integer_coord_records_error():
    diags = []
    row = parse_feature_line("c\ts\tgene\tx\t20\t.\t+\t.\tID=g", 2, diags)
    assert row is None
    assert diags[0].code == "coord"


def test_start_gt_end_warns_but_parses():
    diags = []
    row = parse_feature_line("c\ts\tCDS\t150\t10\t.\t+\t0\tID=g", 1, diags)
    assert row is not None
    assert any(d.code == "start-gt-end" and d.severity == Severity.WARNING for d in diags)


def test_score_and_phase_dot_become_none():
    diags = []
    row = parse_feature_line("c\ts\texon\t1\t9\t.\t-\t.\tID=e", 1, diags)
    assert row.span.score is None and row.span.phase is None
