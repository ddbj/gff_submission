from Bio.Seq import Seq
from ddbj_gff.model import Feature, Span
from ddbj_gff.mss.convert import collect_spans, build_insdc_location, extract_seq


def test_collect_spans_unions_children_of_type():
    mrna = Feature("m", "S", "mRNA", [Span("c", 1, 100, "+")], {}, [])
    e1 = Feature("e1", "S", "exon", [Span("c", 1, 10, "+")], {}, [])
    e2 = Feature("e2", "S", "exon", [Span("c", 20, 30, "+")], {}, [])
    c1 = Feature("cds", "S", "CDS", [Span("c", 3, 10, "+"), Span("c", 20, 28, "+")], {}, [])
    mrna.children = [e1, e2, c1]
    assert len(collect_spans(mrna, "exon")) == 2
    assert len(collect_spans(mrna, "CDS")) == 2  # NCBI-style single feature with 2 spans


def test_plus_strand_join_and_partials():
    spans = [Span("c", 20, 30, "+"), Span("c", 1, 10, "+")]  # unsorted input
    assert build_insdc_location(spans, 10000) == "join(1..10,20..30)"
    assert build_insdc_location([Span("c", 1, 10, "+")], 10000, five_prime_partial=True) == "<1..10"
    assert build_insdc_location([Span("c", 1, 10, "+")], 10000, three_prime_partial=True) == "1..>10"
    assert build_insdc_location(spans, 10000, five_prime_partial=True, three_prime_partial=True) == "join(<1..10,20..>30)"


def test_minus_strand_complement_ascending_inside():
    spans = [Span("c", 1, 10, "-"), Span("c", 20, 30, "-")]
    assert build_insdc_location(spans, 10000) == "complement(join(1..10,20..30))"
    assert build_insdc_location([Span("c", 5, 9, "-")], 10000) == "complement(5..9)"
    # 5' end of a minus feature is the high-coordinate end -> AfterPosition
    assert build_insdc_location([Span("c", 5, 9, "-")], 10000, five_prime_partial=True) == "complement(5..>9)"
    assert build_insdc_location([Span("c", 5, 9, "-")], 10000, three_prime_partial=True) == "complement(<5..9)"


def test_extract_plus_and_minus():
    genome = Seq("ATGAAATAA")            # +: ATGAAATAA
    assert str(extract_seq([Span("c", 1, 9, "+")], genome)) == "ATGAAATAA"
    rc = Seq("TTATTTCAT")                # revcomp of ATGAAATAA
    assert str(extract_seq([Span("c", 1, 9, "-")], rc)) == "ATGAAATAA"
