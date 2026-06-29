from ddbj_gff.model import Directive, Span


def test_span_defaults_and_equality():
    a = Span("chr1", 10, 20)
    b = Span("chr1", 10, 20, strand=".", phase=None, score=None, part=None)
    assert a == b
    assert a.strand == "."


def test_span_sort_key_orders_by_coordinates():
    s1 = Span("chr1", 5, 9, "+", 0)
    s2 = Span("chr1", 10, 20, "+", 1)
    assert s1.sort_key() < s2.sort_key()


def test_span_sort_key_handles_none_phase_and_score():
    # None phase/score must not raise when building the key
    s = Span("chr1", 1, 2)
    key = s.sort_key()
    assert key[0] == "chr1"


def test_directive_holds_kind_and_value():
    d = Directive("##sequence-region chr1 1 100", "sequence-region", ("chr1", 1, 100))
    assert d.kind == "sequence-region"
    assert d.value == ("chr1", 1, 100)
