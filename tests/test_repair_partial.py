from ddbj_gff.repair.partial import is_partial, partial_attrs
from ddbj_gff.model import Feature


def test_is_partial():
    f = Feature("x", "src", "mRNA", [], {"partial": ["true"]})
    assert is_partial(f) is True
    assert is_partial(Feature("y", "src", "mRNA", [], {})) is False


def test_partial_attrs_plus_5prime():
    # + strand, 5' partial -> start_range on col4
    assert partial_attrs(True, False, "+", 100, 500) == {
        "partial": ["true"], "start_range": [".,100"]}


def test_partial_attrs_plus_3prime():
    assert partial_attrs(False, True, "+", 100, 500) == {
        "partial": ["true"], "end_range": ["500,."]}


def test_partial_attrs_minus_5prime_maps_to_end():
    # - strand: 5' is genomic right (col5) -> end_range
    assert partial_attrs(True, False, "-", 100, 500) == {
        "partial": ["true"], "end_range": ["500,."]}


def test_partial_attrs_both_ends():
    assert partial_attrs(True, True, "+", 100, 500) == {
        "partial": ["true"], "start_range": [".,100"], "end_range": ["500,."]}


def test_partial_attrs_none():
    assert partial_attrs(False, False, "+", 100, 500) == {}
