from Bio.SeqFeature import CompoundLocation, FeatureLocation

from ddbj_gff.model import Feature, Span


def test_single_span_location_is_zero_based_halfopen():
    f = Feature("c", "s", "CDS", [Span("chr1", 10, 20, "+", 0)], {}, [])
    loc = f.to_biopython_location()
    assert isinstance(loc, FeatureLocation)
    assert int(loc.start) == 9
    assert int(loc.end) == 20
    assert loc.strand == 1


def test_multi_span_plus_strand_is_compound_in_part_order():
    f = Feature(
        "c", "s", "CDS",
        [Span("chr1", 100, 110, "+", 0, part=2), Span("chr1", 1, 10, "+", 0, part=1)],
        {}, [],
    )
    loc = f.to_biopython_location()
    assert isinstance(loc, CompoundLocation)
    assert [int(p.start) for p in loc.parts] == [0, 99]


def test_minus_strand_location_strand_is_negative():
    f = Feature("c", "s", "CDS", [Span("chr1", 5, 9, "-", 0)], {}, [])
    loc = f.to_biopython_location()
    assert loc.strand == -1


def test_unknown_strand_maps_to_zero():
    f = Feature("c", "s", "CDS", [Span("chr1", 5, 9, "?", 0)], {}, [])
    loc = f.to_biopython_location()
    assert loc.strand == 0
