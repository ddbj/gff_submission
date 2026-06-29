from ddbj_gff.model import Feature, Span


def make_cds(spans):
    return Feature(
        id="cds1",
        source="DDBJ",
        type="CDS",
        spans=spans,
        attributes={"transl_table": ["11"], "Note": ["a", "b"]},
        parent_ids=["gene1"],
    )


def test_feature_property_accessors():
    f = Feature(
        id="g1",
        source="DDBJ",
        type="gene",
        spans=[Span("chr1", 1, 9, "-")],
        attributes={
            "Name": ["rps12"],
            "locus_tag": ["Mp_Cg00010"],
            "gene": ["rps12"],
            "Is_circular": ["true"],
            "is_ordered": ["true"],
            "exception": ["trans-splicing"],
            "number": ["1"],
            "exon_number": ["1"],
            "Dbxref": ["GeneID:1", "X:2"],
        },
        parent_ids=[],
    )
    assert f.name == "rps12"
    assert f.locus_tag == "Mp_Cg00010"
    assert f.gene == "rps12"
    assert f.is_circular is True
    assert f.is_ordered is True
    assert f.is_trans_spliced is True
    assert f.number == 1
    assert f.exon_number == 1
    assert f.dbxref == ["GeneID:1", "X:2"]


def test_trans_spliced_accepts_underscore_or_hyphen():
    f = Feature("x", "s", "CDS", [Span("c", 1, 3, "+", 0)], {"exception": ["trans_splicing"]}, [])
    assert f.is_trans_spliced is True


def test_ordered_spans_plus_strand_by_part():
    spans = [
        Span("c", 100, 110, "+", 0, part=2),
        Span("c", 1, 10, "+", 0, part=1),
    ]
    f = make_cds(spans)
    assert [s.part for s in f.ordered_spans()] == [1, 2]


def test_ordered_spans_minus_strand_descending_when_no_part():
    spans = [
        Span("c", 1, 10, "-", 2),
        Span("c", 100, 110, "-", 0),
    ]
    f = make_cds(spans)
    assert [s.start for s in f.ordered_spans()] == [100, 1]


def test_codon_start_derived_from_first_span_phase():
    spans = [Span("c", 100, 110, "-", 2, part=1), Span("c", 1, 10, "-", 0, part=2)]
    f = make_cds(spans)
    assert f.codon_start == 3  # phase 2 -> codon_start 3


def test_codon_start_none_for_non_cds():
    f = Feature("g", "s", "gene", [Span("c", 1, 9, "+")], {}, [])
    assert f.codon_start is None


def test_ordered_spans_floats_none_part_to_end():
    spans = [Span("c", 1, 10, "+", 0, part=None), Span("c", 100, 110, "+", 0, part=1)]
    f = make_cds(spans)
    ordered = f.ordered_spans()
    assert ordered[0].part == 1
    assert ordered[-1].part is None
