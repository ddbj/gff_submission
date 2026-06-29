from ddbj_gff.model import Directive, Feature, GffDocument, Span


def doc_with(features, directives=None):
    doc = GffDocument()
    doc.directives = directives or [Directive("##gff-version 3", "gff-version", "3")]
    doc.features = features
    doc.feature_index = {f.id: f for f in features if f.id}
    return doc


def test_equal_ignores_feature_and_attribute_order():
    a = doc_with([
        Feature("g1", "S", "gene", [Span("c", 1, 9, "+")], {"Name": ["x"], "gene": ["x"]}, []),
        Feature("g2", "S", "gene", [Span("c", 20, 29, "+")], {}, []),
    ])
    b = doc_with([
        Feature("g2", "S", "gene", [Span("c", 20, 29, "+")], {}, []),
        Feature("g1", "S", "gene", [Span("c", 1, 9, "+")], {"gene": ["x"], "Name": ["x"]}, []),
    ])
    assert a.semantically_equals(b)


def test_equal_ignores_span_order():
    a = doc_with([Feature("c1", "S", "CDS",
                          [Span("c", 1, 9, "+", 0, part=1), Span("c", 20, 29, "+", 0, part=2)], {}, [])])
    b = doc_with([Feature("c1", "S", "CDS",
                          [Span("c", 20, 29, "+", 0, part=2), Span("c", 1, 9, "+", 0, part=1)], {}, [])])
    assert a.semantically_equals(b)


def test_not_equal_when_span_coords_differ():
    a = doc_with([Feature("g1", "S", "gene", [Span("c", 1, 9, "+")], {}, [])])
    b = doc_with([Feature("g1", "S", "gene", [Span("c", 1, 99, "+")], {}, [])])
    assert not a.semantically_equals(b)


def test_not_equal_when_directive_value_differs():
    a = doc_with([], [Directive("x", "species", 4530)])
    b = doc_with([], [Directive("x", "species", 9606)])
    assert not a.semantically_equals(b)
