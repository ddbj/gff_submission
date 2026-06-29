from ddbj_gff.model import Directive, Feature, GffDocument, Span


def build_doc():
    doc = GffDocument()
    doc.directives = [
        Directive("##gff-version 3", "gff-version", "3"),
        Directive("#!insdc-gff-version 1.0.0", "insdc-gff-version", "1.0.0"),
        Directive("##species ...id=4530", "species", 4530),
        Directive("##sequence-region chr1 1 100", "sequence-region", ("chr1", 1, 100)),
        Directive("#!transl_table primary:1", "transl_table", {"primary": 1}),
    ]
    f = Feature("g1", "DDBJ", "gene", [Span("chr1", 1, 9, "+")], {}, [])
    doc.features = [f]
    doc.feature_index = {"g1": f}
    return doc


def test_directive_accessors():
    doc = build_doc()
    assert doc.gff_version == "3"
    assert doc.insdc_gff_version == "1.0.0"
    assert doc.species == 4530
    assert doc.sequence_regions == {"chr1": (1, 100)}
    assert doc.transl_table_map == {"primary": 1}


def test_get_returns_feature_or_none():
    doc = build_doc()
    assert doc.get("g1").id == "g1"
    assert doc.get("missing") is None


def test_defaults_are_independent():
    a = GffDocument()
    b = GffDocument()
    a.features.append("x")
    assert b.features == []
