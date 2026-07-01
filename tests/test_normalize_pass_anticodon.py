from ddbj_gff.model import Feature, Span, GffDocument
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.normalize.passes import NormalizeContext, pass_anticodon


def _ctx():
    return NormalizeContext(vocab=None, seq_lengths=None, config=NormalizeConfig())


def _run(f):
    doc = GffDocument(features=[f], feature_index={f.id: f})
    changes = pass_anticodon(doc, _ctx())
    kids = [x for x in doc.features if x.type == "anticodon"]
    return doc, changes, kids


def test_anticodon_becomes_child():
    f = Feature("t", "S", "tRNA", [Span("chr1", 14674, 14742, "-")],
                {"anticodon": ["(pos:complement(14710..14712),aa:Glu,seq:ttc)"]}, [])
    doc, changes, kids = _run(f)
    assert len(kids) == 1
    k = kids[0]
    assert k.type == "anticodon"
    assert k.spans[0].start == 14710 and k.spans[0].end == 14712 and k.spans[0].strand == "-"
    assert k.attributes["amino_acid"] == ["glutamic acid"]
    assert k.attributes["sequence"] == ["ttc"]
    assert k.attributes["Parent"] == ["t"] and k.parent_ids == ["t"]
    assert "anticodon" not in f.attributes
    assert any(c.action == "add-child-feature" for c in changes)


def test_no_anticodon_is_noop():
    f = Feature("t", "S", "tRNA", [Span("chr1", 1, 70, "+")], {}, [])
    doc, changes, kids = _run(f)
    assert changes == [] and kids == []
