from ddbj_gff.model import Feature, Span, GffDocument
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.normalize.passes import NormalizeContext, pass_transl_except


def _ctx():
    return NormalizeContext(vocab=None, seq_lengths=None, config=NormalizeConfig())


def _cds(attrs):
    # CDS spanning 1..582 on chr1 (+)
    return Feature("c", "S", "CDS", [Span("chr1", 1, 582, "+", 0)], dict(attrs), [])


def _run(f):
    doc = GffDocument(features=[f], feature_index={f.id: f})
    changes = pass_transl_except(doc, _ctx())
    kids = [x for x in doc.features if x.type in ("recoded_codon", "stop_codon")]
    return doc, changes, kids


def test_transl_except_becomes_recoded_codon():
    f = _cds({"transl_except": ["(pos:139..141,aa:Sec)"]})
    doc, changes, kids = _run(f)
    assert len(kids) == 1
    k = kids[0]
    assert k.type == "recoded_codon"
    assert k.spans[0].start == 139 and k.spans[0].end == 141
    assert k.attributes["codon_redefined"] == ["selenocysteine"]
    assert k.attributes["Parent"] == ["c"] and k.parent_ids == ["c"]
    assert k.attributes["ID"] == [k.id]
    assert "transl_except" not in f.attributes
    assert any(c.action == "add-child-feature" for c in changes)
    assert k in f.children and doc.feature_index.get(k.id) is k


def test_complement_pos_sets_minus_strand():
    f = _cds({"transl_except": ["(pos:complement(139..141),aa:Sec)"]})
    _, _, kids = _run(f)
    assert kids[0].spans[0].strand == "-"


def test_stop_aa_becomes_stop_codon():
    f = _cds({"transl_except": ["(pos:580..582,aa:Term)"]})
    _, _, kids = _run(f)
    assert kids[0].type == "stop_codon"


def test_out_of_bounds_is_needs_manual_and_kept():
    f = _cds({"transl_except": ["(pos:9000..9002,aa:Sec)"]})
    doc, changes, kids = _run(f)
    assert kids == []
    assert f.attributes.get("transl_except") == ["(pos:9000..9002,aa:Sec)"]
    assert any(c.action == "needs-manual" for c in changes)


def test_no_transl_except_is_noop():
    f = _cds({})
    doc, changes, kids = _run(f)
    assert changes == [] and kids == []
