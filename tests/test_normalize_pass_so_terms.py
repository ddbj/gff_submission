from ddbj_gff.model import Feature, Span, GffDocument
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.normalize.passes import NormalizeContext, pass_so_terms
from ddbj_gff.validate.vocab import load_vocab


def _ctx():
    return NormalizeContext(vocab=load_vocab(), seq_lengths=None, config=NormalizeConfig())


def _doc(*feats):
    return GffDocument(features=list(feats))


def test_rename_pseudogenic_cds_adds_pseudo_flag():
    f = Feature("c", "S", "pseudogenic_CDS", [Span("chr1", 1, 9, "+", 0)], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "CDS"
    assert f.attributes.get("pseudo") == ["true"]
    assert any(c.action == "rename-type" for c in changes)
    assert any(c.action == "add-qualifier" for c in changes)


def test_rename_processed_pseudogene_adds_keyed_qualifier():
    f = Feature("g", "S", "processed_pseudogene", [Span("chr1", 1, 9, "+")], {}, [])
    pass_so_terms(_doc(f), _ctx())
    assert f.type == "gene"
    assert f.attributes.get("pseudogene") == ["processed"]


def test_same_name_is_noop():
    f = Feature("c", "S", "CDS", [Span("chr1", 1, 9, "+", 0)], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "CDS"
    assert changes == []


def test_unmapped_type_reported_unchanged():
    f = Feature("x", "S", "totally_made_up_type", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "totally_made_up_type"
    assert any(c.action == "unmapped-type" for c in changes)


def test_placeholder_qualifier_not_fabricated():
    # binding_site -> misc_binding with /bound_moiety="<NAME>" (placeholder, non-duplicate SO term)
    f = Feature("b", "S", "binding_site", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "misc_binding"
    assert "bound_moiety" not in f.attributes   # placeholder value NOT added
    assert any(c.action == "needs-manual" for c in changes)


def test_mobile_genetic_element_renames_and_flags_manual():
    f = Feature("m", "S", "mobile_genetic_element", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "mobile_element"
    assert "mobile_element_type" not in f.attributes      # placeholder value NOT fabricated
    assert any(c.action == "needs-manual" for c in changes)


def test_existing_attribute_not_clobbered():
    f = Feature("c", "S", "pseudogenic_CDS", [Span("chr1", 1, 9, "+", 0)], {"pseudo": ["existing"]}, [])
    pass_so_terms(_doc(f), _ctx())
    assert f.attributes["pseudo"] == ["existing"]
