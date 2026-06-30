from ddbj_gff.validate.vocab import load_vocab


def test_feature_qualifiers_loaded():
    v = load_vocab()
    # pseudogenic_CDS -> CDS, qualifier /pseudo
    assert v.insdc_map["pseudogenic_CDS"] == "CDS"
    assert any("pseudo" in q for q in v.feature_qualifiers["pseudogenic_CDS"])


def test_duplicate_so_term_prefers_concrete():
    v = load_vocab()
    # LINE_element appears twice: /mobile_element_type="LINE*" and ="LINE"; concrete (no '*') wins
    quals = v.feature_qualifiers["LINE_element"]
    assert quals
    assert all("*" not in q for q in quals)


def test_placeholder_qualifier_not_displaced_by_empty():
    v = load_vocab()
    # mobile_genetic_element has a placeholder row and an empty row; the placeholder must survive
    quals = v.feature_qualifiers["mobile_genetic_element"]
    assert quals  # not empty
    assert any("<" in q or "*" in q for q in quals)  # the placeholder is retained


def test_insdc_map_unchanged_for_3a():
    v = load_vocab()
    assert "CDS" in v.feature_types
    assert v.insdc_map.get("ncRNA_gene") == "ncRNA"
    # term with no qualifier maps to empty tuple
    assert v.feature_qualifiers.get("CDS", ()) == ()
