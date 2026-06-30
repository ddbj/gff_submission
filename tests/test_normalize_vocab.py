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


def test_insdc_map_unchanged_for_3a():
    v = load_vocab()
    assert "CDS" in v.feature_types
    assert v.insdc_map.get("ncRNA_gene") == "ncRNA"
    # term with no qualifier maps to empty tuple
    assert v.feature_qualifiers.get("CDS", ()) == ()
