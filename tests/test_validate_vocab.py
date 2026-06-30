from ddbj_gff.validate.vocab import Vocab, load_vocab


def test_load_vocab_feature_types():
    v = load_vocab()
    assert isinstance(v, Vocab)
    for term in ("CDS", "mRNA", "gene", "exon"):
        assert term in v.feature_types
    assert "totally_not_a_real_so_term" not in v.feature_types


def test_load_vocab_insdc_map_and_dbxref():
    v = load_vocab()
    assert v.insdc_map.get("CDS") == "CDS"   # SO term CDS maps to INSDC feature CDS
    assert "GenBank" in v.dbxref_dbtags
    assert "GeneID" in v.dbxref_dbtags
