from Bio.Seq import Seq
from ddbj_gff.model import Feature, Span
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import build_cds_feature, build_gene_features
from ddbj_gff.mss.locus_tag import LocusTagAssigner


def cfg():
    return MssConfig(source={}, transl_table=1, product_default="hypothetical protein")


def _gene_with_internal_stop():
    genome = Seq("ATGTAAAAA")  # ATG TAA AAA -> M * K  内部stop
    gene = Feature("g", "S", "gene", [Span("c", 1, 9, "+")], {"ID": ["g"]}, [])
    mrna = Feature("g.t1", "S", "mRNA", [Span("c", 1, 9, "+")], {"ID": ["g.t1"]}, ["g"])
    cds = Feature("cds", "S", "CDS", [Span("c", 1, 9, "+", 0)], {"ID": ["cds"]}, ["g.t1"])
    mrna.children = [cds]
    gene.children = [mrna]
    return gene, mrna, genome


def test_internal_stop_returns_misc_feature():
    gene, mrna, genome = _gene_with_internal_stop()
    f = build_cds_feature(mrna, gene, "L_1", genome, cfg(), [])
    assert f.key == "misc_feature"
    q = {x.key: x.value for x in f.qualifiers}
    assert "product" not in q
    assert q["locus_tag"] == "L_1"
    assert any(x.key == "note" and "internal stop" in x.value.lower() for x in f.qualifiers)


def test_gene_with_internal_stop_emits_only_misc_feature():
    gene, mrna, genome = _gene_with_internal_stop()
    assigner = LocusTagAssigner("L", 6, 10, 10)
    feats = build_gene_features(gene, "nonredundant", assigner, genome, cfg(), [])
    keys = [f.key for f in feats]
    assert keys == ["misc_feature"]  # mRNA/CDS は出さない
