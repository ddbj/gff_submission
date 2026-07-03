from Bio.Seq import Seq
from ddbj_gff import parse
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


def test_emit_mrna_false_suppresses_mrna_keeps_cds():
    # organelle profile: CDS only, no mRNA feature
    gene = Feature("g", "S", "gene", [Span("c", 1, 9, "+")], {"ID": ["g"]}, [])
    mrna = Feature("g.t1", "S", "mRNA", [Span("c", 1, 9, "+")], {"ID": ["g.t1"]}, ["g"])
    cds = Feature("cds", "S", "CDS", [Span("c", 1, 9, "+", 0)], {"ID": ["cds"]}, ["g.t1"])
    mrna.children = [cds]
    gene.children = [mrna]
    cfg2 = MssConfig(source={}, transl_table=1, product_default="hypothetical protein", emit_mrna=False)
    feats = build_gene_features(gene, "nonredundant", LocusTagAssigner("L", 6, 10, 10),
                                Seq("ATGAAATAA"), cfg2, [])
    keys = [f.key for f in feats]
    assert "mRNA" not in keys and "CDS" in keys


def test_gene_with_internal_stop_emits_only_misc_feature():
    gene, mrna, genome = _gene_with_internal_stop()
    assigner = LocusTagAssigner("L", 6, 10, 10)
    feats = build_gene_features(gene, "nonredundant", assigner, genome, cfg(), [])
    keys = [f.key for f in feats]
    assert keys == ["misc_feature"]  # mRNA/CDS は出さない


def test_multi_transcript_all_broken_emits_single_misc_feature():
    GFF = (
        "##gff-version 3\n"
        "c1\tS\tgene\t1\t9\t.\t+\t.\tID=g\n"
        "c1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g.t1;Parent=g\n"
        "c1\tS\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=g.t1\n"
        "c1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g.t2;Parent=g\n"
        "c1\tS\tCDS\t1\t9\t.\t+\t0\tID=cds2;Parent=g.t2\n"
    )
    doc = parse(GFF)
    gene = [f for f in doc.roots if f.type == "gene"][0]
    assigner = LocusTagAssigner("L", 6, 10, 10)
    feats = build_gene_features(gene, "nonredundant", assigner, Seq("ATGTAAAAA"), cfg(), [])
    keys = [f.key for f in feats]
    assert keys == ["misc_feature"]  # 全転写産物 broken -> misc は1つ
    tags = [q.value for f in feats for q in f.qualifiers if q.key == "locus_tag"]
    assert len(tags) == len(set(tags))  # locus_tag 重複なし


def test_mixed_transcript_prefers_valid_cds_no_misc():
    GFF = (
        "##gff-version 3\n"
        "c1\tS\tgene\t1\t12\t.\t+\t.\tID=g\n"
        "c1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g.t1;Parent=g\n"
        "c1\tS\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=g.t1\n"
        "c1\tS\tmRNA\t1\t12\t.\t+\t.\tID=g.t2;Parent=g\n"
        "c1\tS\tCDS\t1\t12\t.\t+\t0\tID=cds2;Parent=g.t2\n"
    )
    doc = parse(GFF)
    gene = [f for f in doc.roots if f.type == "gene"][0]
    assigner = LocusTagAssigner("L", 6, 10, 10)
    # genome ATGAAATAACCC: t1 CDS 1..9 = ATGAAATAA (MK*, clean terminal stop); t2 CDS 1..12 = ATGAAATAACCC (MK*P internal stop)
    feats = build_gene_features(gene, "nonredundant", assigner, Seq("ATGAAATAACCC"), cfg(), [])
    keys = [f.key for f in feats]
    assert "misc_feature" not in keys  # 有効アイソフォームがあるので misc は出さない
    assert "CDS" in keys and "mRNA" in keys
