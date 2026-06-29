from Bio.Seq import Seq
from ddbj_gff.model import Feature, Span
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import build_noncoding_features


def cfg():
    return MssConfig(source={})


def test_mirna_gene_maps_to_precursor_and_ncrna():
    gene = Feature("Mp1g00675", "S", "gene", [Span("chr1", 603255, 603422, "-")],
                   {"gene_biotype": ["miRNA"]}, [])
    pre = Feature("Mp1g00675.pre", "S", "pre_miRNA", [Span("chr1", 603255, 603422, "-")],
                  {"Note": ["Mpo-pre-miR11669"]}, [])
    m1 = Feature("Mp1g00675.1", "S", "miRNA", [Span("chr1", 603382, 603402, "-")],
                 {"Note": ["Mpo-miR11669.1"]}, [])
    m2 = Feature("Mp1g00675.2", "S", "miRNA", [Span("chr1", 603384, 603405, "-")],
                 {"Note": ["Mpo-miR11669.2"]}, [])
    gene.children = [pre, m1, m2]
    feats = build_noncoding_features(gene, "PFX_000010", 700000, cfg())
    assert [f.key for f in feats] == ["precursor_RNA", "ncRNA", "ncRNA"]
    assert feats[0].location == "complement(603255..603422)"
    nc = {q.key: q.value for q in feats[1].qualifiers}
    assert nc["locus_tag"] == "PFX_000010"
    assert nc["ncRNA_class"] == "miRNA"
    assert any(q.key == "note" and q.value == "Mpo-miR11669.1" for q in feats[1].qualifiers)
    assert any(q.key == "note" and "submitter_gene_id: Mp1g00675" in q.value for q in feats[1].qualifiers)


def test_no_recognized_rna_children_returns_empty():
    gene = Feature("g", "S", "gene", [Span("chr1", 1, 9, "+")], {}, [])
    gene.children = []
    assert build_noncoding_features(gene, "PFX_000010", 1000, cfg()) == []
