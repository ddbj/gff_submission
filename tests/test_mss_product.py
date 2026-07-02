from Bio.Seq import Seq
from ddbj_gff.model import Feature, Span
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import _product
from ddbj_gff.mss.product_map import load_product_map


def _mrna(mid, gid, product=None):
    attr = {"product": [product]} if product else {}
    mrna = Feature(mid, "S", "mRNA", [Span("c", 1, 9, "+")], attr, [])
    gene = Feature(gid, "S", "gene", [Span("c", 1, 9, "+")], {}, [])
    return mrna, gene


def test_product_map_hit_by_transcript_id():
    mrna, gene = _mrna("g1.t1", "g1")
    cfg = MssConfig(source={}, product_map={"g1.t1": "tubulin-tyrosine ligase"})
    assert _product(mrna, gene, cfg) == "tubulin-tyrosine ligase"


def test_product_map_hit_by_gene_id_when_transcript_misses():
    mrna, gene = _mrna("g1.t2", "g1")
    cfg = MssConfig(source={}, product_map={"g1": "some kinase"})
    assert _product(mrna, gene, cfg) == "some kinase"


def test_product_falls_back_to_col9():
    mrna, gene = _mrna("g1.t1", "g1", product="50S ribosomal protein L5")
    cfg = MssConfig(source={})
    assert _product(mrna, gene, cfg) == "50S ribosomal protein L5"


def test_product_defaults_to_hypothetical_no_gene_name_fallback():
    mrna, gene = _mrna("g1.t1", "g1")
    mrna.attributes["gene"] = ["MpX"]  # 旧仕様なら "protein MpX" だったが廃止
    cfg = MssConfig(source={}, product_default="hypothetical protein")
    assert _product(mrna, gene, cfg) == "hypothetical protein"


def test_load_product_map_reads_tsv(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("g1.t1\ttubulin\ng1\ttubulin\n", encoding="utf-8")
    m = load_product_map(str(p))
    assert m == {"g1.t1": "tubulin", "g1": "tubulin"}
