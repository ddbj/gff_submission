from ddbj_gff import parse
from ddbj_gff.normalize.passes import pass_wrap_cds_in_mrna, NormalizeContext
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.validate.vocab import load_vocab

GFF = """##gff-version 3
CP\tLiftoff\tgene\t1\t9\t.\t+\t.\tID=gene-x
CP\tLiftoff\tCDS\t1\t9\t.\t+\t0\tID=cds-x;Parent=gene-x;product=foo
MT\tLiftoff\tgene\t1\t9\t.\t+\t.\tID=gene-t
MT\tLiftoff\ttRNA\t1\t9\t.\t+\t.\tID=rna-t;Parent=gene-t
MT\tLiftoff\texon\t1\t9\t.\t+\t.\tID=ex-t;Parent=rna-t
"""


def _ctx():
    return NormalizeContext(vocab=load_vocab(), seq_lengths=None, config=NormalizeConfig())


def test_gene_with_direct_cds_gets_mrna_layer():
    doc = parse(GFF)
    pass_wrap_cds_in_mrna(doc, _ctx())
    gx = doc.get("gene-x")
    assert [c.type for c in gx.children] == ["mRNA"]
    mrna = gx.children[0]
    assert mrna.id == "gene-x.mrna"
    assert [c.type for c in mrna.children] == ["CDS"]
    cds = mrna.children[0]
    assert cds.id == "cds-x" and cds.parent_ids == ["gene-x.mrna"]
    assert doc.get("gene-x.mrna") is mrna  # indexed


def test_trna_gene_left_untouched():
    doc = parse(GFF)
    pass_wrap_cds_in_mrna(doc, _ctx())
    gt = doc.get("gene-t")
    assert [c.type for c in gt.children] == ["tRNA"]  # no mRNA inserted


def test_disabled_by_config():
    doc = parse(GFF)
    cfg = NormalizeConfig(wrap_cds_in_mrna=False)
    ctx = NormalizeContext(vocab=load_vocab(), seq_lengths=None, config=cfg)
    pass_wrap_cds_in_mrna(doc, ctx)
    assert [c.type for c in doc.get("gene-x").children] == ["CDS"]
