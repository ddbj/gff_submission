from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import build_entry_features

GFF = """##gff-version 3
c1\ttRNAscan-SE\tpseudogene\t10\t80\t.\t-\t.\tID=p1;gene_biotype=pseudogene
c1\ttRNAscan-SE\ttRNA\t100\t170\t.\t+\t.\tID=t1;isotype=Ala
"""


def test_pseudogene_excluded():
    doc = parse(GFF)
    diags = []
    per = build_entry_features(doc, {"c1": Seq("A" * 300)},
                               MssConfig(source={}, locus_tag_prefix="P"), diags)
    keys = [f.key for f in per["c1"]]
    assert "pseudogene" not in keys
    assert keys == ["tRNA"]
    assert any(d.code == "pseudogene-skipped" for d in diags)
