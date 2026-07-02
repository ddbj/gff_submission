from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import build_entry_features, convert

GFF = """##gff-version 3
c1\tS\tgene\t1\t9\t.\t+\t.\tID=g1
c1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1
c1\tS\texon\t1\t9\t.\t+\t.\tID=e1;Parent=g1.t1
c1\tS\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=g1.t1
"""


def cfg():
    return MssConfig(source={"organism": "x", "mol_type": "genomic DNA"},
                     locus_tag_prefix="PFX")


def test_build_entry_features_has_no_source():
    doc = parse(GFF)
    per = build_entry_features(doc, {"c1": Seq("ATGAAATAA")}, cfg(), [])
    assert set(per) == {"c1"}
    assert all(f.key != "source" for f in per["c1"])
    assert any(f.key == "CDS" for f in per["c1"])


def test_convert_still_emits_source_first():
    doc = parse(GFF)
    mss, _ = convert(doc, {"c1": Seq("ATGAAATAA")}, cfg(), ["COMMON"])
    assert mss.entries[0].features[0].key == "source"
