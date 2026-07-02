from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import build_entry_features

GFF = """##gff-version 3
c1\tInfernal\tncRNA\t10\t50\t.\t+\t.\tID=n1;Name=U1;Dbxref=RFAM:RF00003;note=hit
c1\ttRNAscan-SE\ttRNA\t60\t90\t.\t-\t.\tID=t1;isotype=Thr;anticodon=CGT
c1\tpybarrnap\trRNA\t100\t200\t.\t+\t.\tID=r1;product=18S ribosomal RNA
"""


def cfg():
    return MssConfig(source={}, locus_tag_prefix="PFX")


def test_parentless_rna_emitted_as_toplevel_with_locus_tag():
    doc = parse(GFF)
    per = build_entry_features(doc, {"c1": Seq("A" * 300)}, cfg(), [])
    keys = [f.key for f in per["c1"]]
    assert keys == ["ncRNA", "tRNA", "rRNA"]  # 位置順
    nc = {q.key: q.value for q in per["c1"][0].qualifiers}
    assert nc["locus_tag"] == "PFX_000010"
    assert nc["ncRNA_class"] == "other"  # RFAM 由来の未知クラスは other
    assert any(q.key == "db_xref" and q.value == "RFAM:RF00003" for q in per["c1"][0].qualifiers)
    rr = {q.key: q.value for q in per["c1"][2].qualifiers}
    assert rr["product"] == "18S ribosomal RNA"
    tr = {q.key: q.value for q in per["c1"][1].qualifiers}
    assert tr["product"] == "tRNA-Thr"  # isotype から導出


def test_parentless_ncrna_lowercase_note_preserved():
    doc = parse(GFF)
    per = build_entry_features(doc, {"c1": Seq("A" * 300)}, cfg(), [])
    nc = per["c1"][0]
    assert any(q.key == "note" and q.value == "hit" for q in nc.qualifiers)
