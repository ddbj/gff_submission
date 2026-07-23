from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.repair.context import RepairContext
from ddbj_gff.repair.registry import get_operation

HDR = "##gff-version 3\n##sequence-region s 1 100000\n"
# CDS 1..12: ATG TAA AAA TAA -> M * K  (internal stop at aa index 1)
GFF = HDR + "\n".join([
    "s\tsrc\tgene\t1\t12\t.\t+\t.\tID=g1;locus_tag=X_0001",
    "s\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=m1;Parent=g1",
    "s\tsrc\tCDS\t1\t12\t.\t+\t0\tID=c1;Parent=m1;transl_table=1",
]) + "\n"
SEQ = {"s": Seq("ATGTAAAAATAA")}


def _ctx():
    return RepairContext(sequences=SEQ, transl_table=1)


def test_detect_internal_stop():
    doc = parse(GFF)
    op = get_operation("internal-stop-to-misc")
    cands = op.detect(doc, _ctx())
    assert len(cands) == 1 and cands[0].feature_id == "c1"


def test_apply_retypes_cds_to_misc_feature_with_note():
    doc = parse(GFF)
    op = get_operation("internal-stop-to-misc")
    op.apply(doc, _ctx(), None)
    c1 = doc.feature_index["c1"]
    assert c1.type == "misc_feature"
    assert any("internal stop" in n for n in c1.note)
    # gene/mRNA and links intact
    g1 = doc.feature_index["g1"]
    m1 = doc.feature_index["m1"]
    assert g1.type == "gene" and m1.type == "mRNA"
    assert c1 in m1.children


def test_no_candidate_for_clean_cds():
    doc = parse(GFF)
    op = get_operation("internal-stop-to-misc")
    ctx = RepairContext(sequences={"s": Seq("ATGAAAGTTTAA")}, transl_table=1)  # M K V
    assert op.detect(doc, ctx) == []
