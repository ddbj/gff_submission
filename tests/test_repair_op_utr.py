from ddbj_gff import parse
from ddbj_gff.repair.context import RepairContext
from ddbj_gff.repair.registry import get_operation

HDR = "##gff-version 3\n##sequence-region s 1 100000\n"

# mRNA 1..1000; exon 1..1000; CDS 1..600  -> 3' UTR present (601..1000),
# 5' UTR absent (exon_lo==cds_lo==1) -> 5' partial on + strand.
GFF = HDR + "\n".join([
    "s\tsrc\tgene\t1\t1000\t.\t+\t.\tID=g1;locus_tag=X_0001",
    "s\tsrc\tmRNA\t1\t1000\t.\t+\t.\tID=m1;Parent=g1",
    "s\tsrc\texon\t1\t1000\t.\t+\t.\tID=e1;Parent=m1",
    "s\tsrc\tCDS\t1\t600\t.\t+\t0\tID=c1;Parent=m1;transl_table=1",
]) + "\n"


def test_detect_finds_five_prime_partial():
    doc = parse(GFF)
    op = get_operation("utr-absent-to-partial-mrna")
    cands = op.detect(doc, RepairContext())
    assert len(cands) == 1
    c = cands[0]
    assert c.feature_id == "m1"
    assert c.payload["five"] is True and c.payload["three"] is False


def test_apply_sets_partial_attrs_on_mrna():
    doc = parse(GFF)
    op = get_operation("utr-absent-to-partial-mrna")
    changes = op.apply(doc, RepairContext(), None)
    m1 = doc.feature_index["m1"]
    assert m1.attributes.get("partial") == ["true"]
    assert m1.attributes.get("start_range") == [".,1"]
    assert "end_range" not in m1.attributes
    assert len(changes) == 1


def test_apply_is_idempotent():
    doc = parse(GFF)
    op = get_operation("utr-absent-to-partial-mrna")
    op.apply(doc, RepairContext(), None)
    assert op.detect(doc, RepairContext()) == []   # already partial -> no candidate
