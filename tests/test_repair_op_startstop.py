# tests/test_repair_op_startstop.py
from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.repair.context import RepairContext
from ddbj_gff.repair.registry import get_operation

HDR = "##gff-version 3\n##sequence-region s 1 100000\n"

# CDS 1..12 on + strand. Sequence ATG AAA GTT AAA -> starts with ATG (start codon),
# does NOT end in a stop codon -> 3' partial only.
GFF = HDR + "\n".join([
    "s\tsrc\tgene\t1\t12\t.\t+\t.\tID=g1;locus_tag=X_0001",
    "s\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=m1;Parent=g1",
    "s\tsrc\tCDS\t1\t12\t.\t+\t0\tID=c1;Parent=m1;transl_table=1",
]) + "\n"
SEQ = {"s": Seq("ATGAAAGTTAAA")}


def _ctx():
    return RepairContext(sequences=SEQ, transl_table=1)


def test_detect_three_prime_partial_missing_stop():
    doc = parse(GFF)
    op = get_operation("missing-start-stop-to-partial-cds")
    cands = op.detect(doc, _ctx())
    assert len(cands) == 1
    assert cands[0].payload["five"] is False    # ATG is a start codon
    assert cands[0].payload["three"] is True     # AAA is not a stop


def test_detect_five_prime_partial_missing_start():
    # sequence starts with CTT (not a start codon) and ends with a stop (TAA)
    doc = parse(GFF)
    op = get_operation("missing-start-stop-to-partial-cds")
    ctx = RepairContext(sequences={"s": Seq("CTTAAAGTTTAA")}, transl_table=1)
    cands = op.detect(doc, ctx)
    assert cands[0].payload["five"] is True      # CTT not a start codon
    assert cands[0].payload["three"] is False    # TAA is a stop


def test_apply_sets_end_range_on_cds():
    doc = parse(GFF)
    op = get_operation("missing-start-stop-to-partial-cds")
    op.apply(doc, _ctx(), None)
    c1 = doc.feature_index["c1"]
    assert c1.attributes.get("partial") == ["true"]
    assert c1.attributes.get("end_range") == ["12,."]
    assert "start_range" not in c1.attributes


def test_requires_sequence_flag():
    assert get_operation("missing-start-stop-to-partial-cds").requires_sequence is True


def test_no_false_three_prime_partial_when_inframe_stop_with_trailing_base():
    # CDS 1..10 on + strand. Sequence ATGAAATAAG -> in-frame codons ATG AAA TAA
    # plus a trailing base G past the stop. Starts with ATG (start codon) and
    # ends in-frame with TAA (stop codon) -> COMPLETE; must not be flagged 3'-partial.
    gff = HDR + "\n".join([
        "s\tsrc\tgene\t1\t10\t.\t+\t.\tID=g2;locus_tag=X_0002",
        "s\tsrc\tmRNA\t1\t10\t.\t+\t.\tID=m2;Parent=g2",
        "s\tsrc\tCDS\t1\t10\t.\t+\t0\tID=c2;Parent=m2;transl_table=1",
    ]) + "\n"
    doc = parse(gff)
    ctx = RepairContext(sequences={"s": Seq("ATGAAATAAG")}, transl_table=1)
    op = get_operation("missing-start-stop-to-partial-cds")
    assert op.detect(doc, ctx) == []
