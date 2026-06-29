import pytest
from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.errors import GffParseError
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import convert

GFF = """\
##gff-version 3
chr1\tS\tgene\t1\t9\t.\t+\t.\tID=g1;gene=MpX
chr1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g1.1;Parent=g1
chr1\tS\texon\t1\t9\t.\t+\t.\tID=e1;Parent=g1.1
chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=g1.1
"""


def cfg():
    return MssConfig(source={"organism": "Foo bar"}, locus_tag_prefix="PFX")


def test_convert_builds_entry_with_source_and_gene():
    doc = parse(GFF)
    seqs = {"chr1": Seq("ATGAAATAA")}
    mss, diags = convert(doc, seqs, cfg(), ["COMMON\tDBLINK"])
    assert mss.common_rows == ["COMMON\tDBLINK"]
    assert len(mss.entries) == 1
    entry = mss.entries[0]
    assert entry.name == "chr1"
    keys = [f.key for f in entry.features]
    assert keys == ["source", "mRNA", "CDS"]
    cds = entry.features[2]
    assert {q.key: q.value for q in cds.qualifiers}["locus_tag"] == "PFX_000010"


def test_missing_sequence_is_error_and_skips_entry():
    doc = parse(GFF)
    mss, diags = convert(doc, {}, cfg(), ["COMMON"])
    assert any(d.code == "missing-sequence" for d in diags)
    assert mss.entries == []
    with pytest.raises(GffParseError):
        convert(doc, {}, cfg(), ["COMMON"], strict=True)


def test_multi_transcript_warns_and_keeps_first():
    gff = GFF + (
        "chr1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g1.2;Parent=g1\n"
        "chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=c2;Parent=g1.2\n"
    )
    doc = parse(gff)
    seqs = {"chr1": Seq("ATGAAATAA")}
    mss, diags = convert(doc, seqs, cfg(), ["COMMON"])
    assert any(d.code == "multi-transcript" for d in diags)
    # only one mRNA + one CDS kept
    keys = [f.key for f in mss.entries[0].features]
    assert keys.count("mRNA") == 1 and keys.count("CDS") == 1


def test_locus_tag_unique_across_entries():
    # two seqids, neither gene has a locus_tag attr -> counter must NOT restart per entry
    gff = (
        "##gff-version 3\n"
        "chr1\tS\tgene\t1\t9\t.\t+\t.\tID=a\n"
        "chr1\tS\tmRNA\t1\t9\t.\t+\t.\tID=a.1;Parent=a\n"
        "chr1\tS\texon\t1\t9\t.\t+\t.\tID=ae;Parent=a.1\n"
        "chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=ac;Parent=a.1\n"
        "chr2\tS\tgene\t1\t9\t.\t+\t.\tID=b\n"
        "chr2\tS\tmRNA\t1\t9\t.\t+\t.\tID=b.1;Parent=b\n"
        "chr2\tS\texon\t1\t9\t.\t+\t.\tID=be;Parent=b.1\n"
        "chr2\tS\tCDS\t1\t9\t.\t+\t0\tID=bc;Parent=b.1\n"
    )
    doc = parse(gff)
    seqs = {"chr1": Seq("ATGAAATAA"), "chr2": Seq("ATGAAATAA")}
    mss, diags = convert(doc, seqs, cfg(), ["COMMON"])
    tags = [q.value for e in mss.entries for f in e.features if f.key == "CDS"
            for q in f.qualifiers if q.key == "locus_tag"]
    assert tags == ["PFX_000010", "PFX_000020"]  # monotonic across entries, no duplicate


def test_no_exon_no_cds_mrna_skipped():
    gff = (
        "##gff-version 3\n"
        "chr1\tS\tgene\t1\t9\t.\t+\t.\tID=g\n"
        "chr1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g.1;Parent=g\n"  # no exon/CDS children
    )
    doc = parse(gff)
    mss, diags = convert(doc, {"chr1": Seq("ATGAAATAA")}, cfg(), ["COMMON"])
    assert any(d.code == "no-exon" for d in diags)
    # no crash; the gene produced no mRNA/CDS features
    assert all(f.key == "source" for e in mss.entries for f in e.features)
