from Bio import SeqIO
from collections import Counter
from ddbj_gff.flatfile.convert import synthesize_features

FIX = "tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk"


def test_synthesis_counts_and_parentage():
    rec = SeqIO.read(FIX, "genbank")
    feats = synthesize_features(rec, "BDQV01000200.1")
    counts = Counter(f.type for f in feats)
    assert counts["gene"] == 2
    assert counts["mRNA"] == 3
    assert counts["exon"] == 7      # 2 + 3 + 2
    assert counts["CDS"] == 3

    by_id = {f.id: f for f in feats}
    genes = [f for f in feats if f.type == "gene"]
    gene_ids = {g.id for g in genes}
    # every mRNA parents a gene; every exon/CDS parents an mRNA
    mrnas = [f for f in feats if f.type == "mRNA"]
    for m in mrnas:
        assert m.parent_ids and m.parent_ids[0] in gene_ids
    mrna_ids = {m.id for m in mrnas}
    for f in feats:
        if f.type in ("exon", "CDS"):
            assert f.parent_ids and f.parent_ids[0] in mrna_ids


def test_altsplice_cds_paired_by_containment():
    rec = SeqIO.read(FIX, "genbank")
    feats = synthesize_features(rec, "BDQV01000200.1")
    mrnas = {f.id: f for f in feats if f.type == "mRNA"}
    # CUMW_191340 has two transcripts: CDS-A (ends 2241) pairs with mRNA spanning ..2597,
    # CDS-B (2245..2762) pairs with the mRNA spanning 2245..2762
    cds = [f for f in feats if f.type == "CDS" and f._first("locus_tag") == "CUMW_191340"]
    assert len(cds) == 2
    for c in cds:
        parent = mrnas[c.parent_ids[0]]
        cs = c.ordered_spans()
        ms = parent.ordered_spans()
        lo_m, hi_m = min(s.start for s in ms), max(s.end for s in ms)
        # CDS fully within its paired mRNA
        assert lo_m <= min(s.start for s in cs) and max(s.end for s in cs) <= hi_m
    # the two CDS are paired to DIFFERENT mRNAs
    assert cds[0].parent_ids[0] != cds[1].parent_ids[0]
