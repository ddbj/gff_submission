from Bio.Seq import Seq
from ddbj_gff.model import Feature, Span
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.locus_tag import LocusTagAssigner
from ddbj_gff.mss.convert import build_gene_features


def cfg():
    return MssConfig(source={}, locus_tag_prefix="PFX")


def assigner():
    return LocusTagAssigner.from_config(cfg())


def two_transcript_gene(cds_same=True):
    # gene with two mRNAs; both share identical CDS (1..9) but differ in exon (UTR)
    gene = Feature("g", "S", "gene", [Span("c", 1, 40, "+")], {}, [])
    t1 = Feature("g.1", "S", "mRNA", [Span("c", 1, 30, "+")], {}, [])
    t1.children = [Feature("e1", "S", "exon", [Span("c", 1, 30, "+")], {}, []),
                   Feature("c1", "S", "CDS", [Span("c", 1, 9, "+", 0)], {}, [])]
    t2 = Feature("g.2", "S", "mRNA", [Span("c", 1, 40, "+")], {}, [])
    cds2 = Span("c", 1, 9, "+", 0) if cds_same else Span("c", 1, 12, "+", 0)
    t2.children = [Feature("e2", "S", "exon", [Span("c", 1, 40, "+")], {}, []),
                   Feature("c2", "S", "CDS", [cds2], {}, [])]
    gene.children = [t1, t2]
    return gene


def keys(feats):
    return [f.key for f in feats]


def test_minimal_keeps_one_transcript():
    g = two_transcript_gene()
    feats = build_gene_features(g, "minimal", assigner(), Seq("ATGAAATAA" + "C" * 31), cfg(), [])
    assert keys(feats) == ["mRNA", "CDS"]


def test_full_keeps_all_transcripts_and_cds():
    g = two_transcript_gene()
    feats = build_gene_features(g, "full", assigner(), Seq("ATGAAATAA" + "C" * 31), cfg(), [])
    assert keys(feats) == ["mRNA", "CDS", "mRNA", "CDS"]


def test_nonredundant_dedupes_shared_cds():
    g = two_transcript_gene(cds_same=True)
    feats = build_gene_features(g, "nonredundant", assigner(), Seq("ATGAAATAA" + "C" * 31), cfg(), [])
    assert keys(feats) == ["mRNA", "mRNA", "CDS"]   # both mRNAs, one shared CDS
    cds = feats[2]
    note = [q.value for q in cds.qualifiers if q.key == "note"][0]
    assert "g.1" in note and "g.2" in note            # note lists both source transcripts


def test_nonredundant_keeps_distinct_cds():
    g = two_transcript_gene(cds_same=False)
    # ATGAAAGGGCCC: t1 CDS 1..9 = ATGAAAGGG (MKG, 3' partial), t2 CDS 1..12 = ATGAAAGGGCCC (MKGP, 3' partial)
    # neither has an internal stop codon
    feats = build_gene_features(g, "nonredundant", assigner(), Seq("ATGAAAGGGCCC" + "C" * 28), cfg(), [])
    assert keys(feats) == ["mRNA", "mRNA", "CDS", "CDS"]  # different CDS -> not deduped


def test_noncoding_gene_dispatch():
    gene = Feature("m", "S", "gene", [Span("c", 1, 50, "-")], {"gene_biotype": ["miRNA"]}, [])
    gene.children = [Feature("m.pre", "S", "pre_miRNA", [Span("c", 1, 50, "-")], {}, []),
                     Feature("m.1", "S", "miRNA", [Span("c", 10, 30, "-")], {}, [])]
    feats = build_gene_features(gene, "nonredundant", assigner(), Seq("A" * 100), cfg(), [])
    assert keys(feats) == ["precursor_RNA", "ncRNA"]


def test_gene_with_no_rna_or_mrna_warns():
    gene = Feature("x", "S", "gene", [Span("c", 1, 9, "+")], {}, [])
    gene.children = []
    diags = []
    feats = build_gene_features(gene, "nonredundant", assigner(), Seq("A" * 9), cfg(), diags)
    assert feats == []
    assert any(d.code == "no-rna" for d in diags)


def test_unrecognized_rna_gene_produces_misc_rna_not_no_rna():
    gene = Feature("g", "S", "gene", [Span("c", 1, 50, "+")], {}, [])
    gene.children = [Feature("g.1", "S", "lnc_RNA", [Span("c", 1, 50, "+")], {}, [])]
    diags = []
    feats = build_gene_features(gene, "nonredundant", assigner(), Seq("A" * 100), cfg(), diags)
    assert [f.key for f in feats] == ["misc_RNA"]
    assert not any(d.code == "no-rna" for d in diags)


def test_transcript_without_exon_or_cds_skipped():
    gene = Feature("g", "S", "gene", [Span("c", 1, 30, "+")], {}, [])
    t = Feature("g.1", "S", "mRNA", [Span("c", 1, 30, "+")], {}, [])
    t.children = []
    gene.children = [t]
    diags = []
    feats = build_gene_features(gene, "full", assigner(), Seq("A" * 30), cfg(), diags)
    assert feats == []
    assert any(d.code == "no-exon" for d in diags)
