from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import Counter
from ddbj_gff.flatfile.convert import synthesize_features

FIX = "tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk"


class _Rec:
    """Minimal stand-in for a BioPython SeqRecord: synthesize_features only
    reads rec.features; each feature needs .type, .qualifiers, and .location
    (with .parts, used by bio_location_to_spans)."""

    def __init__(self, features):
        self.features = features


def _feat(ftype, start, end, strand=1, qualifiers=None):
    return SeqFeature(FeatureLocation(start, end, strand=strand), type=ftype,
                       qualifiers=qualifiers or {})


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


def test_locus_tag_less_features_get_separate_genes():
    """Two independent locus_tag-less CDS (no containing mRNA, no shared
    locus_tag) must NOT collapse into a single mega-gene: each gets its own
    synthetic group key and thus its own gene."""
    c1 = _feat("CDS", 0, 90, qualifiers={"protein_id": ["P1"]})
    c2 = _feat("CDS", 200, 260, qualifiers={"protein_id": ["P2"]})
    feats = synthesize_features(_Rec([c1, c2]), "SEQ1")

    genes = [f for f in feats if f.type == "gene"]
    assert len(genes) == 2
    assert genes[0].id != genes[1].id
    # neither carries a locus_tag attribute, since none was present in the input
    for g in genes:
        assert "locus_tag" not in g.attributes
    # each CDS still lands under its own gene (via its own mRNA), not merged
    cds = [f for f in feats if f.type == "CDS"]
    assert len(cds) == 2
    assert {f._first("protein_id") for f in cds} == {"P1", "P2"}


def test_protein_id_less_cds_on_shared_mrna_get_distinct_ids():
    """Two CDS lacking /protein_id that both pair to the same mRNA must not
    receive the same synthesized CDS GFF ID (would silently span-merge on a
    reparse of the emitted GFF)."""
    mrna = _feat("mRNA", 0, 300, qualifiers={"locus_tag": ["LT1"]})
    c1 = _feat("CDS", 0, 150, qualifiers={"locus_tag": ["LT1"]})    # no protein_id
    c2 = _feat("CDS", 200, 300, qualifiers={"locus_tag": ["LT1"]})  # no protein_id
    feats = synthesize_features(_Rec([mrna, c1, c2]), "SEQ1")

    cds = [f for f in feats if f.type == "CDS"]
    assert len(cds) == 2
    assert cds[0].id != cds[1].id
    assert len({f.id for f in cds}) == 2
