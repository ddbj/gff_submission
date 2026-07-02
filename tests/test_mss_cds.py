from Bio.Seq import Seq
from ddbj_gff.model import Feature, Span
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import build_cds_feature


def cfg():
    return MssConfig(source={}, transl_table=1, product_default="hypothetical protein")


def mrna_with_cds(cds_spans, strand="+", phase0=0, mrna_attr=None):
    gene = Feature("Mp1g1", "S", "gene", [Span("chr1", 1, 999, strand)], {}, [])
    mrna = Feature("Mp1g1.1", "S", "mRNA", [Span("chr1", 1, 999, strand)], mrna_attr or {}, [])
    children = []
    for i, (s, e) in enumerate(cds_spans):
        ph = phase0 if i == 0 else 0
        children.append(Feature(f"cds{i}", "S", "CDS", [Span("chr1", s, e, strand, ph)], {}, []))
    mrna.children = children
    return gene, mrna


def test_complete_cds_plus():
    genome = Seq("ATGAAATAA")  # chr1 1..9
    gene, mrna = mrna_with_cds([(1, 9)], strand="+", phase0=0)
    diags = []
    f = build_cds_feature(mrna, gene, "PFX_000010", genome, cfg(), diags)
    assert f.key == "CDS" and f.location == "1..9"
    q = {x.key: x.value for x in f.qualifiers}
    assert q["locus_tag"] == "PFX_000010"
    assert q["transl_table"] == "1"
    assert q["codon_start"] == "1"
    assert q["product"] == "hypothetical protein"
    assert "submitter_gene_id: Mp1g1" in q["note"]
    assert diags == []  # MK* is a clean translation


def test_five_prime_partial_no_start_codon():
    genome = Seq("GGGAAATAA")  # first codon GGG is not a start codon
    gene, mrna = mrna_with_cds([(1, 9)], strand="+", phase0=0)
    f = build_cds_feature(mrna, gene, "PFX_000010", genome, cfg(), [])
    assert f.location == "<1..9"  # 5' partial


def test_codon_start_from_phase():
    genome = Seq("AAAAAAAAA")
    gene, mrna = mrna_with_cds([(1, 9)], strand="+", phase0=2)  # phase 2 -> codon_start 3
    diags = []
    f = build_cds_feature(mrna, gene, "PFX_000010", genome, cfg(), diags)
    q = {x.key: x.value for x in f.qualifiers}
    assert q["codon_start"] == "3"
    assert f.location.startswith("<")  # codon_start != 1 implies 5' partial
    assert any(d.code == "translation-not-multiple-of-3" for d in diags)


def test_product_defaults_to_hypothetical_when_only_gene_name():
    genome = Seq("ATGAAATAA")
    gene, mrna = mrna_with_cds([(1, 9)], strand="+", phase0=0, mrna_attr={"gene": ["MpX"]})
    f = build_cds_feature(mrna, gene, "PFX_000010", genome, cfg(), [])
    q = {x.key: x.value for x in f.qualifiers}
    assert q["product"] == "hypothetical protein"  # "protein MpX" フォールバックは廃止
    assert q["gene"] == "MpX"


def test_multi_span_join_complete():
    genome = Seq("ATGCCCTAA")  # join(1..3,7..9) -> ATG+TAA = ATGTAA -> M*
    gene, mrna = mrna_with_cds([(1, 3), (7, 9)], strand="+", phase0=0)
    f = build_cds_feature(mrna, gene, "PFX_000010", genome, cfg(), [])
    assert f.location == "join(1..3,7..9)"


def test_internal_stop_diagnostic():
    genome = Seq("ATGTAAAAA")  # ATG TAA AAA -> M*K  internal stop at codon 2
    gene, mrna = mrna_with_cds([(1, 9)], strand="+", phase0=0)
    diags = []
    build_cds_feature(mrna, gene, "PFX_000010", genome, cfg(), diags)
    assert any(d.code == "translation-internal-stop" for d in diags)


def test_translation_no_start_diagnostic():
    genome = Seq("TTGAAATAA")  # TTG in start_codons -> not 5' partial; translate -> "LK*" (no M start)
    gene, mrna = mrna_with_cds([(1, 9)], strand="+", phase0=0)
    diags = []
    build_cds_feature(mrna, gene, "PFX_000010", genome, cfg(), diags)
    assert any(d.code == "translation-no-start" for d in diags)


def test_no_cds_returns_none_and_warns():
    gene = Feature("g", "S", "gene", [Span("chr1", 1, 9, "+")], {}, [])
    mrna = Feature("g.1", "S", "mRNA", [Span("chr1", 1, 9, "+")], {}, [])
    mrna.children = []
    diags = []
    assert build_cds_feature(mrna, gene, "PFX_000010", Seq("AAA"), cfg(), diags) is None
    assert any(d.code == "no-cds" for d in diags)


def test_recoded_codon_child_avoids_internal_stop_warning():
    # CDS ATG TGA AAA TAA on + strand; TGA (4..6) is recoded via a recoded_codon child.
    genome = Seq("ATGTGAAAATAA")
    cds = Feature("c1", "S", "CDS", [Span("s", 1, 12, "+", 0)],
                  {"ID": ["c1"], "transl_table": ["11"]}, ["m1"])
    recoded = Feature("c1_recoded_1", "S", "recoded_codon", [Span("s", 4, 6, "+", 0)],
                      {"ID": ["c1_recoded_1"], "Parent": ["c1"], "codon_redefined": ["selenocysteine"]},
                      ["c1"])
    cds.children = [recoded]
    mrna = Feature("m1", "S", "mRNA", [Span("s", 1, 12, "+")], {"ID": ["m1"]}, ["g1"])
    mrna.children = [cds]
    gene = Feature("g1", "S", "gene", [Span("s", 1, 12, "+")], {"ID": ["g1"]}, [])
    cfg = MssConfig(source={"organism": "x", "mol_type": "genomic DNA"})
    diags: list = []
    feat = build_cds_feature(mrna, gene, "LT_1", genome, cfg, diags)
    assert feat is not None
    assert not any(d.code == "translation-internal-stop" for d in diags)
