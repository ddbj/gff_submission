from ddbj_gff.model import Feature, Span, GffDocument
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.normalize.passes import NormalizeContext, pass_so_terms
from ddbj_gff.validate.vocab import load_vocab


def _ctx():
    return NormalizeContext(vocab=load_vocab(), seq_lengths=None, config=NormalizeConfig())


def _doc(*feats):
    return GffDocument(features=list(feats))


# --- collapse to a core whitelist type (gene/mRNA/CDS/exon/intron) ---

def test_pseudogenic_cds_collapses_to_cds_with_pseudo_flag():
    f = Feature("c", "S", "pseudogenic_CDS", [Span("chr1", 1, 9, "+", 0)], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "CDS"
    assert f.attributes.get("pseudo") == ["true"]
    assert any(c.action == "rename-type" for c in changes)
    assert any(c.action == "add-qualifier" for c in changes)


def test_processed_pseudogene_collapses_to_gene_with_qualifier():
    f = Feature("g", "S", "processed_pseudogene", [Span("chr1", 1, 9, "+")], {}, [])
    pass_so_terms(_doc(f), _ctx())
    assert f.type == "gene"
    assert f.attributes.get("pseudogene") == ["processed"]


def test_coding_exon_collapses_to_exon():
    f = Feature("e", "S", "coding_exon", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "exon"
    assert any(c.action == "rename-type" for c in changes)


def test_spliceosomal_intron_collapses_to_intron():
    f = Feature("i", "S", "spliceosomal_intron", [Span("chr1", 1, 9, "+")], {}, [])
    pass_so_terms(_doc(f), _ctx())
    assert f.type == "intron"


def test_same_name_is_noop():
    f = Feature("c", "S", "CDS", [Span("chr1", 1, 9, "+", 0)], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "CDS"
    assert changes == []


def test_existing_attribute_not_clobbered():
    f = Feature("c", "S", "pseudogenic_CDS", [Span("chr1", 1, 9, "+", 0)], {"pseudo": ["existing"]}, [])
    pass_so_terms(_doc(f), _ctx())
    assert f.attributes["pseudo"] == ["existing"]


# --- non-core targets are LEFT for Phase 2 (Finding A/B fix) ---

def test_mirna_left_for_phase2():
    f = Feature("r", "S", "miRNA", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "miRNA"        # unchanged: Phase 2 maps it to ncRNA[ncRNA_class=miRNA]
    assert changes == []


def test_pre_mirna_left_for_phase2():
    f = Feature("r", "S", "pre_miRNA", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "pre_miRNA"    # unchanged: Phase 2 maps it to precursor_RNA
    assert changes == []


def test_five_prime_utr_left_alone():
    f = Feature("u", "S", "five_prime_UTR", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "five_prime_UTR"
    assert changes == []


def test_non_core_insdc_target_left_alone():
    # binding_site -> misc_binding, mobile_genetic_element -> mobile_element:
    # targets not in the core whitelist -> left unchanged, no fabricated qualifier, no Change
    for t in ("binding_site", "mobile_genetic_element"):
        f = Feature("x", "S", t, [Span("chr1", 1, 9, "+")], {}, [])
        changes = pass_so_terms(_doc(f), _ctx())
        assert f.type == t
        assert f.attributes == {}
        assert changes == []


def test_unmapped_type_reported_unchanged():
    f = Feature("x", "S", "totally_made_up_type", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "totally_made_up_type"
    assert any(c.action == "unmapped-type" for c in changes)
