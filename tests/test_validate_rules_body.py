from ddbj_gff.model import Feature, Span, GffDocument
from ddbj_gff.validate.vocab import Vocab
from ddbj_gff.validate import rules

V = Vocab(frozenset({"gene", "mRNA", "CDS", "exon"}), {"CDS": "CDS"}, frozenset({"GenBank", "GeneID"}))


def codes(diags):
    return {d.code for d in diags}


def test_feature_type_not_insdc():
    f = Feature("a", "S", "weird_type", [Span("c", 1, 9, "+")], {}, [])
    assert "feature-type-not-insdc" in codes(rules.rule_feature_type(GffDocument(features=[f]), V))
    f2 = Feature("b", "S", "gene", [Span("c", 1, 9, "+")], {}, [])
    assert "feature-type-not-insdc" not in codes(rules.rule_feature_type(GffDocument(features=[f2]), V))


def test_multiple_parents_and_dangling():
    parent = Feature("g", "S", "gene", [Span("c", 1, 99, "+")], {}, [])
    multi = Feature("m", "S", "mRNA", [Span("c", 1, 99, "+")], {}, ["g", "h"])   # 2 parents
    dangling = Feature("m2", "S", "mRNA", [Span("c", 1, 99, "+")], {}, ["ghost"])
    doc = GffDocument(features=[parent, multi, dangling],
                      feature_index={"g": parent, "m": multi, "m2": dangling})
    c = codes(rules.rule_parents(doc, V))
    assert "multiple-parents" in c
    assert "dangling-parent" in c


def test_cds_missing_transl_table():
    cds = Feature("c1", "S", "CDS", [Span("c", 1, 9, "+", 0)], {}, [])
    doc = GffDocument(features=[cds])  # no transl_table attr, no file-level #!transl_table
    assert "cds-missing-transl-table" in codes(rules.rule_cds(doc, V))


def test_cds_transl_table_satisfied_by_file_pragma():
    from ddbj_gff.model import Directive
    cds = Feature("c1", "S", "CDS", [Span("c", 1, 9, "+", 0)], {}, [])
    doc = GffDocument(directives=[Directive("#!transl_table primary:1", "transl_table", {"primary": 1})],
                      features=[cds])
    assert "cds-missing-transl-table" not in codes(rules.rule_cds(doc, V))


def test_cds_invalid_phase():
    cds = Feature("c1", "S", "CDS", [Span("c", 1, 9, "+", None)], {"transl_table": ["11"]}, [])  # phase None
    assert "cds-invalid-phase" in codes(rules.rule_cds(doc=GffDocument(features=[cds]), vocab=V))


def test_gene_missing_locus_tag():
    g = Feature("g", "S", "gene", [Span("c", 1, 9, "+")], {}, [])
    assert "gene-missing-locus-tag" in codes(rules.rule_gene_locus_tag(GffDocument(features=[g]), V))
    g2 = Feature("g2", "S", "gene", [Span("c", 1, 9, "+")], {"locus_tag": ["X_1"]}, [])
    assert "gene-missing-locus-tag" not in codes(rules.rule_gene_locus_tag(GffDocument(features=[g2]), V))


def test_dbxref_unknown_dbtag():
    f = Feature("a", "S", "gene", [Span("c", 1, 9, "+")], {"Dbxref": ["GenBank:X1", "WeirdDB:9"]}, [])
    c = codes(rules.rule_dbxref(GffDocument(features=[f]), V))
    assert "dbxref-unknown-dbtag" in c   # WeirdDB unknown


def test_special_case_detection():
    ts = Feature("g", "S", "gene", [Span("c", 1, 9, "-")], {"exception": ["trans-splicing"]}, [])
    te = Feature("c1", "S", "CDS", [Span("c", 1, 9, "+", 0)],
                 {"transl_except": ["(pos:1..3,aa:Sec)"], "transl_table": ["11"]}, [])
    c = codes(rules.rule_special_case(GffDocument(features=[ts, te]), V))
    assert "noncanonical-special-case" in c


def test_all_rules_list_complete():
    assert len(rules.ALL_RULES) == 10
