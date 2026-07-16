from ddbj_gff.normalize.config import NormalizeConfig, load_normalize_config


def test_reparent_config_default_on():
    assert NormalizeConfig().reparent_gene_children is True


def test_reparent_config_loader_reads_flag(tmp_path):
    p = tmp_path / "n.toml"
    p.write_text("[normalize]\nreparent_gene_children = false\n")
    assert load_normalize_config(str(p)).reparent_gene_children is False


def test_reparent_config_loader_defaults_on(tmp_path):
    p = tmp_path / "n.toml"
    p.write_text("[normalize]\ntaxid = 1\n")
    assert load_normalize_config(str(p)).reparent_gene_children is True


from ddbj_gff import parse
from ddbj_gff.normalize.normalize import normalize

HDR = "##gff-version 3\n##sequence-region c 1 100000\n"


def _norm(gff, **cfg):
    doc = parse(gff)
    work, report = normalize(doc, config=NormalizeConfig(taxid=1, **cfg))
    return work, report


def _feat(doc, fid):
    return doc.feature_index.get(fid)


def _validate(doc):
    from ddbj_gff.validate import validate
    return validate(doc)


def test_misparented_cds_exon_reparented_to_mrna():
    gff = HDR + (
        "c\tx\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "c\tx\tmRNA\t100\t500\t.\t+\t.\tID=g1.1;Parent=g1\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g1.cds1;Parent=g1\n"
        "c\tx\texon\t100\t200\t.\t+\t.\tID=g1.ex1;Parent=g1\n"
    )
    doc, report = _norm(gff)
    gene, mrna = _feat(doc, "g1"), _feat(doc, "g1.1")
    cds, exon = _feat(doc, "g1.cds1"), _feat(doc, "g1.ex1")
    assert cds.parent_ids == ["g1.1"] and exon.parent_ids == ["g1.1"]
    assert cds.attributes["Parent"] == ["g1.1"]
    assert {c.id for c in mrna.children} == {"g1.cds1", "g1.ex1"}
    assert {c.id for c in gene.children} == {"g1.1"}
    assert any(ch.action == "reparent-to-mrna" for ch in report.applied)
    assert "dangling-parent" not in {d.code for d in _validate(doc)}


def test_wellformed_three_level_is_noop():
    gff = HDR + (
        "c\tx\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "c\tx\tmRNA\t100\t500\t.\t+\t.\tID=g1.1;Parent=g1\n"
        "c\tx\tCDS\t100\t500\t.\t+\t0\tID=g1.cds1;Parent=g1.1\n"
        "c\tx\texon\t100\t500\t.\t+\t.\tID=g1.ex1;Parent=g1.1\n"
    )
    doc, report = _norm(gff)
    assert not any(ch.action == "reparent-to-mrna" for ch in report.applied)
    assert _feat(doc, "g1.cds1").parent_ids == ["g1.1"]
    assert {c.id for c in _feat(doc, "g1").children} == {"g1.1"}


def test_two_mrnas_gene_level_cds_not_reparented():
    gff = HDR + (
        "c\tx\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "c\tx\tmRNA\t100\t500\t.\t+\t.\tID=g1.1;Parent=g1\n"
        "c\tx\tmRNA\t100\t500\t.\t+\t.\tID=g1.2;Parent=g1\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g1.cds1;Parent=g1\n"
    )
    doc, report = _norm(gff)
    assert _feat(doc, "g1.cds1").parent_ids == ["g1"]
    assert not any(ch.action == "reparent-to-mrna" for ch in report.applied)
    assert any(ch.action == "needs-manual" for ch in report.unresolved)


def test_no_transcript_left_for_wrap():
    gff = HDR + (
        "c\tx\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "c\tx\tCDS\t100\t500\t.\t+\t0\tID=g1.cds1;Parent=g1\n"
        "c\tx\texon\t100\t500\t.\t+\t.\tID=g1.ex1;Parent=g1\n"
    )
    doc, report = _norm(gff)
    mrnas = [f for f in doc.features if f.type == "mRNA"]
    assert len(mrnas) == 1
    assert {c.id for c in mrnas[0].children} == {"g1.cds1", "g1.ex1"}
    assert any(ch.action == "add-child-feature" for ch in report.applied)
    assert not any(ch.action == "reparent-to-mrna" for ch in report.applied)


def test_sole_trna_transcript_not_reparented():
    gff = HDR + (
        "c\tx\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "c\tx\ttRNA\t100\t500\t.\t+\t.\tID=g1.t;Parent=g1\n"
        "c\tx\texon\t100\t500\t.\t+\t.\tID=g1.ex1;Parent=g1\n"
    )
    doc, report = _norm(gff)
    assert _feat(doc, "g1.ex1").parent_ids == ["g1"]
    assert not any(ch.action == "reparent-to-mrna" for ch in report.applied)
    assert any(ch.action == "needs-manual" for ch in report.unresolved)


def test_multiexon_misparent_reparented_span_covers():
    gff = HDR + (
        "c\tx\tgene\t100\t900\t.\t+\t.\tID=g1\n"
        "c\tx\tmRNA\t100\t900\t.\t+\t.\tID=g1.1;Parent=g1\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g1.c1;Parent=g1\n"
        "c\tx\tCDS\t800\t900\t.\t+\t2\tID=g1.c2;Parent=g1\n"
        "c\tx\texon\t100\t200\t.\t+\t.\tID=g1.e1;Parent=g1\n"
        "c\tx\texon\t800\t900\t.\t+\t.\tID=g1.e2;Parent=g1\n"
    )
    doc, report = _norm(gff)
    mrna = _feat(doc, "g1.1")
    assert {c.id for c in mrna.children} == {"g1.c1", "g1.c2", "g1.e1", "g1.e2"}
    assert mrna.spans[0].start == 100 and mrna.spans[0].end == 900
    for cid in ("g1.c1", "g1.c2", "g1.e1", "g1.e2"):
        assert _feat(doc, cid).parent_ids == ["g1.1"]


def test_flag_off_no_reparent():
    gff = HDR + (
        "c\tx\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "c\tx\tmRNA\t100\t500\t.\t+\t.\tID=g1.1;Parent=g1\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g1.cds1;Parent=g1\n"
    )
    doc, report = _norm(gff, reparent_gene_children=False)
    assert _feat(doc, "g1.cds1").parent_ids == ["g1"]
    assert not any(ch.action == "reparent-to-mrna" for ch in report.applied)


def test_trans_spliced_mrna_gene_level_intron_not_reparented():
    # An mRNA flagged trans-spliced (exception=trans-splicing) with a gene-level
    # intron sibling: the guard must skip the gene entirely, leaving the intron
    # parented to the gene and the mRNA's own children untouched. This protects
    # the trans-splicing structure flatfile2gff builds (introns intentionally
    # gene-level) from the pass's single-Span recompute.
    gff = HDR + (
        "c\tx\tgene\t100\t900\t.\t+\t.\tID=g1\n"
        "c\tx\tmRNA\t100\t900\t.\t+\t.\tID=g1.1;Parent=g1;exception=trans-splicing\n"
        "c\tx\tCDS\t100\t200\t.\t+\t0\tID=g1.c1;Parent=g1.1\n"
        "c\tx\tintron\t201\t799\t.\t+\t.\tID=g1.in1;Parent=g1\n"
    )
    doc, report = _norm(gff)
    assert _feat(doc, "g1.in1").parent_ids == ["g1"]              # intron NOT reparented
    assert not any(ch.action == "reparent-to-mrna" for ch in report.applied)
    assert {c.id for c in _feat(doc, "g1.1").children} == {"g1.c1"}  # mRNA keeps only its own CDS
