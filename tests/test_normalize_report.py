from ddbj_gff.normalize.report import Change, NormalizationReport


def test_report_render_and_counts():
    r = NormalizationReport(
        applied=[Change("rename-type", "c1", "pseudogenic_CDS -> CDS")],
        unresolved=[Change("no-taxid", "species", "no taxid provided")],
    )
    text = r.render()
    assert "1 applied, 1 need attention" in text
    assert "rename-type" in text
    assert "no-taxid" in text


def test_report_defaults_empty():
    r = NormalizationReport()
    assert r.applied == [] and r.unresolved == []
    assert "0 applied, 0 need attention" in r.render()
