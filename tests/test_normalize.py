from ddbj_gff import parse
from ddbj_gff.errors import Severity
from ddbj_gff.validate import validate
from ddbj_gff.normalize import normalize, NormalizeConfig

GFF_MESSY = (
    "##gff-version 3\n"
    "chr1\tS\tgene\t100\t900\t.\t+\t.\tID=g;locus_tag=ABC_1\n"
    "chr1\tS\tCDS\t130\t870\t.\t+\t0\tID=c;Parent=g\n"
)  # missing insdc-gff-version, species, sequence-region, transl_table


def test_normalize_clears_targeted_errors():
    norm, report = normalize(parse(GFF_MESSY),
                             seq_lengths={"chr1": 10000}, config=NormalizeConfig(taxid=3702))
    codes = {d.code for d in validate(norm) if d.severity == Severity.ERROR}
    assert "missing-insdc-gff-version" not in codes
    assert "missing-species-taxid" not in codes
    assert "missing-sequence-region" not in codes
    assert "cds-missing-transl-table" not in codes


def test_normalize_does_not_mutate_input():
    doc = parse(GFF_MESSY)
    before = len(doc.directives)
    normalize(doc, config=NormalizeConfig(taxid=3702))
    assert len(doc.directives) == before  # input untouched (works on a copy)


def test_report_separates_applied_and_unresolved():
    # no taxid -> species unresolved; an unmapped feature -> unresolved
    g = "##gff-version 3\nchr1\tS\tmade_up\t1\t9\t.\t+\t.\tID=x\n"
    _, report = normalize(parse(g))
    assert any(c.action == "no-taxid" for c in report.unresolved)
    assert any(c.action == "unmapped-type" for c in report.unresolved)
    assert all(c.action in ("add-directive",) for c in report.applied)


def test_idempotent():
    norm1, _ = normalize(parse(GFF_MESSY), seq_lengths={"chr1": 10000}, config=NormalizeConfig(taxid=3702))
    norm2, report2 = normalize(norm1, seq_lengths={"chr1": 10000}, config=NormalizeConfig(taxid=3702))
    assert len(norm2.directives) == len(norm1.directives)
    assert report2.applied == []  # nothing left to change
