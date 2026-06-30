from ddbj_gff import parse
from ddbj_gff.errors import Severity
from ddbj_gff.validate import validate

GFF_BAD = (
    "##gff-version 3\n"
    "chr1\tS\tgene\t1\t9\t.\t+\t.\tID=g\n"   # no #!insdc-gff-version, no ##species, no ##sequence-region, gene no locus_tag
)


def test_validate_runs_all_and_finds_errors():
    diags = validate(parse(GFF_BAD))
    codes = {d.code for d in diags}
    assert "missing-insdc-gff-version" in codes
    assert "missing-sequence-region" in codes
    assert "gene-missing-locus-tag" in codes


def test_severity_override_off_and_promote():
    base = {d.code for d in validate(parse(GFF_BAD))}
    assert "gene-missing-locus-tag" in base
    # off removes it
    off = validate(parse(GFF_BAD), severity_overrides={"gene-missing-locus-tag": "off"})
    assert "gene-missing-locus-tag" not in {d.code for d in off}
    # promote to error
    promoted = validate(parse(GFF_BAD), severity_overrides={"gene-missing-locus-tag": "error"})
    g = [d for d in promoted if d.code == "gene-missing-locus-tag"][0]
    assert g.severity == Severity.ERROR


def test_validate_sorted_and_does_not_mutate_doc():
    doc = parse(GFF_BAD)
    before = len(doc.diagnostics)
    diags = validate(doc)
    keys = [(d.line_no if d.line_no is not None else -1, d.code) for d in diags]
    assert keys == sorted(keys)
    assert len(doc.diagnostics) == before  # detect-only: validate must not append to doc.diagnostics
