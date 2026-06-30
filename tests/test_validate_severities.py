import pytest
from ddbj_gff.errors import Severity
from ddbj_gff.validate.severities import DEFAULT_SEVERITIES, make_diagnostic, resolve_level


def test_default_severities_cover_known_codes():
    for code in ("missing-insdc-gff-version", "feature-type-not-insdc",
                 "noncanonical-special-case", "gene-missing-locus-tag"):
        assert code in DEFAULT_SEVERITIES


def test_make_diagnostic_uses_default_severity():
    d = make_diagnostic("missing-insdc-gff-version", "no insdc version")
    assert d.severity == Severity.ERROR
    assert d.code == "missing-insdc-gff-version"
    d2 = make_diagnostic("feature-type-not-insdc", "x")
    assert d2.severity == Severity.WARNING
    d3 = make_diagnostic("noncanonical-special-case", "x")
    assert d3.severity == Severity.INFO


def test_resolve_level():
    assert resolve_level("error") == Severity.ERROR
    assert resolve_level("warning") == Severity.WARNING
    assert resolve_level("info") == Severity.INFO
    assert resolve_level("off") is None
    with pytest.raises(ValueError):
        resolve_level("bogus")
