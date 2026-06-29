import pytest

from ddbj_gff.errors import Diagnostic, GffParseError, Severity


def test_severity_members():
    assert Severity.ERROR.value == "ERROR"
    assert {s.name for s in Severity} == {"ERROR", "WARNING", "INFO"}


def test_diagnostic_is_frozen_and_equal():
    d1 = Diagnostic(Severity.WARNING, 12, "dangling-parent", "Parent x not found")
    d2 = Diagnostic(Severity.WARNING, 12, "dangling-parent", "Parent x not found")
    assert d1 == d2
    with pytest.raises(Exception):
        d1.code = "other"  # frozen dataclass


def test_parse_error_carries_diagnostic():
    d = Diagnostic(Severity.ERROR, 3, "col-count", "expected 9 columns, got 7")
    err = GffParseError(d)
    assert err.diagnostic is d
    assert "col-count" in str(err)
    assert "line 3" in str(err)
