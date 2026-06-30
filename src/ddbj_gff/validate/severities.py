from __future__ import annotations

from ..errors import Diagnostic, Severity

DEFAULT_SEVERITIES: dict[str, Severity] = {
    "missing-gff-version": Severity.ERROR,
    "missing-insdc-gff-version": Severity.ERROR,
    "missing-species-taxid": Severity.ERROR,
    "missing-sequence-region": Severity.ERROR,
    "duplicate-sequence-region": Severity.ERROR,
    "non-ascii": Severity.ERROR,
    "undefined-seqid": Severity.ERROR,
    "feature-outside-region": Severity.ERROR,
    "start-gt-end": Severity.ERROR,
    "feature-type-not-insdc": Severity.WARNING,
    "multiple-parents": Severity.ERROR,
    "dangling-parent": Severity.ERROR,
    "cds-missing-transl-table": Severity.ERROR,
    "cds-invalid-phase": Severity.ERROR,
    "gene-missing-locus-tag": Severity.WARNING,
    "dbxref-unknown-dbtag": Severity.WARNING,
    "noncanonical-special-case": Severity.INFO,
}

_LEVELS = {"error": Severity.ERROR, "warning": Severity.WARNING, "info": Severity.INFO}


def make_diagnostic(code: str, message: str, line_no: int | None = None) -> Diagnostic:
    return Diagnostic(DEFAULT_SEVERITIES.get(code, Severity.WARNING), line_no, code, message)


def resolve_level(name: str) -> Severity | None:
    key = name.lower()
    if key == "off":
        return None
    if key in _LEVELS:
        return _LEVELS[key]
    raise ValueError(f"invalid severity level: {name!r} (use error/warning/info/off)")
