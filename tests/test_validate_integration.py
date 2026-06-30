import gzip
from pathlib import Path
import pytest
from ddbj_gff import parse
from ddbj_gff.errors import Severity
from ddbj_gff.validate import validate

FIX = Path(__file__).parent / "validate_fixtures"
ROOT = Path(__file__).resolve().parents[1]


def test_valid_insdc_fixture_has_no_errors():
    diags = validate(parse((FIX / "valid_insdc.gff3").read_text()))
    errors = [d for d in diags if d.severity == Severity.ERROR]
    assert errors == [], f"unexpected errors: {[(d.code, d.message) for d in errors]}"


@pytest.mark.slow
def test_rice_cp_flags_insdc_violations():
    p = ROOT / "examples" / "rice_cp" / "rice_cp.gff3"
    if not p.exists():
        pytest.skip(f"missing {p}")
    codes = {d.code for d in validate(parse(p.read_text(errors="replace")))}
    # NCBI-style file lacks the INSDC version directive
    assert "missing-insdc-gff-version" in codes
    # rps12 uses exception=trans-splicing + part= (non-canonical)
    assert "noncanonical-special-case" in codes


@pytest.mark.slow
def test_ecoli_flags_transl_except_noncanonical():
    p = ROOT / "examples" / "ecoli" / "GCF_000005845.2_ASM584v2_genomic.gff.gz"
    if not p.exists():
        pytest.skip(f"missing {p}")
    text = gzip.decompress(p.read_bytes()).decode("ascii", errors="replace")
    codes = {d.code for d in validate(parse(text))}
    assert "missing-insdc-gff-version" in codes
    assert "noncanonical-special-case" in codes   # transl_except attribute present
