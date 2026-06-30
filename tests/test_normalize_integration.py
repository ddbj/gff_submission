from pathlib import Path
import pytest
from ddbj_gff import parse
from ddbj_gff.errors import Severity
from ddbj_gff.validate import validate
from ddbj_gff.writer import write
from ddbj_gff.normalize import normalize, NormalizeConfig

FIX = Path(__file__).parent / "normalize_fixtures"
ROOT = Path(__file__).resolve().parents[1]


def test_messy_fixture_normalizes_and_validates():
    doc = parse((FIX / "messy_input.gff3").read_text())
    norm, report = normalize(doc, seq_lengths={"chr1": 5000}, config=NormalizeConfig(taxid=3702))
    # SO terms normalized to INSDC names
    types = {f.type for f in norm.features}
    assert "exon" in types and "CDS" in types          # coding_exon->exon, pseudogenic_CDS->CDS
    assert "coding_exon" not in types and "pseudogenic_CDS" not in types
    # targeted ERRORs cleared by directive completion
    errors = {d.code for d in validate(norm) if d.severity == Severity.ERROR}
    assert "missing-insdc-gff-version" not in errors
    assert "missing-species-taxid" not in errors
    assert "missing-sequence-region" not in errors
    # output is writable GFF text and round-trips through parse
    text = write(norm)
    assert "#!insdc-gff-version 1.0.0" in text
    reparsed = parse(text)
    assert reparsed.insdc_gff_version == "1.0.0"


@pytest.mark.slow
def test_rice_cp_normalize_clears_version_keeps_special_case():
    p = ROOT / "examples" / "rice_cp" / "rice_cp.gff3"
    if not p.exists():
        pytest.skip(f"missing {p}")
    norm, _ = normalize(parse(p.read_text(errors="replace")), config=NormalizeConfig(taxid=39947))
    diags = validate(norm)
    errors = {d.code for d in diags if d.severity == Severity.ERROR}
    codes = {d.code for d in diags}
    assert "missing-insdc-gff-version" not in errors          # 3B added the directive
    assert "noncanonical-special-case" in codes               # trans-splicing still flagged (deferred to 3B-full)
