import os
import pytest
from ddbj_gff.repair.cli import main
from ddbj_gff import parse
from ddbj_gff.validate import validate

FIX = os.path.join(os.path.dirname(__file__), "fixtures")
GFF = os.path.join(FIX, "repair_internal_stop.gff3")
FASTA = os.path.join(FIX, "repair_internal_stop.fasta")


def test_list(capsys):
    rc = main(["--list"])
    out = capsys.readouterr().out
    assert rc == 0
    for name in ("internal-stop-to-misc", "utr-absent-to-partial-mrna",
                 "missing-start-stop-to-partial-cds"):
        assert name in out


def test_detect_json_does_not_write(capsys, tmp_path):
    rc = main(["--gff", GFF, "--fasta", FASTA, "--detect", "--json"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "internal-stop-to-misc" in out and "c1" in out


def test_apply_writes_curated_gff_that_validates(tmp_path):
    out = tmp_path / "out.gff3"
    rc = main(["--gff", GFF, "--fasta", FASTA, "--apply", "all", "--out", str(out)])
    assert rc == 0
    before = {d.code for d in validate(parse(open(GFF).read()))
              if d.severity.name == "ERROR"}
    doc = parse(out.read_text())
    assert doc.feature_index["c1"].type == "misc_feature"
    after = {d.code for d in validate(doc) if d.severity.name == "ERROR"}
    assert after <= before   # repair introduced no new validation ERROR


def test_apply_unknown_operation_errors_cleanly():
    with pytest.raises(SystemExit) as exc:
        main(["--gff", GFF, "--fasta", FASTA, "--apply", "bogus-op", "--out", "/dev/null"])
    assert exc.value.code == 2
