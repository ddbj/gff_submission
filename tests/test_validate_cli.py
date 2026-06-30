import pytest

from ddbj_gff.validate.cli import main

VALID = (
    "##gff-version 3\n"
    "#!insdc-gff-version 1.0.0\n"
    "##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=3702\n"
    "##sequence-region chr1 1 1000\n"
    "chr1\tS\tgene\t1\t99\t.\t+\t.\tID=g;locus_tag=ABC_1\n"
    "chr1\tS\tmRNA\t1\t99\t.\t+\t.\tID=m;Parent=g\n"
    "chr1\tS\texon\t1\t99\t.\t+\t.\tID=e;Parent=m\n"
    "chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=c;Parent=m;transl_table=1\n"
)
BAD = "##gff-version 3\nchr1\tS\tgene\t1\t9\t.\t+\t.\tID=g\n"


def test_cli_valid_returns_zero(tmp_path):
    p = tmp_path / "v.gff"; p.write_text(VALID)
    assert main(["--gff", str(p)]) == 0


def test_cli_bad_returns_one_and_reports(tmp_path, capsys):
    p = tmp_path / "b.gff"; p.write_text(BAD)
    rc = main(["--gff", str(p)])
    assert rc == 1
    err = capsys.readouterr().err
    assert "ERROR" in err
    assert "missing-insdc-gff-version" in err


def test_cli_bad_severity_level_exits_with_2(tmp_path):
    # --severity CODE=bogus must produce a clean argparse error (SystemExit 2), not a traceback.
    p = tmp_path / "b.gff"; p.write_text(BAD)
    with pytest.raises(SystemExit) as exc_info:
        main(["--gff", str(p), "--severity", "gene-missing-locus-tag=bogus"])
    assert exc_info.value.code == 2


def test_cli_severity_override(tmp_path):
    # turning the only ERRORs off should make rc 0 — promote nothing, just off the blockers
    p = tmp_path / "b.gff"; p.write_text(BAD)
    rc = main(["--gff", str(p),
               "--severity", "missing-insdc-gff-version=off",
               "--severity", "missing-species-taxid=off",
               "--severity", "missing-sequence-region=off",
               "--severity", "undefined-seqid=off"])
    assert rc == 0
