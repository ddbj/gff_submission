import pytest
from ddbj_gff.normalize.cli import main

GFF = "##gff-version 3\nchr1\tS\tgene\t1\t9\t.\t+\t.\tID=g;locus_tag=X_1\n"


def test_cli_normalizes_to_stdout(tmp_path, capsys):
    p = tmp_path / "in.gff"
    p.write_text(GFF)
    rc = main(["--gff", str(p), "--taxid", "3702"])
    assert rc == 0
    out = capsys.readouterr()
    assert "#!insdc-gff-version 1.0.0" in out.out
    assert "##species" in out.out and "id=3702" in out.out
    assert "##sequence-region chr1 1" in out.out
    assert "[applied]" in out.err  # report to stderr


def test_cli_writes_out_file(tmp_path):
    p = tmp_path / "in.gff"
    p.write_text(GFF)
    outp = tmp_path / "out.gff"
    rc = main(["--gff", str(p), "--taxid", "3702", "--out", str(outp)])
    assert rc == 0
    assert "#!insdc-gff-version" in outp.read_text()


def test_cli_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        main(["--gff", str(tmp_path / "nope.gff")])
