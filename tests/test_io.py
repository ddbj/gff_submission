import gzip
from ddbj_gff.io import open_text


def test_open_text_reads_gzip(tmp_path):
    p = tmp_path / "x.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">a\nACGT\n")
    with open_text(str(p)) as fh:
        assert fh.read() == ">a\nACGT\n"


def test_open_text_reads_plain(tmp_path):
    p = tmp_path / "y.txt"
    p.write_text("hello", encoding="utf-8")
    with open_text(str(p)) as fh:
        assert fh.read() == "hello"
