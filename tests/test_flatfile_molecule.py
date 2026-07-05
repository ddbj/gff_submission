from Bio import SeqIO
from ddbj_gff.flatfile.molecule import detect_molecule

FIX = "tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk"


def test_detect_nuclear_plant():
    rec = SeqIO.read(FIX, "genbank")
    m = detect_molecule(rec)
    assert m.taxid == 55188
    assert m.organism == "Citrus unshiu"
    assert m.division == "PLN"
    assert m.topology == "linear"
    assert m.compartment == "nuclear"     # source has no /organelle
    assert m.hierarchy == "three_level"
    assert m.transl_table == 1            # from CDS /transl_table
