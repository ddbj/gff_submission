from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from ddbj_gff.mss.translate import translate_cds_with_transl_except


def test_translate_applies_selenocysteine_exception():
    # ATG TGA AAA TAA : codon2 TGA is a stop under table 11, recoded to Sec (U) via transl_except
    parent = Seq("ATGTGAAAATAA")
    feat = SeqFeature(SimpleLocation(0, 12, strand=1), type="CDS",
                      qualifiers={"transl_table": ["11"], "codon_start": ["1"],
                                  "transl_except": ["(pos:4..6,aa:Sec)"]})
    protein = str(translate_cds_with_transl_except(feat, parent))
    assert "U" in protein          # TGA recoded to selenocysteine
    assert "*" not in protein      # no internal stop after applying the exception
    assert protein.startswith("M")
