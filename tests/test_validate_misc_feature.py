from ddbj_gff import parse
from ddbj_gff.validate import validate

HDR = "##gff-version 3\n##sequence-region s 1 100000\n"
GFF = HDR + "s\tsrc\tmisc_feature\t1\t12\t.\t+\t.\tID=x;Note=nonfunctional\n"


def test_misc_feature_is_accepted_type():
    diags = validate(parse(GFF))
    assert not any(d.code == "feature-type-not-insdc" for d in diags)
