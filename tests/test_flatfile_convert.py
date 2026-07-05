from Bio import SeqIO
from ddbj_gff.flatfile.convert import flatfile_to_gff
from ddbj_gff.writer import write
from ddbj_gff.validate import validate

FIX = "tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk"


def test_flatfile_to_gff_validates_and_has_hierarchy():
    rec = SeqIO.read(FIX, "genbank")
    doc = flatfile_to_gff(rec)
    types = [f.type for f in doc.features]
    assert types.count("gene") == 2 and types.count("mRNA") == 3
    assert types.count("exon") == 7 and types.count("CDS") == 3
    assert any(f.type == "region" for f in doc.features)
    # directives added by normalize
    text = write(doc)
    assert "##sequence-region" in text
    assert "id=55188" in text            # ##species taxid
    assert "##gff-version 3" in text
    # canonical GFF validates with no ERROR-level diagnostics
    diags = validate(doc)
    errors = [d for d in diags if getattr(d, "severity", None) and d.severity.name == "ERROR"]
    assert errors == [], f"unexpected validate errors: {[d.code for d in errors]}"
