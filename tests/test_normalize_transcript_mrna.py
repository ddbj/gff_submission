from ddbj_gff import parse
from ddbj_gff.normalize.passes import pass_coerce_transcript_to_mrna, NormalizeContext
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.validate.vocab import load_vocab

GFF = """##gff-version 3
c1\tS\tgene\t1\t9\t.\t+\t.\tID=g1
c1\tS\ttranscript\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1
c1\tS\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=g1.t1
c1\tS\tgene\t20\t30\t.\t+\t.\tID=g2
c1\tS\ttranscript\t20\t30\t.\t+\t.\tID=nc1;Parent=g2
"""


def test_coding_transcript_becomes_mrna_noncoding_untouched():
    doc = parse(GFF)
    ctx = NormalizeContext(vocab=load_vocab(), seq_lengths=None, config=NormalizeConfig())
    pass_coerce_transcript_to_mrna(doc, ctx)
    types = {f.id: f.type for f in doc.features}
    assert types["g1.t1"] == "mRNA"     # CDS を持つ transcript
    assert types["nc1"] == "transcript"  # CDS なしは不変
