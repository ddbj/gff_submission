from ddbj_gff.normalize.config import NormalizeConfig


def test_merge_config_defaults_off():
    c = NormalizeConfig()
    assert c.merge_overlapping_loci is False
    assert c.merge_overlap_min_fraction == 0.0


from ddbj_gff import parse
from ddbj_gff.normalize.normalize import normalize
from ddbj_gff.normalize.config import NormalizeConfig

HDR = "##gff-version 3\n##sequence-region c 1 100000\n"


def _gff(rows):
    return HDR + "".join(rows)


def _gene(gid, s, e, strand="+"):
    m = f"{gid}.m"
    return (f"c\tx\tgene\t{s}\t{e}\t.\t{strand}\t.\tID={gid}\n"
            f"c\tx\tmRNA\t{s}\t{e}\t.\t{strand}\t.\tID={m};Parent={gid}\n"
            f"c\tx\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID={gid}.c;Parent={m}\n")


def _norm(gff, **cfg):
    doc = parse(gff)
    work, _ = normalize(doc, config=NormalizeConfig(taxid=1, **cfg))
    return work


def _genes(doc):
    return [f for f in doc.features if f.type == "gene"]


def _mrnas(doc):
    return [f for f in doc.features if f.type == "mRNA"]


def test_merge_two_overlapping_same_strand():
    gff = _gff([_gene("gA", 100, 500), _gene("gB", 300, 700)])   # overlap 300..500
    doc = _norm(gff, merge_overlapping_loci=True)
    genes = _genes(doc)
    assert len(genes) == 1 and genes[0].id == "gA"               # rep = lowest start
    assert genes[0].spans[0].start == 100 and genes[0].spans[0].end == 700   # union
    mrnas = _mrnas(doc)
    assert len(mrnas) == 2 and all(m.parent_ids == ["gA"] for m in mrnas)     # both under gA
    assert "gB" not in doc.feature_index                          # merged-away gene removed
    codes = {d.code for d in validate_ok(doc)}
    assert "dangling-parent" not in codes


def validate_ok(doc):
    from ddbj_gff.validate import validate
    return validate(doc)


def test_flag_off_no_change():
    gff = _gff([_gene("gA", 100, 500), _gene("gB", 300, 700)])
    doc = _norm(gff)                                             # flag default off
    assert len(_genes(doc)) == 2 and len(_mrnas(doc)) == 2


def test_opposite_strand_not_merged():
    gff = _gff([_gene("gA", 100, 500, "+"), _gene("gB", 300, 700, "-")])
    doc = _norm(gff, merge_overlapping_loci=True)
    assert len(_genes(doc)) == 2                                 # antisense stays separate


def test_transitive_chain_merged():
    gff = _gff([_gene("gA", 100, 500), _gene("gB", 400, 800), _gene("gC", 700, 1000)])
    doc = _norm(gff, merge_overlapping_loci=True)                # A~B, B~C, A!~C
    genes = _genes(doc)
    assert len(genes) == 1 and genes[0].id == "gA"
    assert genes[0].spans[0].end == 1000
    assert len(_mrnas(doc)) == 3


def test_threshold_below_not_merged():
    # gA 100..500 (len401), gB 300..700 (len401); overlap 201; 201/401 = 0.50
    gff = _gff([_gene("gA", 100, 500), _gene("gB", 300, 700)])
    doc = _norm(gff, merge_overlapping_loci=True, merge_overlap_min_fraction=0.9)
    assert len(_genes(doc)) == 2                                 # 0.50 < 0.90 -> not merged


def test_trans_spliced_exempt():
    # a trans-spliced CDS (two parts) whose mRNA extent 100..900 overlaps a normal gene gN 200..600
    trans = ("c\tx\tgene\t100\t900\t.\t+\t.\tID=gT\n"
             "c\tx\tmRNA\t100\t900\t.\t+\t.\tID=gT.m;Parent=gT\n"
             "c\tx\tCDS\t100\t200\t.\t+\t0\tID=gT.c;Parent=gT.m;exception=trans-splicing;part=1\n"
             "c\tx\tCDS\t800\t900\t.\t+\t0\tID=gT.c;Parent=gT.m;exception=trans-splicing;part=2\n")
    gff = _gff([trans, _gene("gN", 200, 600)])
    doc = _norm(gff, merge_overlapping_loci=True)
    ids = {g.id for g in _genes(doc)}
    assert ids == {"gT", "gN"}                                  # trans-spliced gT exempt; gN untouched
