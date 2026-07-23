from ddbj_gff.repair.context import RepairContext
from ddbj_gff.repair.report import Candidate, candidates_to_json, render_candidates
from ddbj_gff.repair.registry import Operation, register, get_operation, list_operations, REGISTRY


def test_context_defaults():
    ctx = RepairContext()
    assert ctx.sequences is None
    assert ctx.transl_table == 1


def test_candidate_json_and_render():
    c = Candidate("op-x", "gene1", "chr1", "would do X", payload={"five": True})
    import json
    parsed = json.loads(candidates_to_json([c]))
    assert parsed == [{"operation": "op-x", "feature_id": "gene1",
                       "seqid": "chr1", "detail": "would do X",
                       "payload": {"five": True}}]
    text = render_candidates([c])
    assert "op-x" in text and "gene1" in text


def test_register_and_lookup():
    op = Operation("op-test", "test op", requires_sequence=False,
                   detect=lambda doc, ctx: [], apply=lambda doc, ctx, sel: [])
    register(op)
    assert get_operation("op-test") is op
    assert "op-test" in {o.name for o in list_operations()}
    del REGISTRY["op-test"]
