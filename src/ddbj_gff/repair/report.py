from __future__ import annotations

import json
from dataclasses import dataclass, field


@dataclass
class Candidate:
    operation: str
    feature_id: str | None
    seqid: str
    detail: str
    payload: dict = field(default_factory=dict)


def candidates_to_json(cands: list[Candidate]) -> str:
    return json.dumps([
        {"operation": c.operation, "feature_id": c.feature_id, "seqid": c.seqid,
         "detail": c.detail, "payload": c.payload}
        for c in cands
    ])


def render_candidates(cands: list[Candidate]) -> str:
    lines = [f"repair: {len(cands)} candidate(s)"]
    for c in cands:
        lines.append(f"  [{c.operation}] {c.feature_id} ({c.seqid}): {c.detail}")
    return "\n".join(lines) + "\n"


def render_changes(changes: list) -> str:
    lines = [f"repair: {len(changes)} change(s) applied"]
    for c in changes:
        lines.append(f"  [{c.action}] {c.target}: {c.message}")
    return "\n".join(lines) + "\n"
