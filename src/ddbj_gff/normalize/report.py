from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Change:
    action: str    # add-directive | rename-type | add-qualifier | approx-region
                   # | unmapped-type | needs-manual | no-taxid
    target: str    # feature id / seqid / directive 名
    message: str


@dataclass
class NormalizationReport:
    applied: list = field(default_factory=list)
    unresolved: list = field(default_factory=list)

    def render(self) -> str:
        lines = [f"normalization: {len(self.applied)} applied, {len(self.unresolved)} need attention"]
        for c in self.applied:
            lines.append(f"  [applied]   {c.action} {c.target}: {c.message}")
        for c in self.unresolved:
            lines.append(f"  [attention] {c.action} {c.target}: {c.message}")
        return "\n".join(lines) + "\n"
