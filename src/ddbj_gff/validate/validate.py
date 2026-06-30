from __future__ import annotations

import dataclasses

from .rules import ALL_RULES
from .severities import resolve_level
from .vocab import load_vocab


def validate(doc, *, severity_overrides: dict[str, str] | None = None) -> list:
    overrides = severity_overrides or {}
    vocab = load_vocab()
    diags: list = []
    for rule in ALL_RULES:
        diags.extend(rule(doc, vocab))

    out: list = []
    for d in diags:
        if d.code in overrides:
            try:
                level = resolve_level(overrides[d.code])
            except ValueError:
                out.append(d)        # unknown level → keep diagnostic, severity unchanged (lenient, §5)
                continue
            if level is None:        # "off"
                continue
            d = dataclasses.replace(d, severity=level)
        out.append(d)

    out.sort(key=lambda d: (d.line_no if d.line_no is not None else -1, d.code))
    return out
