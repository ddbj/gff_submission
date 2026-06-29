from __future__ import annotations


class LocusTagAssigner:
    def __init__(self, prefix: str, width: int = 6, start: int = 10, step: int = 10):
        self._prefix = prefix
        self._width = width
        self._next = start
        self._step = step

    @classmethod
    def from_config(cls, cfg) -> "LocusTagAssigner":
        return cls(cfg.locus_tag_prefix, cfg.locus_tag_width,
                   cfg.locus_tag_start, cfg.locus_tag_step)

    def assign(self, feature) -> str:
        existing = feature.locus_tag
        if existing:
            return existing
        tag = f"{self._prefix}_{self._next:0{self._width}d}"
        self._next += self._step
        return tag
