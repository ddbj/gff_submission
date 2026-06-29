from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class MssQualifier:
    key: str
    value: str


@dataclass
class MssFeature:
    key: str
    location: str
    qualifiers: list[MssQualifier] = field(default_factory=list)


@dataclass
class MssEntry:
    name: str
    features: list[MssFeature] = field(default_factory=list)


@dataclass
class MssDocument:
    common_rows: list[str] = field(default_factory=list)
    entries: list[MssEntry] = field(default_factory=list)
