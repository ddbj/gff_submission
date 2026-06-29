# GFF3 パーサ＋オブジェクトモデル（フェーズ1）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** INSDC/SO GFF3 を自作オブジェクトモデルへ読み込み、再び GFF3 へ書き出し、意味的ラウンドトリップが成立する Python ライブラリ `ddbj_gff` を作る。

**Architecture:** GFF の「行」固有情報（座標/strand/phase/score/part）を `Span`、「feature」固有情報（ID/type/source/属性）を `Feature` に分離。同一 ID の複数行を1つの `Feature`（複数 span）へ集約し、Parent で親子グラフを張る。パーサは 2 パス・インメモリ、既定 lenient（診断収集）。Biopython は座標演算・配列 I/O のみに使用。

**Tech Stack:** Python 3.11+ / Biopython（依存はこれのみ）/ pytest / uv / linux-amd64 Docker。

## Global Constraints

- Python `>=3.11`。型注釈は `from __future__ import annotations` 前提で `X | None` を使ってよい。
- 実行依存は **biopython のみ**。bcbio-gff は使用しない。dev 依存は pytest。
- パッケージは `src/` レイアウト、import 名は `ddbj_gff`。
- 開発・テストは `linux/amd64` Ubuntu コンテナ内 `uv run pytest` を基準とする。
- 座標は GFF 流（1-based・inclusive）を真実の源として保持。Biopython へ渡す時のみ 0-based half-open へ変換。
- パーサ既定は lenient（`document.diagnostics` に記録）。`strict=True` で最初の ERROR 診断にて `GffParseError` を送出。
- 文字エンコード: ASCII を想定（非 ASCII は WARNING 記録のうえ受理）。
- 各タスクは「失敗するテスト→失敗確認→最小実装→成功確認→コミット」の順で進める。

---

## File Structure

| ファイル | 責務 |
|---|---|
| `pyproject.toml` / `Dockerfile` / `.devcontainer/devcontainer.json` / `.gitignore` | 環境・パッケージ定義 |
| `src/ddbj_gff/__init__.py` | 公開 API の再エクスポート |
| `src/ddbj_gff/errors.py` | `Severity` / `Diagnostic` / `GffParseError` |
| `src/ddbj_gff/attributes.py` | 列9の `%XX` エンコード/デコード・parse/serialize |
| `src/ddbj_gff/model.py` | `Span` / `Directive` / `Feature` / `GffDocument`（データ構造のみ） |
| `src/ddbj_gff/parser.py` | テキスト → `GffDocument`（ディレクティブ/行/集約/グラフ/FASTA/診断） |
| `src/ddbj_gff/writer.py` | `GffDocument` → GFF3 テキスト |
| `tests/…` | 各モジュールの単体・ラウンドトリップ・統合テスト |
| `tests/fixtures/*.gff3` | 難所を凝縮した小フィクスチャ |

---

## Task 1: プロジェクト雛形と実行環境

**Files:**
- Create: `pyproject.toml`
- Create: `src/ddbj_gff/__init__.py`
- Create: `Dockerfile`
- Create: `.devcontainer/devcontainer.json`
- Create: `.gitignore`
- Create: `tests/__init__.py`
- Create: `tests/test_smoke.py`

**Interfaces:**
- Consumes: なし
- Produces: import 可能なパッケージ `ddbj_gff`、`uv run pytest` が動く環境、`pytest -m slow` マーカー定義。

- [ ] **Step 1: `pyproject.toml` を作成**

```toml
[project]
name = "ddbj-gff"
version = "0.1.0"
description = "INSDC/SO GFF3 parser and object model (DDBJ GFF tools, phase 1)"
requires-python = ">=3.11"
dependencies = ["biopython>=1.83"]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.wheel]
packages = ["src/ddbj_gff"]

[dependency-groups]
dev = ["pytest>=8"]

[tool.pytest.ini_options]
markers = ["slow: integration/perf tests against large example files (deselect with '-m \"not slow\"')"]
testpaths = ["tests"]
addopts = "-m 'not slow'"
```

- [ ] **Step 2: パッケージとテストの雛形を作成**

`src/ddbj_gff/__init__.py`:
```python
"""ddbj_gff: INSDC/SO GFF3 parser and object model (phase 1)."""

__all__ = []
```

`tests/__init__.py`:
```python
```

`tests/test_smoke.py`:
```python
def test_package_imports():
    import ddbj_gff

    assert ddbj_gff is not None
```

- [ ] **Step 3: コンテナ定義を作成**

`Dockerfile`:
```dockerfile
FROM --platform=linux/amd64 ubuntu:24.04

RUN apt-get update && apt-get install -y --no-install-recommends \
        curl ca-certificates git build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:${PATH}"

WORKDIR /workspace
COPY pyproject.toml ./
RUN uv python install 3.11

COPY . .
RUN uv sync

CMD ["bash"]
```

`.devcontainer/devcontainer.json`:
```json
{
  "name": "ddbj-gff",
  "build": { "dockerfile": "../Dockerfile", "context": ".." },
  "runArgs": ["--platform=linux/amd64"],
  "postCreateCommand": "uv sync"
}
```

- [ ] **Step 4: `.gitignore` を作成（巨大データを除外）**

```gitignore
.DS_Store
__pycache__/
*.pyc
.venv/
.pytest_cache/
*.egg-info/
uv.lock

# Large example data (kept out of git; small GFF/GB fixtures are tracked)
examples/**/*.fa
examples/**/*.fasta
examples/**/*.gz
examples/marchantia/MpTak1_v7.1.marpolbase.gff
examples/arabidopsis/AT_chr1.gff3
```

- [ ] **Step 5: 依存を解決しスモークテストが通ることを確認**

Run: `uv sync && uv run pytest tests/test_smoke.py -v`
Expected: `test_package_imports PASSED`

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml src tests Dockerfile .devcontainer .gitignore
git commit -m "chore: scaffold ddbj_gff package, uv env, and amd64 container"
```

---

## Task 2: errors.py（診断と例外）

**Files:**
- Create: `src/ddbj_gff/errors.py`
- Test: `tests/test_errors.py`

**Interfaces:**
- Consumes: なし
- Produces:
  - `class Severity(Enum)` with members `ERROR`, `WARNING`, `INFO` (values are the same strings).
  - `@dataclass(frozen=True) class Diagnostic` fields: `severity: Severity`, `line_no: int | None`, `code: str`, `message: str`.
  - `class GffParseError(Exception)` with attribute `.diagnostic: Diagnostic`, constructed as `GffParseError(diagnostic)`.

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_errors.py`:
```python
import pytest

from ddbj_gff.errors import Diagnostic, GffParseError, Severity


def test_severity_members():
    assert Severity.ERROR.value == "ERROR"
    assert {s.name for s in Severity} == {"ERROR", "WARNING", "INFO"}


def test_diagnostic_is_frozen_and_equal():
    d1 = Diagnostic(Severity.WARNING, 12, "dangling-parent", "Parent x not found")
    d2 = Diagnostic(Severity.WARNING, 12, "dangling-parent", "Parent x not found")
    assert d1 == d2
    with pytest.raises(Exception):
        d1.code = "other"  # frozen dataclass


def test_parse_error_carries_diagnostic():
    d = Diagnostic(Severity.ERROR, 3, "col-count", "expected 9 columns, got 7")
    err = GffParseError(d)
    assert err.diagnostic is d
    assert "col-count" in str(err)
    assert "line 3" in str(err)
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `uv run pytest tests/test_errors.py -v`
Expected: FAIL（`ModuleNotFoundError: No module named 'ddbj_gff.errors'`）

- [ ] **Step 3: 最小実装を書く**

`src/ddbj_gff/errors.py`:
```python
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Severity(Enum):
    ERROR = "ERROR"
    WARNING = "WARNING"
    INFO = "INFO"


@dataclass(frozen=True)
class Diagnostic:
    severity: Severity
    line_no: int | None
    code: str
    message: str


class GffParseError(Exception):
    def __init__(self, diagnostic: Diagnostic):
        self.diagnostic = diagnostic
        line = "?" if diagnostic.line_no is None else diagnostic.line_no
        super().__init__(
            f"{diagnostic.severity.value} (line {line}) "
            f"[{diagnostic.code}] {diagnostic.message}"
        )
```

- [ ] **Step 4: テストが通ることを確認**

Run: `uv run pytest tests/test_errors.py -v`
Expected: 3 passed

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/errors.py tests/test_errors.py
git commit -m "feat(errors): add Severity, Diagnostic, GffParseError"
```

---

## Task 3: attributes.py（列9のエンコード/デコードと parse/serialize）

**Files:**
- Create: `src/ddbj_gff/attributes.py`
- Test: `tests/test_attributes.py`

**Interfaces:**
- Consumes: なし
- Produces:
  - `encode_value(text: str) -> str` — 構造文字 `; = & , %` と制御文字/タブ/改行を `%XX` 化、他は逐語。
  - `decode_value(text: str) -> str` — 任意の `%XX` を復号。
  - `parse_attributes(col9: str) -> dict[str, list[str]]` — 順序保持、多値は `,` 分割、`.`/空は `{}`。
  - `serialize_attributes(attrs: dict[str, list[str]]) -> str` — 逆変換。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_attributes.py`:
```python
from ddbj_gff.attributes import (
    decode_value,
    encode_value,
    parse_attributes,
    serialize_attributes,
)


def test_encode_decode_reserved_roundtrip():
    raw = "a;b=c,d&e%f\tg"
    enc = encode_value(raw)
    assert ";" not in enc and "=" not in enc and "," not in enc
    assert "%3B" in enc and "%3D" in enc and "%2C" in enc and "%25" in enc
    assert decode_value(enc) == raw


def test_parse_simple_and_order_preserved():
    attrs = parse_attributes("ID=gene1;Name=psbA;locus_tag=AKK66")
    assert attrs == {"ID": ["gene1"], "Name": ["psbA"], "locus_tag": ["AKK66"]}
    assert list(attrs.keys()) == ["ID", "Name", "locus_tag"]


def test_parse_multivalue_dbxref():
    attrs = parse_attributes("Dbxref=GenBank:YP_1,GeneID:29")
    assert attrs["Dbxref"] == ["GenBank:YP_1", "GeneID:29"]


def test_parse_decodes_percent_and_keeps_literal_tilde():
    attrs = parse_attributes("Note=LSC%3B~large single-copy region")
    assert attrs["Note"] == ["LSC;~large single-copy region"]


def test_parse_empty_or_dot():
    assert parse_attributes("") == {}
    assert parse_attributes(".") == {}


def test_serialize_roundtrip():
    col9 = "ID=cds1;Note=has%3Bsemicolon;Dbxref=A:1,B:2"
    assert serialize_attributes(parse_attributes(col9)) == col9
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `uv run pytest tests/test_attributes.py -v`
Expected: FAIL（`ModuleNotFoundError: No module named 'ddbj_gff.attributes'`）

- [ ] **Step 3: 最小実装を書く**

`src/ddbj_gff/attributes.py`:
```python
from __future__ import annotations

from urllib.parse import unquote

_RESERVED = {
    "%": "%25",  # must come logically first when encoding
    ";": "%3B",
    "=": "%3D",
    "&": "%26",
    ",": "%2C",
    "\t": "%09",
    "\n": "%0A",
    "\r": "%0D",
}


def encode_value(text: str) -> str:
    out: list[str] = []
    for ch in text:
        if ch in _RESERVED:
            out.append(_RESERVED[ch])
        elif ord(ch) < 0x20:
            out.append("%%%02X" % ord(ch))
        else:
            out.append(ch)
    return "".join(out)


def decode_value(text: str) -> str:
    return unquote(text)


def parse_attributes(col9: str) -> dict[str, list[str]]:
    attrs: dict[str, list[str]] = {}
    if col9 in ("", "."):
        return attrs
    for pair in col9.split(";"):
        if pair == "":
            continue
        if "=" not in pair:
            key = decode_value(pair)
            attrs.setdefault(key, [])
            continue
        raw_key, raw_val = pair.split("=", 1)
        key = decode_value(raw_key)
        values = [decode_value(v) for v in raw_val.split(",")]
        if key in attrs:
            attrs[key].extend(values)
        else:
            attrs[key] = values
    return attrs


def serialize_attributes(attrs: dict[str, list[str]]) -> str:
    parts: list[str] = []
    for key, values in attrs.items():
        ek = encode_value(key)
        ev = ",".join(encode_value(v) for v in values)
        parts.append(f"{ek}={ev}")
    return ";".join(parts)
```

- [ ] **Step 4: テストが通ることを確認**

Run: `uv run pytest tests/test_attributes.py -v`
Expected: 6 passed

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/attributes.py tests/test_attributes.py
git commit -m "feat(attributes): GFF3 column 9 encode/decode and parse/serialize"
```

---

## Task 4: model.py — Span と Directive（値オブジェクト）

**Files:**
- Create: `src/ddbj_gff/model.py`
- Test: `tests/test_model_span.py`

**Interfaces:**
- Consumes: なし
- Produces:
  - `@dataclass class Span` fields: `seqid: str`, `start: int`, `end: int`, `strand: str = "."`, `phase: int | None = None`, `score: float | None = None`, `part: int | None = None`. メソッド `sort_key() -> tuple`.
  - `@dataclass class Directive` fields: `raw: str`, `kind: str`, `value: object = None`.

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_model_span.py`:
```python
from ddbj_gff.model import Directive, Span


def test_span_defaults_and_equality():
    a = Span("chr1", 10, 20)
    b = Span("chr1", 10, 20, strand=".", phase=None, score=None, part=None)
    assert a == b
    assert a.strand == "."


def test_span_sort_key_orders_by_coordinates():
    s1 = Span("chr1", 5, 9, "+", 0)
    s2 = Span("chr1", 10, 20, "+", 1)
    assert s1.sort_key() < s2.sort_key()


def test_span_sort_key_handles_none_phase_and_score():
    # None phase/score must not raise when building the key
    s = Span("chr1", 1, 2)
    key = s.sort_key()
    assert key[0] == "chr1"


def test_directive_holds_kind_and_value():
    d = Directive("##sequence-region chr1 1 100", "sequence-region", ("chr1", 1, 100))
    assert d.kind == "sequence-region"
    assert d.value == ("chr1", 1, 100)
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `uv run pytest tests/test_model_span.py -v`
Expected: FAIL（`ModuleNotFoundError: No module named 'ddbj_gff.model'`）

- [ ] **Step 3: 最小実装を書く**

`src/ddbj_gff/model.py`:
```python
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Span:
    seqid: str
    start: int  # 1-based inclusive
    end: int
    strand: str = "."  # one of: + - . ?
    phase: int | None = None  # 0/1/2 for CDS
    score: float | None = None
    part: int | None = None

    def sort_key(self) -> tuple:
        return (
            self.seqid,
            self.start,
            self.end,
            self.strand,
            -1 if self.phase is None else self.phase,
            float("-inf") if self.score is None else self.score,
            -1 if self.part is None else self.part,
        )


@dataclass
class Directive:
    raw: str
    kind: str
    value: object = None
```

- [ ] **Step 4: テストが通ることを確認**

Run: `uv run pytest tests/test_model_span.py -v`
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/model.py tests/test_model_span.py
git commit -m "feat(model): add Span and Directive value objects"
```

---

## Task 5: model.py — Feature（フィールド・属性アクセサ・ordered_spans・codon_start）

**Files:**
- Modify: `src/ddbj_gff/model.py`
- Test: `tests/test_model_feature.py`

**Interfaces:**
- Consumes: `Span`（Task 4）
- Produces:
  - `@dataclass class Feature` fields: `id: str | None`, `source: str`, `type: str`, `spans: list[Span]`, `attributes: dict[str, list[str]]`, `parent_ids: list[str]`, `children: list[Feature]`, `parents: list[Feature]`（後者2つは `repr=False`）。
  - properties: `name`, `locus_tag`, `gene`, `product`, `protein_id`, `transl_table:int|None`, `note:list[str]`, `dbxref:list[str]`, `number:int|None`, `exon_number:int|None`, `is_circular:bool`, `is_ordered:bool`, `is_trans_spliced:bool`, `codon_start:int|None`。
  - method: `ordered_spans() -> list[Span]`（part 順、無ければ strand 別の座標順）。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_model_feature.py`:
```python
from ddbj_gff.model import Feature, Span


def make_cds(spans):
    return Feature(
        id="cds1",
        source="DDBJ",
        type="CDS",
        spans=spans,
        attributes={"transl_table": ["11"], "Note": ["a", "b"]},
        parent_ids=["gene1"],
    )


def test_feature_property_accessors():
    f = Feature(
        id="g1",
        source="DDBJ",
        type="gene",
        spans=[Span("chr1", 1, 9, "-")],
        attributes={
            "Name": ["rps12"],
            "locus_tag": ["Mp_Cg00010"],
            "gene": ["rps12"],
            "Is_circular": ["true"],
            "is_ordered": ["true"],
            "exception": ["trans-splicing"],
            "number": ["1"],
            "exon_number": ["1"],
            "Dbxref": ["GeneID:1", "X:2"],
        },
        parent_ids=[],
    )
    assert f.name == "rps12"
    assert f.locus_tag == "Mp_Cg00010"
    assert f.gene == "rps12"
    assert f.is_circular is True
    assert f.is_ordered is True
    assert f.is_trans_spliced is True
    assert f.number == 1
    assert f.exon_number == 1
    assert f.dbxref == ["GeneID:1", "X:2"]


def test_trans_spliced_accepts_underscore_or_hyphen():
    f = Feature("x", "s", "CDS", [Span("c", 1, 3, "+", 0)], {"exception": ["trans_splicing"]}, [])
    assert f.is_trans_spliced is True


def test_ordered_spans_plus_strand_by_part():
    spans = [
        Span("c", 100, 110, "+", 0, part=2),
        Span("c", 1, 10, "+", 0, part=1),
    ]
    f = make_cds(spans)
    assert [s.part for s in f.ordered_spans()] == [1, 2]


def test_ordered_spans_minus_strand_descending_when_no_part():
    spans = [
        Span("c", 1, 10, "-", 2),
        Span("c", 100, 110, "-", 0),
    ]
    f = make_cds(spans)
    assert [s.start for s in f.ordered_spans()] == [100, 1]


def test_codon_start_derived_from_first_span_phase():
    spans = [Span("c", 100, 110, "-", 2, part=1), Span("c", 1, 10, "-", 0, part=2)]
    f = make_cds(spans)
    assert f.codon_start == 3  # phase 2 -> codon_start 3


def test_codon_start_none_for_non_cds():
    f = Feature("g", "s", "gene", [Span("c", 1, 9, "+")], {}, [])
    assert f.codon_start is None
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `uv run pytest tests/test_model_feature.py -v`
Expected: FAIL（`ImportError: cannot import name 'Feature'`）

- [ ] **Step 3: 最小実装を書く（`model.py` に追記）**

`src/ddbj_gff/model.py` に以下を追加（`Span`/`Directive` の下、ファイル冒頭の import 群は変更不要）:
```python
@dataclass
class Feature:
    id: str | None
    source: str
    type: str
    spans: list[Span] = field(default_factory=list)
    attributes: dict[str, list[str]] = field(default_factory=dict)
    parent_ids: list[str] = field(default_factory=list)
    children: list["Feature"] = field(default_factory=list, repr=False)
    parents: list["Feature"] = field(default_factory=list, repr=False)

    def _first(self, key: str) -> str | None:
        vals = self.attributes.get(key)
        return vals[0] if vals else None

    def _int(self, key: str) -> int | None:
        v = self._first(key)
        return int(v) if v is not None and v != "" else None

    @property
    def name(self) -> str | None:
        return self._first("Name")

    @property
    def locus_tag(self) -> str | None:
        return self._first("locus_tag")

    @property
    def gene(self) -> str | None:
        return self._first("gene")

    @property
    def product(self) -> str | None:
        return self._first("product")

    @property
    def protein_id(self) -> str | None:
        return self._first("protein_id")

    @property
    def transl_table(self) -> int | None:
        return self._int("transl_table")

    @property
    def number(self) -> int | None:
        return self._int("number")

    @property
    def exon_number(self) -> int | None:
        return self._int("exon_number")

    @property
    def note(self) -> list[str]:
        return self.attributes.get("Note", [])

    @property
    def dbxref(self) -> list[str]:
        return self.attributes.get("Dbxref", [])

    @property
    def is_circular(self) -> bool:
        return self._first("Is_circular") == "true"

    @property
    def is_ordered(self) -> bool:
        return self._first("is_ordered") == "true"

    @property
    def is_trans_spliced(self) -> bool:
        v = self._first("exception")
        return v is not None and v.replace("_", "-") == "trans-splicing"

    def ordered_spans(self) -> list[Span]:
        spans = list(self.spans)
        if not spans:
            return spans
        if any(s.part is not None for s in spans):
            return sorted(spans, key=lambda s: (s.part is None, s.part or 0))
        if spans[0].strand == "-":
            return sorted(spans, key=lambda s: s.start, reverse=True)
        return sorted(spans, key=lambda s: s.start)

    @property
    def codon_start(self) -> int | None:
        if self.type != "CDS" or not self.spans:
            return None
        phase = self.ordered_spans()[0].phase
        return None if phase is None else phase + 1
```

- [ ] **Step 4: テストが通ることを確認**

Run: `uv run pytest tests/test_model_feature.py -v`
Expected: 6 passed

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/model.py tests/test_model_feature.py
git commit -m "feat(model): add Feature with accessors, ordered_spans, codon_start"
```

---

## Task 6: model.py — Feature.to_biopython_location()

**Files:**
- Modify: `src/ddbj_gff/model.py`
- Test: `tests/test_model_location.py`

**Interfaces:**
- Consumes: `Feature`, `Span`（Task 5）、Biopython `FeatureLocation`/`CompoundLocation`。
- Produces: `Feature.to_biopython_location() -> FeatureLocation | CompoundLocation`。1-based→0-based half-open 変換。`ordered_spans()` 順。strand: `+`→1, `-`→-1, それ以外→0。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_model_location.py`:
```python
from Bio.SeqFeature import CompoundLocation, FeatureLocation

from ddbj_gff.model import Feature, Span


def test_single_span_location_is_zero_based_halfopen():
    f = Feature("c", "s", "CDS", [Span("chr1", 10, 20, "+", 0)], {}, [])
    loc = f.to_biopython_location()
    assert isinstance(loc, FeatureLocation)
    assert int(loc.start) == 9
    assert int(loc.end) == 20
    assert loc.strand == 1


def test_multi_span_plus_strand_is_compound_in_part_order():
    f = Feature(
        "c", "s", "CDS",
        [Span("chr1", 100, 110, "+", 0, part=2), Span("chr1", 1, 10, "+", 0, part=1)],
        {}, [],
    )
    loc = f.to_biopython_location()
    assert isinstance(loc, CompoundLocation)
    assert [int(p.start) for p in loc.parts] == [0, 99]


def test_minus_strand_location_strand_is_negative():
    f = Feature("c", "s", "CDS", [Span("chr1", 5, 9, "-", 0)], {}, [])
    loc = f.to_biopython_location()
    assert loc.strand == -1


def test_unknown_strand_maps_to_zero():
    f = Feature("c", "s", "CDS", [Span("chr1", 5, 9, "?", 0)], {}, [])
    loc = f.to_biopython_location()
    assert loc.strand == 0
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `uv run pytest tests/test_model_location.py -v`
Expected: FAIL（`AttributeError: 'Feature' object has no attribute 'to_biopython_location'`）

- [ ] **Step 3: 最小実装を書く**

`src/ddbj_gff/model.py` の import 行を更新（ファイル先頭付近）:
```python
from Bio.SeqFeature import CompoundLocation, FeatureLocation
```

`Feature` クラスに以下メソッドを追加:
```python
    _STRAND_MAP = {"+": 1, "-": -1}

    def to_biopython_location(self):
        parts = []
        for s in self.ordered_spans():
            strand = self._STRAND_MAP.get(s.strand, 0)
            parts.append(FeatureLocation(s.start - 1, s.end, strand=strand))
        if len(parts) == 1:
            return parts[0]
        return CompoundLocation(parts)
```

- [ ] **Step 4: テストが通ることを確認**

Run: `uv run pytest tests/test_model_location.py -v`
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/model.py tests/test_model_location.py
git commit -m "feat(model): Feature.to_biopython_location bridge to Biopython"
```

---

## Task 7: model.py — GffDocument（コンテナとディレクティブアクセサ）

**Files:**
- Modify: `src/ddbj_gff/model.py`
- Test: `tests/test_model_document.py`

**Interfaces:**
- Consumes: `Feature`, `Directive`, `Diagnostic`（errors）。
- Produces:
  - `@dataclass class GffDocument` fields: `directives: list[Directive]`, `features: list[Feature]`, `feature_index: dict[str, Feature]`, `roots: list[Feature]`, `fasta: dict | None = None`, `sequences: dict | None = None`, `diagnostics: list[Diagnostic]`。
  - properties: `gff_version:str|None`, `insdc_gff_version:str|None`, `species:int|None`, `sequence_regions:dict[str,tuple[int,int]]`, `transl_table_map:dict|None`。
  - method: `get(feature_id) -> Feature | None`。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_model_document.py`:
```python
from ddbj_gff.model import Directive, Feature, GffDocument, Span


def build_doc():
    doc = GffDocument()
    doc.directives = [
        Directive("##gff-version 3", "gff-version", "3"),
        Directive("#!insdc-gff-version 1.0.0", "insdc-gff-version", "1.0.0"),
        Directive("##species ...id=4530", "species", 4530),
        Directive("##sequence-region chr1 1 100", "sequence-region", ("chr1", 1, 100)),
        Directive("#!transl_table primary:1", "transl_table", {"primary": 1}),
    ]
    f = Feature("g1", "DDBJ", "gene", [Span("chr1", 1, 9, "+")], {}, [])
    doc.features = [f]
    doc.feature_index = {"g1": f}
    return doc


def test_directive_accessors():
    doc = build_doc()
    assert doc.gff_version == "3"
    assert doc.insdc_gff_version == "1.0.0"
    assert doc.species == 4530
    assert doc.sequence_regions == {"chr1": (1, 100)}
    assert doc.transl_table_map == {"primary": 1}


def test_get_returns_feature_or_none():
    doc = build_doc()
    assert doc.get("g1").id == "g1"
    assert doc.get("missing") is None


def test_defaults_are_independent():
    a = GffDocument()
    b = GffDocument()
    a.features.append("x")
    assert b.features == []
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `uv run pytest tests/test_model_document.py -v`
Expected: FAIL（`ImportError: cannot import name 'GffDocument'`）

- [ ] **Step 3: 最小実装を書く**

`src/ddbj_gff/model.py` の import に追加:
```python
from .errors import Diagnostic
```

ファイル末尾に追加:
```python
@dataclass
class GffDocument:
    directives: list[Directive] = field(default_factory=list)
    features: list[Feature] = field(default_factory=list)
    feature_index: dict[str, Feature] = field(default_factory=dict)
    roots: list[Feature] = field(default_factory=list)
    fasta: dict | None = None
    sequences: dict | None = None
    diagnostics: list[Diagnostic] = field(default_factory=list)

    def _directive(self, kind: str) -> Directive | None:
        for d in self.directives:
            if d.kind == kind:
                return d
        return None

    @property
    def gff_version(self) -> str | None:
        d = self._directive("gff-version")
        return d.value if d else None

    @property
    def insdc_gff_version(self) -> str | None:
        d = self._directive("insdc-gff-version")
        return d.value if d else None

    @property
    def species(self) -> int | None:
        d = self._directive("species")
        return d.value if d else None

    @property
    def sequence_regions(self) -> dict[str, tuple[int, int]]:
        out: dict[str, tuple[int, int]] = {}
        for d in self.directives:
            if d.kind == "sequence-region" and d.value:
                seqid, start, end = d.value
                out[seqid] = (start, end)
        return out

    @property
    def transl_table_map(self) -> dict | None:
        d = self._directive("transl_table")
        return d.value if d else None

    def get(self, feature_id: str) -> Feature | None:
        return self.feature_index.get(feature_id)
```

- [ ] **Step 4: テストが通ることを確認**

Run: `uv run pytest tests/test_model_document.py -v`
Expected: 3 passed

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/model.py tests/test_model_document.py
git commit -m "feat(model): add GffDocument container and directive accessors"
```

---

## Task 8: model.py — semantically_equals()（ラウンドトリップのオラクル）

**Files:**
- Modify: `src/ddbj_gff/model.py`
- Test: `tests/test_model_equality.py`

**Interfaces:**
- Consumes: `GffDocument`, `Feature`, `Span`, `Directive`。
- Produces: `GffDocument.semantically_equals(other: GffDocument) -> bool`。ディレクティブ集合（kind＋正規化 value）と feature 集合（id/type/source/ソート済 spans/属性多重集合/parent_ids/親子の id 集合）の一致で判定。出現順・字面は不問。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_model_equality.py`:
```python
from ddbj_gff.model import Directive, Feature, GffDocument, Span


def doc_with(features, directives=None):
    doc = GffDocument()
    doc.directives = directives or [Directive("##gff-version 3", "gff-version", "3")]
    doc.features = features
    doc.feature_index = {f.id: f for f in features if f.id}
    return doc


def test_equal_ignores_feature_and_attribute_order():
    a = doc_with([
        Feature("g1", "S", "gene", [Span("c", 1, 9, "+")], {"Name": ["x"], "gene": ["x"]}, []),
        Feature("g2", "S", "gene", [Span("c", 20, 29, "+")], {}, []),
    ])
    b = doc_with([
        Feature("g2", "S", "gene", [Span("c", 20, 29, "+")], {}, []),
        Feature("g1", "S", "gene", [Span("c", 1, 9, "+")], {"gene": ["x"], "Name": ["x"]}, []),
    ])
    assert a.semantically_equals(b)


def test_equal_ignores_span_order():
    a = doc_with([Feature("c1", "S", "CDS",
                          [Span("c", 1, 9, "+", 0, part=1), Span("c", 20, 29, "+", 0, part=2)], {}, [])])
    b = doc_with([Feature("c1", "S", "CDS",
                          [Span("c", 20, 29, "+", 0, part=2), Span("c", 1, 9, "+", 0, part=1)], {}, [])])
    assert a.semantically_equals(b)


def test_not_equal_when_span_coords_differ():
    a = doc_with([Feature("g1", "S", "gene", [Span("c", 1, 9, "+")], {}, [])])
    b = doc_with([Feature("g1", "S", "gene", [Span("c", 1, 99, "+")], {}, [])])
    assert not a.semantically_equals(b)


def test_not_equal_when_directive_value_differs():
    a = doc_with([], [Directive("x", "species", 4530)])
    b = doc_with([], [Directive("x", "species", 9606)])
    assert not a.semantically_equals(b)
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `uv run pytest tests/test_model_equality.py -v`
Expected: FAIL（`AttributeError: 'GffDocument' object has no attribute 'semantically_equals'`）

- [ ] **Step 3: 最小実装を書く**

`GffDocument` に追加:
```python
    @staticmethod
    def _directive_key(d: "Directive"):
        v = d.value
        if isinstance(v, dict):
            v = tuple(sorted(v.items()))
        elif isinstance(v, list):
            v = tuple(v)
        return (d.kind, v)

    @staticmethod
    def _feature_key(f: "Feature"):
        return (
            f.id,
            f.type,
            f.source,
            tuple(sorted(s.sort_key() for s in f.spans)),
            tuple(sorted((k, tuple(v)) for k, v in f.attributes.items())),
            tuple(sorted(f.parent_ids)),
            tuple(sorted(c.id for c in f.children if c.id)),
            tuple(sorted(p.id for p in f.parents if p.id)),
        )

    def semantically_equals(self, other: "GffDocument") -> bool:
        if {self._directive_key(d) for d in self.directives} != {
            other._directive_key(d) for d in other.directives
        }:
            return False
        return {self._feature_key(f) for f in self.features} == {
            other._feature_key(f) for f in other.features
        }
```

- [ ] **Step 4: テストが通ることを確認**

Run: `uv run pytest tests/test_model_equality.py -v`
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/model.py tests/test_model_equality.py
git commit -m "feat(model): semantically_equals roundtrip oracle"
```

---

## Task 9: parser.py — parse_directive()

**Files:**
- Create: `src/ddbj_gff/parser.py`
- Test: `tests/test_parser_directive.py`

**Interfaces:**
- Consumes: `Directive`（model）。
- Produces: `parse_directive(line: str) -> Directive`。kind と構造化 value:
  - `gff-version`/`gff-spec-version`/`insdc-gff-version` → `str`
  - `sequence-region` → `(seqid, int, int)`（不正なら `None`）
  - `species` → `int`（taxid。`id=NNN` 抽出。無ければ raw 文字列）
  - `transl_table` → `dict[str,int]`
  - `FASTA` → `None`、`###` → kind `resolution-boundary`
  - その他 → kind `unknown`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_parser_directive.py`:
```python
from ddbj_gff.parser import parse_directive


def test_gff_version():
    d = parse_directive("##gff-version 3")
    assert d.kind == "gff-version" and d.value == "3"


def test_sequence_region():
    d = parse_directive("##sequence-region NC_031333.1 1 134502")
    assert d.kind == "sequence-region"
    assert d.value == ("NC_031333.1", 1, 134502)


def test_species_extracts_taxid():
    d = parse_directive("##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4530")
    assert d.kind == "species" and d.value == 4530


def test_insdc_version_and_spec_version():
    assert parse_directive("#!insdc-gff-version 1.2.3").value == "1.2.3"
    assert parse_directive("#!gff-spec-version 1.21").kind == "gff-spec-version"


def test_transl_table_map():
    d = parse_directive("#!transl_table primary:1,chloroplast:11")
    assert d.kind == "transl_table"
    assert d.value == {"primary": 1, "chloroplast": 11}


def test_fasta_and_boundary_and_unknown():
    assert parse_directive("##FASTA").kind == "FASTA"
    assert parse_directive("###").kind == "resolution-boundary"
    assert parse_directive("##something-else foo").kind == "unknown"
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `uv run pytest tests/test_parser_directive.py -v`
Expected: FAIL（`ModuleNotFoundError: No module named 'ddbj_gff.parser'`）

- [ ] **Step 3: 最小実装を書く**

`src/ddbj_gff/parser.py`:
```python
from __future__ import annotations

import re

from .model import Directive

_TAXID_RE = re.compile(r"id=(\d+)")


def parse_directive(line: str) -> Directive:
    raw = line.rstrip("\r\n")
    if raw.strip() == "###":
        return Directive(raw, "resolution-boundary", None)

    content = raw[2:].strip() if raw.startswith(("##", "#!")) else raw.lstrip("#").strip()
    parts = content.split(None, 1)
    name = parts[0] if parts else ""
    rest = parts[1] if len(parts) > 1 else ""

    if name in ("gff-version", "gff-spec-version", "insdc-gff-version"):
        return Directive(raw, name, rest.strip())
    if name == "sequence-region":
        fields = rest.split()
        if len(fields) >= 3:
            return Directive(raw, "sequence-region", (fields[0], int(fields[1]), int(fields[2])))
        return Directive(raw, "sequence-region", None)
    if name == "species":
        m = _TAXID_RE.search(rest)
        return Directive(raw, "species", int(m.group(1)) if m else rest.strip())
    if name == "transl_table":
        table: dict[str, int] = {}
        for item in re.split(r"[,\s]+", rest.strip()):
            if ":" in item:
                k, v = item.split(":", 1)
                table[k] = int(v)
        return Directive(raw, "transl_table", table)
    if name == "FASTA":
        return Directive(raw, "FASTA", None)
    return Directive(raw, "unknown", rest if rest else None)
```

- [ ] **Step 4: テストが通ることを確認**

Run: `uv run pytest tests/test_parser_directive.py -v`
Expected: 6 passed

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/parser.py tests/test_parser_directive.py
git commit -m "feat(parser): structured directive parsing"
```

---

## Task 10: parser.py — parse_feature_line()

**Files:**
- Modify: `src/ddbj_gff/parser.py`
- Test: `tests/test_parser_line.py`

**Interfaces:**
- Consumes: `Span`（model）, `parse_attributes`（attributes）, `Diagnostic`/`Severity`（errors）。
- Produces:
  - `@dataclass class ParsedRow` fields: `id: str | None`, `source: str`, `type: str`, `span: Span`, `attributes: dict[str, list[str]]`, `parent_ids: list[str]`, `line_no: int`。
  - `parse_feature_line(line: str, line_no: int, diagnostics: list[Diagnostic]) -> ParsedRow | None`。`part` は `span.part` に取り込み、`attributes` からは除去。列数不正/座標非整数で `None`＋ERROR、start>end と非 ASCII は WARNING。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_parser_line.py`:
```python
from ddbj_gff.errors import Severity
from ddbj_gff.parser import parse_feature_line


def test_parse_basic_cds_line():
    line = "chr1\tDDBJ\tCDS\t10\t20\t.\t+\t0\tID=cds1;Parent=g1;part=1;product=PsbA"
    diags = []
    row = parse_feature_line(line, 5, diags)
    assert row.id == "cds1"
    assert row.type == "CDS"
    assert row.parent_ids == ["g1"]
    assert row.span.start == 10 and row.span.end == 20
    assert row.span.strand == "+" and row.span.phase == 0
    assert row.span.part == 1
    assert "part" not in row.attributes  # part lifted onto the span
    assert row.attributes["product"] == ["PsbA"]
    assert diags == []


def test_wrong_column_count_records_error_and_returns_none():
    diags = []
    row = parse_feature_line("a\tb\tc", 3, diags)
    assert row is None
    assert diags[0].severity == Severity.ERROR
    assert diags[0].code == "col-count"


def test_non_integer_coord_records_error():
    diags = []
    row = parse_feature_line("c\ts\tgene\tx\t20\t.\t+\t.\tID=g", 2, diags)
    assert row is None
    assert diags[0].code == "coord"


def test_start_gt_end_warns_but_parses():
    diags = []
    row = parse_feature_line("c\ts\tCDS\t150\t10\t.\t+\t0\tID=g", 1, diags)
    assert row is not None
    assert any(d.code == "start-gt-end" and d.severity == Severity.WARNING for d in diags)


def test_score_and_phase_dot_become_none():
    diags = []
    row = parse_feature_line("c\ts\texon\t1\t9\t.\t-\t.\tID=e", 1, diags)
    assert row.span.score is None and row.span.phase is None
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `uv run pytest tests/test_parser_line.py -v`
Expected: FAIL（`ImportError: cannot import name 'parse_feature_line'`）

- [ ] **Step 3: 最小実装を書く**

`src/ddbj_gff/parser.py` の import を更新:
```python
from dataclasses import dataclass

from .attributes import parse_attributes
from .errors import Diagnostic, Severity
from .model import Directive, Span
```

ファイル末尾に追加:
```python
@dataclass
class ParsedRow:
    id: str | None
    source: str
    type: str
    span: Span
    attributes: dict[str, list[str]]
    parent_ids: list[str]
    line_no: int


def parse_feature_line(
    line: str, line_no: int, diagnostics: list[Diagnostic]
) -> ParsedRow | None:
    cols = line.rstrip("\r\n").split("\t")
    if len(cols) != 9:
        diagnostics.append(
            Diagnostic(Severity.ERROR, line_no, "col-count", f"expected 9 columns, got {len(cols)}")
        )
        return None
    seqid, source, ftype, start_s, end_s, score_s, strand, phase_s, attr_s = cols
    try:
        start = int(start_s)
        end = int(end_s)
    except ValueError:
        diagnostics.append(
            Diagnostic(Severity.ERROR, line_no, "coord", f"non-integer start/end: {start_s!r},{end_s!r}")
        )
        return None
    score = None if score_s == "." else float(score_s)
    phase = None if phase_s == "." else int(phase_s)

    if start > end:
        diagnostics.append(
            Diagnostic(Severity.WARNING, line_no, "start-gt-end",
                       f"start>end ({start}>{end}); possible origin-spanning feature")
        )
    if not attr_s.isascii():
        diagnostics.append(
            Diagnostic(Severity.WARNING, line_no, "non-ascii", "non-ASCII characters in attributes")
        )

    attrs = parse_attributes(attr_s)
    part = None
    if attrs.get("part"):
        part = int(attrs["part"][0])
        attrs.pop("part", None)

    span = Span(seqid, start, end, strand, phase, score, part)
    fid = attrs.get("ID", [None])[0]
    parent_ids = list(attrs.get("Parent", []))
    return ParsedRow(fid, source, ftype, span, attrs, parent_ids, line_no)
```

- [ ] **Step 4: テストが通ることを確認**

Run: `uv run pytest tests/test_parser_line.py -v`
Expected: 5 passed

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/parser.py tests/test_parser_line.py
git commit -m "feat(parser): parse a single feature line into ParsedRow"
```

---

## Task 11: parser.py — parse()（行分類・集約・グラフ・FASTA・診断・strict）

**Files:**
- Modify: `src/ddbj_gff/parser.py`
- Modify: `src/ddbj_gff/__init__.py`
- Test: `tests/test_parser_parse.py`

**Interfaces:**
- Consumes: `parse_directive`, `parse_feature_line`, `ParsedRow`（同モジュール）, `GffDocument`, `Feature`（model）, `GffParseError`（errors）, Biopython `SeqIO`。
- Produces: `parse(text: str, *, strict: bool = False) -> GffDocument`。同 ID 集約、parent グラフ解決、`##FASTA` ペプチド読込、診断収集、`strict=True` で最初の ERROR にて `GffParseError`。`__init__` から `parse`, `GffDocument`, `Feature`, `Span`, `Directive` を再エクスポート。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_parser_parse.py`:
```python
import pytest

from ddbj_gff.errors import GffParseError, Severity
from ddbj_gff.parser import parse

CANONICAL = """\
##gff-version 3
##sequence-region chr1 1 100
chr1\tS\tgene\t1\t99\t.\t+\t.\tID=g1;Name=x
chr1\tS\tmRNA\t1\t99\t.\t+\t.\tID=m1;Parent=g1
chr1\tS\tCDS\t1\t10\t.\t+\t0\tID=c1;Parent=m1;part=1
chr1\tS\tCDS\t50\t99\t.\t+\t2\tID=c1;Parent=m1;part=2
"""


def test_directives_and_features_parsed():
    doc = parse(CANONICAL)
    assert doc.gff_version == "3"
    assert doc.sequence_regions == {"chr1": (1, 100)}
    assert {f.id for f in doc.features} == {"g1", "m1", "c1"}


def test_same_id_rows_aggregate_into_one_feature_with_spans():
    doc = parse(CANONICAL)
    cds = doc.get("c1")
    assert len(cds.spans) == 2
    assert [s.part for s in cds.ordered_spans()] == [1, 2]


def test_parent_child_graph_resolved():
    doc = parse(CANONICAL)
    g1 = doc.get("g1")
    assert [c.id for c in g1.children] == ["m1"]
    assert doc.get("m1").parents[0].id == "g1"
    assert [r.id for r in doc.roots] == ["g1"]


def test_forward_parent_reference_resolves():
    text = (
        "chr1\tS\tmRNA\t1\t99\t.\t+\t.\tID=m1;Parent=g1\n"
        "chr1\tS\tgene\t1\t99\t.\t+\t.\tID=g1\n"
    )
    doc = parse(text)
    assert doc.get("g1").children[0].id == "m1"


def test_dangling_parent_records_warning():
    doc = parse("chr1\tS\tmRNA\t1\t9\t.\t+\t.\tID=m1;Parent=ghost\n")
    assert any(d.code == "dangling-parent" and d.severity == Severity.WARNING
               for d in doc.diagnostics)


def test_no_id_rows_become_standalone_features():
    doc = parse("chr1\tS\tregion\t1\t9\t.\t+\t.\tNote=x\n")
    assert len(doc.features) == 1
    assert doc.features[0].id is None


def test_peptide_fasta_loaded():
    text = CANONICAL + "##FASTA\n>c1\nMAAA\n"
    doc = parse(text)
    assert str(doc.fasta["c1"]) == "MAAA"


def test_strict_raises_on_error_diagnostic():
    with pytest.raises(GffParseError):
        parse("too\tfew\tcols\n", strict=True)
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `uv run pytest tests/test_parser_parse.py -v`
Expected: FAIL（`ImportError: cannot import name 'parse'`）

- [ ] **Step 3: 最小実装を書く**

`src/ddbj_gff/parser.py` の import に追加:
```python
from io import StringIO

from Bio import SeqIO

from .errors import GffParseError
from .model import Feature, GffDocument
```

ファイル末尾に追加:
```python
def _add_row(doc: GffDocument, row: ParsedRow) -> None:
    if row.id is not None and row.id in doc.feature_index:
        feat = doc.feature_index[row.id]
        if feat.type != row.type:
            doc.diagnostics.append(
                Diagnostic(
                    Severity.ERROR, row.line_no, "id-type-mismatch",
                    f"ID {row.id!r} reused with different type ({feat.type} vs {row.type})",
                )
            )
            doc.features.append(
                Feature(row.id, row.source, row.type, [row.span], dict(row.attributes), list(row.parent_ids))
            )
            return
        feat.spans.append(row.span)
        return
    feat = Feature(row.id, row.source, row.type, [row.span], dict(row.attributes), list(row.parent_ids))
    doc.features.append(feat)
    if row.id is not None:
        doc.feature_index[row.id] = feat


def _resolve_graph(doc: GffDocument) -> None:
    for feat in doc.features:
        for pid in feat.parent_ids:
            parent = doc.feature_index.get(pid)
            if parent is None:
                doc.diagnostics.append(
                    Diagnostic(Severity.WARNING, None, "dangling-parent",
                               f"Parent {pid!r} not found for feature {feat.id!r}")
                )
                continue
            parent.children.append(feat)
            feat.parents.append(parent)
    doc.roots = [f for f in doc.features if not f.parent_ids]


def _parse_fasta(lines: list[str]) -> dict:
    handle = StringIO("\n".join(lines))
    return {rec.id: rec.seq for rec in SeqIO.parse(handle, "fasta")}


def parse(text: str, *, strict: bool = False) -> GffDocument:
    doc = GffDocument()
    in_fasta = False
    fasta_lines: list[str] = []
    dropped_comments = 0

    for line_no, line in enumerate(text.splitlines(), start=1):
        if in_fasta:
            fasta_lines.append(line)
            continue
        if line == "":
            continue
        if line.startswith("#"):
            if line.startswith(("##", "#!")) or line.strip() == "###":
                directive = parse_directive(line)
                doc.directives.append(directive)
                if directive.kind == "FASTA":
                    in_fasta = True
            else:
                dropped_comments += 1
            continue
        row = parse_feature_line(line, line_no, doc.diagnostics)
        if row is not None:
            _add_row(doc, row)

    if fasta_lines:
        doc.fasta = _parse_fasta(fasta_lines)
    if dropped_comments:
        doc.diagnostics.append(
            Diagnostic(Severity.INFO, None, "dropped-comments", f"{dropped_comments} bare comment line(s) ignored")
        )

    _resolve_graph(doc)

    if strict:
        for d in doc.diagnostics:
            if d.severity == Severity.ERROR:
                raise GffParseError(d)
    return doc
```

`src/ddbj_gff/__init__.py` を更新:
```python
"""ddbj_gff: INSDC/SO GFF3 parser and object model (phase 1)."""

from .model import Directive, Feature, GffDocument, Span
from .parser import parse

__all__ = ["parse", "GffDocument", "Feature", "Span", "Directive"]
```

- [ ] **Step 4: テストが通ることを確認**

Run: `uv run pytest tests/test_parser_parse.py -v`
Expected: 8 passed

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/parser.py src/ddbj_gff/__init__.py tests/test_parser_parse.py
git commit -m "feat(parser): full parse with aggregation, graph, FASTA, diagnostics"
```

---

## Task 12: writer.py — write()

**Files:**
- Create: `src/ddbj_gff/writer.py`
- Modify: `src/ddbj_gff/__init__.py`
- Test: `tests/test_writer.py`

**Interfaces:**
- Consumes: `GffDocument`, `Feature`, `Span`（model）, `serialize_attributes`（attributes）。
- Produces: `write(doc: GffDocument, *, canonical_sort: bool = False) -> str`。ヘッダディレクティブ（`FASTA` 以外）→ feature（`document.features` 順、`canonical_sort=True` で seqid→start ソート）を 1 span=1 行で出力（`part` は span から注入）→ 末尾に `##FASTA`＋ペプチド。`__init__` から `write` を再エクスポート。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_writer.py`:
```python
from ddbj_gff.model import Directive, Feature, GffDocument, Span
from ddbj_gff.writer import write


def test_write_one_line_per_span_and_part_injected():
    f = Feature(
        "c1", "S", "CDS",
        [Span("chr1", 1, 10, "+", 0, part=1), Span("chr1", 50, 99, "+", 2, part=2)],
        {"ID": ["c1"], "Parent": ["m1"]},
        ["m1"],
    )
    doc = GffDocument(directives=[Directive("##gff-version 3", "gff-version", "3")], features=[f])
    text = write(doc)
    rows = [l for l in text.splitlines() if not l.startswith("#")]
    assert len(rows) == 2
    assert rows[0].split("\t")[3:5] == ["1", "10"]
    assert "part=1" in rows[0]
    assert "part=2" in rows[1]
    assert rows[1].split("\t")[7] == "2"  # phase of span 2


def test_header_directives_emitted_fasta_moved_to_end():
    f = Feature("g1", "S", "gene", [Span("c", 1, 9, "+")], {"ID": ["g1"]}, [])
    doc = GffDocument(
        directives=[
            Directive("##gff-version 3", "gff-version", "3"),
            Directive("##FASTA", "FASTA", None),
        ],
        features=[f],
        fasta={"g1": "MAAA"},
    )
    text = write(doc)
    lines = text.splitlines()
    assert lines[0] == "##gff-version 3"
    assert "##FASTA" in lines
    assert lines.index("##gff-version 3") < lines.index("##FASTA")
    assert lines.index("##FASTA") > lines.index([l for l in lines if l.startswith("c\t")][0])
    assert ">g1" in lines


def test_score_and_phase_dot_when_none():
    f = Feature("e", "S", "exon", [Span("c", 1, 9, "-")], {"ID": ["e"]}, [])
    doc = GffDocument(features=[f])
    row = [l for l in write(doc).splitlines() if l.startswith("c\t")][0]
    cols = row.split("\t")
    assert cols[5] == "." and cols[7] == "."
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `uv run pytest tests/test_writer.py -v`
Expected: FAIL（`ModuleNotFoundError: No module named 'ddbj_gff.writer'`）

- [ ] **Step 3: 最小実装を書く**

`src/ddbj_gff/writer.py`:
```python
from __future__ import annotations

from .attributes import serialize_attributes
from .model import Feature, GffDocument, Span


def _fmt_score(score: float | None) -> str:
    return "." if score is None else ("%g" % score)


def _format_row(feat: Feature, span: Span) -> str:
    attrs = dict(feat.attributes)
    if span.part is not None:
        attrs["part"] = [str(span.part)]
    col9 = serialize_attributes(attrs)
    phase = "." if span.phase is None else str(span.phase)
    return "\t".join(
        [
            span.seqid,
            feat.source,
            feat.type,
            str(span.start),
            str(span.end),
            _fmt_score(span.score),
            span.strand,
            phase,
            col9,
        ]
    )


def _format_fasta(fasta: dict) -> str:
    out: list[str] = []
    for fid, seq in fasta.items():
        out.append(f">{fid}")
        s = str(seq)
        for i in range(0, len(s), 60):
            out.append(s[i : i + 60])
    return "\n".join(out) + "\n"


def _canonical_sort(features: list[Feature]) -> list[Feature]:
    def key(f: Feature):
        s = f.ordered_spans()
        first = s[0] if s else Span("", 0, 0)
        return (first.seqid, first.start, first.end)

    return sorted(features, key=key)


def write(doc: GffDocument, *, canonical_sort: bool = False) -> str:
    lines: list[str] = []
    for d in doc.directives:
        if d.kind == "FASTA":
            continue
        lines.append(d.raw)

    features = _canonical_sort(doc.features) if canonical_sort else doc.features
    for feat in features:
        for span in feat.ordered_spans():
            lines.append(_format_row(feat, span))

    text = "\n".join(lines)
    if text:
        text += "\n"
    if doc.fasta:
        text += "##FASTA\n" + _format_fasta(doc.fasta)
    return text
```

`src/ddbj_gff/__init__.py` を更新（`write` を追加）:
```python
"""ddbj_gff: INSDC/SO GFF3 parser and object model (phase 1)."""

from .model import Directive, Feature, GffDocument, Span
from .parser import parse
from .writer import write

__all__ = ["parse", "write", "GffDocument", "Feature", "Span", "Directive"]
```

- [ ] **Step 4: テストが通ることを確認**

Run: `uv run pytest tests/test_writer.py -v`
Expected: 3 passed

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/writer.py src/ddbj_gff/__init__.py tests/test_writer.py
git commit -m "feat(writer): serialize GffDocument back to GFF3"
```

---

## Task 13: ラウンドトリップテストとキュレートフィクスチャ

**Files:**
- Create: `tests/fixtures/canonical_gene.gff3`
- Create: `tests/fixtures/discontinuous_cds.gff3`
- Create: `tests/fixtures/trans_splicing.gff3`
- Create: `tests/fixtures/transl_except.gff3`
- Create: `tests/fixtures/circular.gff3`
- Create: `tests/fixtures/attributes_escaping.gff3`
- Create: `tests/fixtures/peptide_fasta.gff3`
- Test: `tests/test_roundtrip.py`

**Interfaces:**
- Consumes: `parse`, `write`（公開 API）。
- Produces: 各フィクスチャで `parse(write(parse(t))).semantically_equals(parse(t))` を保証。

- [ ] **Step 1: フィクスチャを作成**

`tests/fixtures/canonical_gene.gff3`:
```text
##gff-version 3
##sequence-region chr1 1 6000
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9823
chr1	S	gene	456	5212	.	+	.	ID=gene_FABP5;gene=FABP5
chr1	S	mRNA	456	616	.	+	.	ID=mRNA_FABP5;Parent=gene_FABP5
chr1	S	mRNA	3103	3275	.	+	.	ID=mRNA_FABP5;Parent=gene_FABP5
chr1	S	exon	456	616	.	+	.	ID=exon_FABP5;Parent=mRNA_FABP5;Note=number%3D1
chr1	S	CDS	538	616	.	+	0	ID=CDS_FABP5;Parent=mRNA_FABP5;protein_id=ACA05023.1
chr1	S	CDS	3103	3275	.	+	2	ID=CDS_FABP5;Parent=mRNA_FABP5;protein_id=ACA05023.1
```

`tests/fixtures/discontinuous_cds.gff3`:
```text
##gff-version 3
##sequence-region NC_031333.1 1 134502
NC_031333.1	RefSeq	CDS	43658	43789	.	-	0	ID=cds-ycf3;Parent=gene-ycf3;product=Ycf3
NC_031333.1	RefSeq	CDS	42694	42921	.	-	0	ID=cds-ycf3;Parent=gene-ycf3;product=Ycf3
NC_031333.1	RefSeq	CDS	41811	41969	.	-	0	ID=cds-ycf3;Parent=gene-ycf3;product=Ycf3
```

`tests/fixtures/trans_splicing.gff3`:
```text
##gff-version 3
##sequence-region AP025455.1 1 120306
AP025455.1	DDBJ	gene	65903	66802	.	-	.	ID=gene-Mp_Cg00010;Name=rps12;exception=trans-splicing;is_ordered=true;locus_tag=Mp_Cg00010;part=1
AP025455.1	DDBJ	gene	1	854	.	?	.	ID=gene-Mp_Cg00010;Name=rps12;exception=trans-splicing;is_ordered=true;locus_tag=Mp_Cg00010;part=2
AP025455.1	DDBJ	intron	65903	66688	.	-	.	ID=id-Mp_Cg00010;Parent=gene-Mp_Cg00010;exception=trans-splicing;exon_number=1;number=1;part=1
AP025455.1	DDBJ	intron	1	92	.	?	.	ID=id-Mp_Cg00010;Parent=gene-Mp_Cg00010;exception=trans-splicing;exon_number=1;number=1;part=2
```

`tests/fixtures/transl_except.gff3`:
```text
##gff-version 3
##sequence-region NC_000913.3 1 4641652
NC_000913.3	RefSeq	CDS	1547401	1550448	.	+	0	ID=cds-fdnG;Parent=gene-b1474;product=formate dehydrogenase N subunit alpha;transl_except=(pos:1547986..1547988%2Caa:Sec);transl_table=11
NC_000913.3	RefSeq	CDS	492092	493375	.	+	0	ID=cds-dnaX;Parent=gene-b0470;exception=ribosomal slippage;transl_table=11
NC_000913.3	RefSeq	CDS	493375	493386	.	+	0	ID=cds-dnaX;Parent=gene-b0470;exception=ribosomal slippage;transl_table=11
```

`tests/fixtures/circular.gff3`:
```text
##gff-version 3
##sequence-region NC_000913.3 1 5000
NC_000913.3	RefSeq	region	1	5000	.	+	.	ID=NC_000913.3:1..5000;Is_circular=true
NC_000913.3	RefSeq	CDS	4900	5200	.	+	0	ID=cds-origin;Parent=gene-origin;product=origin spanning
```

`tests/fixtures/attributes_escaping.gff3`:
```text
##gff-version 3
##sequence-region chr1 1 1000
chr1	S	mRNA	1	100	.	+	.	ID=m1;Note=has%3Bsemicolon and%2Ccomma;Dbxref=GenBank:NM_1,GeneID:29;product=protein %28X%29
```

`tests/fixtures/peptide_fasta.gff3`:
```text
##gff-version 3
##sequence-region chr1 1 1000
chr1	S	CDS	1	30	.	+	0	ID=c1;Parent=m1;protein_id=ACA05023.1
##FASTA
>c1
MASIQQLVGRWRLV
```

- [ ] **Step 2: 失敗するテストを書く**

`tests/test_roundtrip.py`:
```python
from pathlib import Path

import pytest

from ddbj_gff import parse, write

FIXTURES = Path(__file__).parent / "fixtures"
FILES = [
    "canonical_gene.gff3",
    "discontinuous_cds.gff3",
    "trans_splicing.gff3",
    "transl_except.gff3",
    "circular.gff3",
    "attributes_escaping.gff3",
    "peptide_fasta.gff3",
]


@pytest.mark.parametrize("name", FILES)
def test_semantic_roundtrip(name):
    text = (FIXTURES / name).read_text()
    once = parse(text)
    twice = parse(write(once))
    assert twice.semantically_equals(once)


def test_trans_splicing_structure_preserved():
    text = (FIXTURES / "trans_splicing.gff3").read_text()
    doc = parse(text)
    gene = doc.get("gene-Mp_Cg00010")
    assert len(gene.spans) == 2
    assert gene.is_trans_spliced and gene.is_ordered
    intron = doc.get("id-Mp_Cg00010")
    assert intron.number == 1  # number is the qualifier, not the part order
    assert [s.part for s in intron.ordered_spans()] == [1, 2]


def test_discontinuous_cds_has_three_spans():
    text = (FIXTURES / "discontinuous_cds.gff3").read_text()
    assert len(parse(text).get("cds-ycf3").spans) == 3
```

- [ ] **Step 3: テストを実行して結果を確認**

Run: `uv run pytest tests/test_roundtrip.py -v`
Expected: 全 9 件（parametrize 7＋構造2）が PASS。FAIL する場合は `semantically_equals` が落とす差分（多くは `part`/属性の扱い）を `write`/`parse` 側で修正してから再実行。

- [ ] **Step 4: Commit**

```bash
git add tests/fixtures tests/test_roundtrip.py
git commit -m "test: curated fixtures and semantic roundtrip coverage"
```

---

## Task 14: 実ファイル統合テスト（slow・存在時のみ）

**Files:**
- Create: `tests/test_integration.py`

**Interfaces:**
- Consumes: `parse`, `write`（公開 API）, `examples/` の実ファイル。
- Produces: 実データでの往復＋不変条件。`@pytest.mark.slow`。ファイルが無ければ skip。`-m slow` で実行。

- [ ] **Step 1: 統合テストを書く**

`tests/test_integration.py`:
```python
import gzip
from pathlib import Path

import pytest

from ddbj_gff import parse, write

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"

pytestmark = pytest.mark.slow


def _read_text(path: Path) -> str:
    if path.suffix == ".gz":
        return gzip.decompress(path.read_bytes()).decode("ascii", errors="replace")
    return path.read_text(errors="replace")


def _require(path: Path) -> str:
    if not path.exists():
        pytest.skip(f"example file not present: {path}")
    return _read_text(path)


def test_rice_cp_roundtrip_and_trans_splicing():
    text = _require(EXAMPLES / "rice_cp" / "rice_cp.gff3")
    doc = parse(text)
    # rps12 trans-spliced CDS has 3 segments sharing one ID
    cds = doc.get("cds-YP_009305283.1")
    assert cds is not None and len(cds.spans) == 3
    assert parse(write(doc)).semantically_equals(doc)


def test_chloroplast_ddbj_roundtrip():
    text = _require(EXAMPLES / "marchantia" / "chloroplast.gff3")
    doc = parse(text)
    assert parse(write(doc)).semantically_equals(doc)


def test_ecoli_transl_except_present_and_roundtrip():
    text = _require(EXAMPLES / "ecoli" / "GCF_000005845.2_ASM584v2_genomic.gff.gz")
    doc = parse(text)
    has_transl_except = any("transl_except" in f.attributes for f in doc.features)
    assert has_transl_except
    assert parse(write(doc)).semantically_equals(doc)


def test_arabidopsis_parses_within_budget():
    import time

    text = _require(EXAMPLES / "arabidopsis" / "AT_chr1.gff3")
    t0 = time.perf_counter()
    doc = parse(text)
    elapsed = time.perf_counter() - t0
    assert len(doc.features) > 1000
    assert elapsed < 120  # loose ceiling: 121MB must parse in-memory in reasonable time
```

- [ ] **Step 2: 統合テストを実行**

Run: `uv run pytest tests/test_integration.py -m slow -v`
Expected: 存在するファイルは PASS、未配置（gz/121MB は `.gitignore` 対象）は SKIP。少なくとも `rice_cp.gff3` と `chloroplast.gff3`（git 追跡対象）は PASS。

- [ ] **Step 3: 全テストの確認（slow 除外の通常実行）**

Run: `uv run pytest -v`
Expected: Task 2–13 の全テストが PASS（slow は addopts によりデフォルト除外）。

- [ ] **Step 4: Commit**

```bash
git add tests/test_integration.py
git commit -m "test: slow integration tests against real example files"
```

---

## Self-Review

**1. Spec coverage**（spec §章 → 実装タスク）:
- §3 構成・コンテナ・巨大ファイル → Task 1
- §4 Span/Directive/Feature/GffDocument/属性3層/Biopython 境界 → Task 4–7（属性アクセサ Task 5、`number`≠`part` は Task 5/10/13、`to_biopython_location` Task 6）
- §4.1 集約ルール → Task 11（`_add_row`）、§4.2 親子グラフ → Task 11（`_resolve_graph`）
- §5.1 2 パス・行分類・FASTA・`###` → Task 9/11、§5.2 ライター（初出順＋canonical オプション）→ Task 12
- §5.3 意味的ラウンドトリップ・`semantically_equals` → Task 8/13
- §5.4 lenient＋診断・strict → Task 10/11
- §6 テスト戦略（attributes/model/parser/writer/roundtrip/integration）→ Task 2–14
- §7 スコープ境界（解釈・変換はしない）→ 計画は表現のみ実装（正規化・MSS・flat→GFF・翻訳は対象外）
- §8 検証データ → Task 13/14

ギャップ: なし（正規化・MSS 出力等は仕様で Phase2/3 と明記、本計画の対象外）。

**2. Placeholder scan**: 各ステップに実コードを記載。"TBD"/"後で"/"適宜" は不使用。

**3. Type consistency**: `parse(text, *, strict)`, `write(doc, *, canonical_sort)`, `Span(seqid,start,end,strand,phase,score,part)`, `Feature(id,source,type,spans,attributes,parent_ids,children,parents)`, `ordered_spans()`, `to_biopython_location()`, `semantically_equals()`, `ParsedRow`, `parse_directive`, `parse_feature_line`, `parse_attributes`/`serialize_attributes`/`encode_value`/`decode_value`, `Severity`/`Diagnostic`/`GffParseError` — タスク間で名称・シグネチャ一致を確認済み。
