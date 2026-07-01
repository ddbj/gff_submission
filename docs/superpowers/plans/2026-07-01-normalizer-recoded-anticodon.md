# 特殊ケース① recoded_codon / anticodon ＋ transl_except 翻訳 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 3-B に `pass_transl_except`（CDS `transl_except` → `recoded_codon`/`stop_codon` 子feature）と `pass_anticodon`（tRNA `anticodon` → `anticodon` 子feature）を追加し、Phase 2 CDS 翻訳を transl_except 対応のテスト済み関数に置換する。

**Architecture:** `normalize/passes.py` に stdlib のみの2パス＋pos-spec パーサ、共有 `aa_names` モジュール。`mss/translate.py` に翻訳関数を verbatim vendoring し、`mss/convert.py` の `build_cds_feature` が transl_except 属性＋recoded_codon 子から例外を集めて適用。

**Tech Stack:** Python 3.11+ / 既存 `ddbj_gff` / Biopython（Phase2 翻訳のみ） / pytest / dev コンテナ。

## Global Constraints

- 変更は `src/ddbj_gff/normalize/passes.py`・`normalize.py`、新規 `src/ddbj_gff/aa_names.py`・`src/ddbj_gff/mss/translate.py`、`src/ddbj_gff/mss/convert.py`、`src/ddbj_gff/validate/rules.py`（`rule_feature_type` に特殊型 accept-set を小追加・Task 5 のみ）、及びテスト。`vocab.py`・`pass_directives`・`pass_so_terms`・CLI・バンドルデータは不変。
- normalize パス（passes.py）は **stdlib のみ**（pos は正規表現、子feature 生成は stdlib）。Biopython は `mss/translate.py`・`convert.py` に閉じる。
- 子feature は `Feature(id, source, type, spans, attributes, parent_ids)`。writer は col9 を `attributes` から再構築するので、子の `attributes` に **`ID` と `Parent` を必ず入れる**（かつ `parent_ids` にも親IDを入れる）。
- recoded_codon 子: col3=`recoded_codon`、`codon_redefined=<full name>`、Span.phase=0。stop_codon 子（aa が終止）: col3=`stop_codon`。anticodon 子: col3=`anticodon`、`amino_acid=<full name>`・`sequence=<seq>`、Span.phase=None。
- pos は単一 range 前提。パース不能/`join`含む/親境界外 → 変換せず `needs-manual` report、当該 spec は属性に残す。
- `transl_except` 属性値は Phase1 が URL-decode 済み（`['(pos:139..141,aa:Sec)']` の形）。
- `ALL_PASSES` 末尾に2パスを追加。`normalize._APPLIED` に `"add-child-feature"` を追加。
- vendored `translate.py` は `nigyta/translate_with_exception` @ `d3c382242f1372afb2b49b47a245ba8dcf548cf4`（NIG 自前コード・再利用承認済み）を verbatim＋provenance ヘッダ。
- **テストは dev コンテナ内**: `docker exec ddbj-gff-dev uv run pytest …`。`git`/`curl` は host 可。TDD。

## 既存インターフェース（参照）

- `model.Feature(id, source, type, spans=[], attributes={}, parent_ids=[], children=[], parents=[])`; `Feature.type` 可変; `.attributes: dict[str,list[str]]`; `.children`; `.to_biopython_location()`（`ordered_spans` を Bio Location に）。`model.Span(seqid, start, end, strand=".", phase=None, score=None, part=None)`。`GffDocument.features`/`.feature_index`。
- `normalize.passes`: `NormalizeContext(vocab, seq_lengths, config)`、`from ..model import Directive`、`from .report import Change`。`report.Change(action, target, message)`。
- `normalize.normalize`: `ALL_PASSES = [pass_directives, pass_so_terms]`、`_APPLIED = {"add-directive","rename-type","add-qualifier"}`。
- `mss/convert.py build_cds_feature(mrna, gene, locus_tag, genome_seq, cfg, diagnostics)`: `spans=collect_spans(mrna,"CDS")`; `codon_start` from phase; `table_id` from CDS `transl_table` else cfg; 現状 `protein = str(Seq(coding_full).translate(table=table_id))`（line 154）＋ `translation-internal-stop`(line 157)/`translation-no-start`(line 160) 診断。`from Bio.Seq import Seq` 既存。
- vendored 公開 API: `translate_cds_with_transl_except(feature: Bio.SeqFeature.SeqFeature, parent_seq: Seq|SeqRecord, stop_symbol="*") -> Seq`。`feature.qualifiers` の `transl_except`(list of `(pos:..,aa:..)`)・`transl_table`・`codon_start` を読み、initiator を M 強制・末尾 stop 除去。

## File Structure

| ファイル | 責務 | 変更 |
|---|---|---|
| `src/ddbj_gff/aa_names.py` | 略号↔full name↔stop 判定（共有） | Create (Task 1) |
| `src/ddbj_gff/normalize/passes.py` | `_parse_pos_spec` / `pass_transl_except` / `pass_anticodon` | Modify (T1,T2) |
| `src/ddbj_gff/normalize/normalize.py` | `ALL_PASSES`＋`_APPLIED` 更新 | Modify (T1,T2) |
| `src/ddbj_gff/mss/translate.py` | vendored 翻訳関数 | Create (T3) |
| `src/ddbj_gff/mss/convert.py` | `build_cds_feature` 翻訳統合 | Modify (T4) |
| `src/ddbj_gff/validate/rules.py` | `rule_feature_type` に特殊型 accept-set | Modify (T5) |
| `tests/test_normalize_pass_transl_except.py`/`_anticodon.py`/`test_mss_translate.py`/`test_mss_cds.py`(追記)/`test_validate_rules_body.py`(追記)/`test_normalize_integration.py`(追記) | テスト | T1–T6 |

---

### Task 1: aa_names ＋ _parse_pos_spec ＋ pass_transl_except

**Files:**
- Create: `src/ddbj_gff/aa_names.py`
- Modify: `src/ddbj_gff/normalize/passes.py`, `src/ddbj_gff/normalize/normalize.py`
- Test: `tests/test_normalize_pass_transl_except.py`

**Interfaces:**
- Produces: `aa_names.full_name(a)`, `aa_names.to_abbrev(n)`, `aa_names.is_stop(a)`, `aa_names.ABBREV_TO_NAME`. `passes._parse_pos_spec(spec) -> dict|None` (keys start,end,strand,aa,seq). `passes.pass_transl_except(doc, ctx) -> list[Change]`.

- [ ] **Step 1: aa_names.py を作成**

`src/ddbj_gff/aa_names.py`:
```python
"""Amino-acid abbreviation <-> INSDC full-name mappings for special-case canonicalization."""
from __future__ import annotations

ABBREV_TO_NAME: dict[str, str] = {
    "Ala": "alanine", "Arg": "arginine", "Asn": "asparagine", "Asp": "aspartic acid",
    "Cys": "cysteine", "Gln": "glutamine", "Glu": "glutamic acid", "Gly": "glycine",
    "His": "histidine", "Ile": "isoleucine", "Leu": "leucine", "Lys": "lysine",
    "Met": "methionine", "Phe": "phenylalanine", "Pro": "proline", "Ser": "serine",
    "Thr": "threonine", "Trp": "tryptophan", "Tyr": "tyrosine", "Val": "valine",
    "Sec": "selenocysteine", "Pyl": "pyrrolysine",
}
NAME_TO_ABBREV: dict[str, str] = {v: k for k, v in ABBREV_TO_NAME.items()}
_STOP = {"Term", "TERM", "*"}


def full_name(abbrev: str) -> str:
    return ABBREV_TO_NAME.get(abbrev, abbrev)


def to_abbrev(name: str) -> str:
    return NAME_TO_ABBREV.get(name, name)


def is_stop(aa: str) -> bool:
    return aa in _STOP
```

- [ ] **Step 2: 失敗するテストを書く**

`tests/test_normalize_pass_transl_except.py`:
```python
from ddbj_gff.model import Feature, Span, GffDocument
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.normalize.passes import NormalizeContext, pass_transl_except


def _ctx():
    return NormalizeContext(vocab=None, seq_lengths=None, config=NormalizeConfig())


def _cds(attrs):
    # CDS spanning 1..582 on chr1 (+)
    return Feature("c", "S", "CDS", [Span("chr1", 1, 582, "+", 0)], dict(attrs), [])


def _run(f):
    doc = GffDocument(features=[f], feature_index={f.id: f})
    changes = pass_transl_except(doc, _ctx())
    kids = [x for x in doc.features if x.type in ("recoded_codon", "stop_codon")]
    return doc, changes, kids


def test_transl_except_becomes_recoded_codon():
    f = _cds({"transl_except": ["(pos:139..141,aa:Sec)"]})
    doc, changes, kids = _run(f)
    assert len(kids) == 1
    k = kids[0]
    assert k.type == "recoded_codon"
    assert k.spans[0].start == 139 and k.spans[0].end == 141
    assert k.attributes["codon_redefined"] == ["selenocysteine"]
    assert k.attributes["Parent"] == ["c"] and k.parent_ids == ["c"]
    assert k.attributes["ID"] == [k.id]
    assert "transl_except" not in f.attributes
    assert any(c.action == "add-child-feature" for c in changes)
    assert k in f.children and doc.feature_index.get(k.id) is k


def test_complement_pos_sets_minus_strand():
    f = _cds({"transl_except": ["(pos:complement(139..141),aa:Sec)"]})
    _, _, kids = _run(f)
    assert kids[0].spans[0].strand == "-"


def test_stop_aa_becomes_stop_codon():
    f = _cds({"transl_except": ["(pos:580..582,aa:Term)"]})
    _, _, kids = _run(f)
    assert kids[0].type == "stop_codon"


def test_out_of_bounds_is_needs_manual_and_kept():
    f = _cds({"transl_except": ["(pos:9000..9002,aa:Sec)"]})
    doc, changes, kids = _run(f)
    assert kids == []
    assert f.attributes.get("transl_except") == ["(pos:9000..9002,aa:Sec)"]
    assert any(c.action == "needs-manual" for c in changes)


def test_no_transl_except_is_noop():
    f = _cds({})
    doc, changes, kids = _run(f)
    assert changes == [] and kids == []
```

- [ ] **Step 3: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_pass_transl_except.py -v` → FAIL（`pass_transl_except` 未定義）。

- [ ] **Step 4: passes.py に実装** — import に `Feature, Span` と `re`、`aa_names` を追加し、末尾に追記:
```python
import re
from ..model import Directive, Feature, Span   # ← Feature, Span を追加
from .. import aa_names                          # ← 追加


def _parse_pos_spec(spec: str) -> dict | None:
    """Parse '(pos:139..141,aa:Sec[,seq:ttc])' (already URL-decoded). Single range only."""
    if "join" in spec.lower():
        return None
    m = re.search(r"pos:(complement\()?(\d+)(?:\.\.(\d+))?\)?", spec)
    if not m:
        return None
    start = int(m.group(2))
    end = int(m.group(3)) if m.group(3) else start
    strand = "-" if m.group(1) else "+"
    aa_m = re.search(r"aa:([A-Za-z*]+)", spec)
    seq_m = re.search(r"seq:([A-Za-z]+)", spec)
    return {"start": start, "end": end, "strand": strand,
            "aa": aa_m.group(1) if aa_m else None,
            "seq": seq_m.group(1) if seq_m else None}


def _attach_children(doc, pending) -> None:
    for child, parent in pending:
        doc.features.append(child)
        if child.id:
            doc.feature_index[child.id] = child
        parent.children.append(child)


def pass_transl_except(doc, ctx) -> list:
    changes: list = []
    pending: list = []
    for f in list(doc.features):
        if f.type != "CDS":
            continue
        specs = f.attributes.get("transl_except")
        if not specs or not f.spans:
            continue
        lo = min(s.start for s in f.spans)
        hi = max(s.end for s in f.spans)
        seqid = f.spans[0].seqid
        kept: list = []
        made = 0
        for spec in specs:
            p = _parse_pos_spec(spec)
            if p is None or p["aa"] is None or not (lo <= p["start"] and p["end"] <= hi):
                changes.append(Change("needs-manual", f.id or "?",
                                      f"CDS {f.id!r} transl_except {spec!r} not a single in-bounds pos; kept as attribute"))
                kept.append(spec)
                continue
            made += 1
            child_id = f"{f.id}_recoded_{made}"
            if aa_names.is_stop(p["aa"]):
                ctype = "stop_codon"
                attrs = {"ID": [child_id], "Parent": [f.id], "Note": ["stop codon completed"]}
            else:
                ctype = "recoded_codon"
                attrs = {"ID": [child_id], "Parent": [f.id],
                         "codon_redefined": [aa_names.full_name(p["aa"])]}
            child = Feature(child_id, f.source, ctype,
                            [Span(seqid, p["start"], p["end"], p["strand"], 0)], attrs, [f.id])
            pending.append((child, f))
            changes.append(Change("add-child-feature", child_id,
                                  f"CDS {f.id!r}: transl_except -> {ctype} ({p['start']}..{p['end']})"))
        if made:
            if kept:
                f.attributes["transl_except"] = kept
            else:
                del f.attributes["transl_except"]
    _attach_children(doc, pending)
    return changes
```

- [ ] **Step 5: normalize.py を更新** — `ALL_PASSES` と `_APPLIED`:
```python
from .passes import NormalizeContext, pass_directives, pass_so_terms, pass_transl_except

ALL_PASSES = [pass_directives, pass_so_terms, pass_transl_except]

_APPLIED = {"add-directive", "rename-type", "add-qualifier", "add-child-feature"}
```

- [ ] **Step 6: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_pass_transl_except.py -v` → 5 passed。全体 `docker exec ddbj-gff-dev uv run pytest -q`（回帰なし）。

- [ ] **Step 7: Commit**
```bash
git add src/ddbj_gff/aa_names.py src/ddbj_gff/normalize/passes.py src/ddbj_gff/normalize/normalize.py tests/test_normalize_pass_transl_except.py
git commit -m "feat(normalize): pass_transl_except -> recoded_codon/stop_codon child features"
```

---

### Task 2: pass_anticodon

**Files:**
- Modify: `src/ddbj_gff/normalize/passes.py`, `src/ddbj_gff/normalize/normalize.py`
- Test: `tests/test_normalize_pass_anticodon.py`

**Interfaces:**
- Consumes: `_parse_pos_spec`, `_attach_children`, `aa_names`（Task 1）。
- Produces: `pass_anticodon(doc, ctx) -> list[Change]`。`ALL_PASSES` に追加。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_normalize_pass_anticodon.py`:
```python
from ddbj_gff.model import Feature, Span, GffDocument
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.normalize.passes import NormalizeContext, pass_anticodon


def _ctx():
    return NormalizeContext(vocab=None, seq_lengths=None, config=NormalizeConfig())


def _run(f):
    doc = GffDocument(features=[f], feature_index={f.id: f})
    changes = pass_anticodon(doc, _ctx())
    kids = [x for x in doc.features if x.type == "anticodon"]
    return doc, changes, kids


def test_anticodon_becomes_child():
    f = Feature("t", "S", "tRNA", [Span("chr1", 14674, 14742, "-")],
                {"anticodon": ["(pos:complement(14710..14712),aa:Glu,seq:ttc)"]}, [])
    doc, changes, kids = _run(f)
    assert len(kids) == 1
    k = kids[0]
    assert k.type == "anticodon"
    assert k.spans[0].start == 14710 and k.spans[0].end == 14712 and k.spans[0].strand == "-"
    assert k.attributes["amino_acid"] == ["glutamic acid"]
    assert k.attributes["sequence"] == ["ttc"]
    assert k.attributes["Parent"] == ["t"] and k.parent_ids == ["t"]
    assert "anticodon" not in f.attributes
    assert any(c.action == "add-child-feature" for c in changes)


def test_no_anticodon_is_noop():
    f = Feature("t", "S", "tRNA", [Span("chr1", 1, 70, "+")], {}, [])
    doc, changes, kids = _run(f)
    assert changes == [] and kids == []
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_pass_anticodon.py -v` → FAIL。

- [ ] **Step 3: passes.py に追記**
```python
def pass_anticodon(doc, ctx) -> list:
    changes: list = []
    pending: list = []
    for f in list(doc.features):
        if f.type != "tRNA":
            continue
        specs = f.attributes.get("anticodon")
        if not specs or not f.spans:
            continue
        lo = min(s.start for s in f.spans)
        hi = max(s.end for s in f.spans)
        seqid = f.spans[0].seqid
        kept: list = []
        made = 0
        for spec in specs:
            p = _parse_pos_spec(spec)
            if p is None or not (lo <= p["start"] and p["end"] <= hi):
                changes.append(Change("needs-manual", f.id or "?",
                                      f"tRNA {f.id!r} anticodon {spec!r} not a single in-bounds pos; kept as attribute"))
                kept.append(spec)
                continue
            made += 1
            child_id = f"{f.id}_anticodon_{made}"
            attrs = {"ID": [child_id], "Parent": [f.id]}
            if p["aa"]:
                attrs["amino_acid"] = [aa_names.full_name(p["aa"])]
            if p["seq"]:
                attrs["sequence"] = [p["seq"]]
            child = Feature(child_id, f.source, "anticodon",
                            [Span(seqid, p["start"], p["end"], p["strand"], None)], attrs, [f.id])
            pending.append((child, f))
            changes.append(Change("add-child-feature", child_id,
                                  f"tRNA {f.id!r}: anticodon -> anticodon child ({p['start']}..{p['end']})"))
        if made:
            if kept:
                f.attributes["anticodon"] = kept
            else:
                del f.attributes["anticodon"]
    _attach_children(doc, pending)
    return changes
```

- [ ] **Step 4: normalize.py を更新**
```python
from .passes import (NormalizeContext, pass_directives, pass_so_terms,
                     pass_transl_except, pass_anticodon)

ALL_PASSES = [pass_directives, pass_so_terms, pass_transl_except, pass_anticodon]
```

- [ ] **Step 5: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_pass_anticodon.py -v` → 2 passed。全体 `docker exec ddbj-gff-dev uv run pytest -q`。

- [ ] **Step 6: Commit**
```bash
git add src/ddbj_gff/normalize/passes.py src/ddbj_gff/normalize/normalize.py tests/test_normalize_pass_anticodon.py
git commit -m "feat(normalize): pass_anticodon -> anticodon child feature"
```

---

### Task 3: 翻訳関数の vendoring

**Files:**
- Create: `src/ddbj_gff/mss/translate.py`
- Test: `tests/test_mss_translate.py`

**Interfaces:**
- Produces: `mss.translate.translate_cds_with_transl_except(feature, parent_seq, stop_symbol="*") -> Seq`。

- [ ] **Step 1: ファイルを取得（host, curl）**
```bash
curl -fsSL https://raw.githubusercontent.com/nigyta/translate_with_exception/d3c382242f1372afb2b49b47a245ba8dcf548cf4/translate_with_transl_except.py -o src/ddbj_gff/mss/translate.py
```
先頭に provenance ヘッダを挿入（既存 import 行の前）:
```python
# Vendored verbatim from https://github.com/nigyta/translate_with_exception
#   commit d3c382242f1372afb2b49b47a245ba8dcf548cf4, file translate_with_transl_except.py
# NIG (National Institute of Genetics) own code, reused with authorization.
# Public API used here: translate_cds_with_transl_except(feature, parent_seq, stop_symbol="*") -> Seq
```

- [ ] **Step 2: smoke テストを書く**

`tests/test_mss_translate.py`:
```python
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from ddbj_gff.mss.translate import translate_cds_with_transl_except


def test_translate_applies_selenocysteine_exception():
    # ATG TGA AAA TAA : codon2 TGA is a stop under table 11, recoded to Sec (U) via transl_except
    parent = Seq("ATGTGAAAATAA")
    feat = SeqFeature(SimpleLocation(0, 12, strand=1), type="CDS",
                      qualifiers={"transl_table": ["11"], "codon_start": ["1"],
                                  "transl_except": ["(pos:4..6,aa:Sec)"]})
    protein = str(translate_cds_with_transl_except(feat, parent))
    assert "U" in protein          # TGA recoded to selenocysteine
    assert "*" not in protein      # no internal stop after applying the exception
    assert protein.startswith("M")
```

- [ ] **Step 3: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_mss_translate.py -v` → 1 passed。（もし pos の解釈で失敗する場合は、vendored 関数の期待入力に合わせ pos を CDS 相対に調整して原因を報告。関数本体は変更しない。）

- [ ] **Step 4: Commit**
```bash
git add src/ddbj_gff/mss/translate.py tests/test_mss_translate.py
git commit -m "vendor(mss): translate_cds_with_transl_except (nigyta/translate_with_exception @ d3c3822)"
```

---

### Task 4: Phase 2 build_cds_feature 翻訳統合

**Files:**
- Modify: `src/ddbj_gff/mss/convert.py`
- Test: `tests/test_mss_cds.py`（追記）

**Interfaces:**
- Consumes: `mss.translate.translate_cds_with_transl_except`（Task 3）、`aa_names`（Task 1）。
- Produces: `_collect_transl_excepts(cds_feat) -> list[str]`。`build_cds_feature` が例外を適用して翻訳。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_mss_cds.py` に追記（既存の import・ヘルパを利用。無ければ以下を先頭に追加: `from ddbj_gff.model import Feature, Span; from Bio.Seq import Seq; from ddbj_gff.mss.config import MssConfig; from ddbj_gff.mss.convert import build_cds_feature`）:
```python
def test_recoded_codon_child_avoids_internal_stop_warning():
    # CDS ATG TGA AAA TAA on + strand; TGA (4..6) is recoded via a recoded_codon child.
    genome = Seq("ATGTGAAAATAA")
    cds = Feature("c1", "S", "CDS", [Span("s", 1, 12, "+", 0)],
                  {"ID": ["c1"], "transl_table": ["11"]}, ["m1"])
    recoded = Feature("c1_recoded_1", "S", "recoded_codon", [Span("s", 4, 6, "+", 0)],
                      {"ID": ["c1_recoded_1"], "Parent": ["c1"], "codon_redefined": ["selenocysteine"]},
                      ["c1"])
    cds.children = [recoded]
    mrna = Feature("m1", "S", "mRNA", [Span("s", 1, 12, "+")], {"ID": ["m1"]}, ["g1"])
    mrna.children = [cds]
    gene = Feature("g1", "S", "gene", [Span("s", 1, 12, "+")], {"ID": ["g1"]}, [])
    cfg = MssConfig(source={"organism": "x", "mol_type": "genomic DNA"})
    diags: list = []
    feat = build_cds_feature(mrna, gene, "LT_1", genome, cfg, diags)
    assert feat is not None
    assert not any(d.code == "translation-internal-stop" for d in diags)
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_mss_cds.py::test_recoded_codon_child_avoids_internal_stop_warning -v` → FAIL（現状 plain translate で TGA が `*` → internal-stop 診断）。

- [ ] **Step 3: convert.py を実装** — ファイル上部（他ヘルパ付近）に追加:
```python
from .translate import translate_cds_with_transl_except
from Bio.SeqFeature import SeqFeature
from .. import aa_names


def _collect_transl_excepts(cds_feat) -> list:
    """Gather transl_except specs from the CDS attribute (raw) and recoded_codon/stop_codon children."""
    specs = list(cds_feat.attributes.get("transl_except", []))
    for child in cds_feat.children:
        if child.type not in ("recoded_codon", "stop_codon"):
            continue
        sp = child.spans[0]
        loc = f"{sp.start}..{sp.end}" if sp.start != sp.end else f"{sp.start}"
        if sp.strand == "-":
            loc = f"complement({loc})"
        if child.type == "stop_codon":
            aa = "Term"
        else:
            aa = aa_names.to_abbrev((child.attributes.get("codon_redefined") or [""])[0])
        specs.append(f"(pos:{loc},aa:{aa})")
    return specs
```
`build_cds_feature` 内の翻訳部（現行 line 153-154 の `coding_full = ...` / `protein = str(Seq(coding_full).translate(table=table_id))`）を以下に置換:
```python
    coding_full = coding[: len(coding) - len(coding) % 3]
    cds_feat = next((c for c in mrna.children if c.type == "CDS"), None)
    excepts = _collect_transl_excepts(cds_feat) if cds_feat is not None else []
    if excepts:
        sf = SeqFeature(cds_feat.to_biopython_location(), type="CDS",
                        qualifiers={"transl_table": [str(table_id)],
                                    "codon_start": [str(codon_start)],
                                    "transl_except": excepts})
        protein = str(translate_cds_with_transl_except(sf, genome_seq))
    else:
        protein = str(Seq(coding_full).translate(table=table_id))
```
（後続の `body`/internal-stop/no-start 判定は不変。例外 CDS は関数が initiator を M 強制するため no-start は発火しないが、これは許容。）

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_mss_cds.py -v` → 追加テスト含め pass。全体 `docker exec ddbj-gff-dev uv run pytest -q`（回帰なし。既存 CDS テストは excepts 空で従来 plain translate 経路を通る）。

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/mss/convert.py tests/test_mss_cds.py
git commit -m "feat(mss): apply transl_except/recoded_codon in CDS translation via vendored function"
```

---

### Task 5: 検証器が INSDC-GFF3 特殊feature型を受入

**Files:**
- Modify: `src/ddbj_gff/validate/rules.py`（`rule_feature_type`）
- Test: `tests/test_validate_rules_body.py`（追記）

**Interfaces:**
- Produces: `rules._INSDC_GFF3_SPECIAL` 定数。`rule_feature_type` が特殊型を `feature-type-not-insdc` から除外。3-B の canonical 子feature（recoded_codon/anticodon/stop_codon）が誤フラグされなくなる。

- [ ] **Step 1: 失敗するテストを書く** — `tests/test_validate_rules_body.py` 末尾に追記（既存 `V`/`codes`/`rules`/`Feature`/`Span`/`GffDocument` を利用）:
```python
def test_insdc_gff3_special_types_not_flagged():
    for t in ("recoded_codon", "anticodon", "stop_codon", "start_codon"):
        f = Feature("x", "S", t, [Span("c", 1, 3, "+")], {}, [])
        c = codes(rules.rule_feature_type(GffDocument(features=[f]), V))
        assert "feature-type-not-insdc" not in c, t
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate_rules_body.py::test_insdc_gff3_special_types_not_flagged -v` → FAIL（recoded_codon 等が V.feature_types に無く発火）。

- [ ] **Step 3: rules.py を修正** — `rule_feature_type` の直前に定数を追加し、条件に `and f.type not in _INSDC_GFF3_SPECIAL` を足す（**メッセージ本文・他は不変**）:
```python
# canonical INSDC GFF3 special-case feature types (spec v0.5): valid col3 values not present
# in the SO-term column of feature-mapping.tsv
_INSDC_GFF3_SPECIAL = {"recoded_codon", "anticodon", "stop_codon", "start_codon"}
```
既存 `rule_feature_type` の判定行を:
```python
        if f.type not in vocab.feature_types and f.type not in _INSDC_GFF3_SPECIAL:
```
に変更（`make_diagnostic("feature-type-not-insdc", …)` の本文はそのまま）。

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate_rules_body.py -v` → 追加テスト含め pass。全体 `docker exec ddbj-gff-dev uv run pytest -q`（回帰なし）。

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/validate/rules.py tests/test_validate_rules_body.py
git commit -m "feat(validate): accept INSDC-GFF3 special feature types (recoded_codon/anticodon/stop_codon/start_codon)"
```

---

### Task 6: ecoli end-to-end 統合（slow）＋ anticodon フィクスチャ

**Files:**
- Test: `tests/test_normalize_integration.py`（追記）

**Interfaces:** Consumes `parse`/`normalize`/`validate`/`convert` or `build_cds_feature`、実 ecoli（存在時）。

- [ ] **Step 1: テストを追記**

`tests/test_normalize_integration.py` 末尾（import は既存の `parse`/`validate`/`normalize`/`NormalizeConfig`/`Severity`/`Path`/`ROOT`/`pytest` を利用）:
```python
def test_anticodon_fixture_becomes_child():
    from ddbj_gff.normalize import normalize as _normalize  # already imported as normalize
    gff = (
        "##gff-version 3\n"
        "##sequence-region chr1 1 20000\n"
        "chr1\tS\ttRNA\t14674\t14742\t.\t-\t.\tID=t1;anticodon=(pos:complement(14710..14712)%2Caa:Glu%2Cseq:ttc)\n"
    )
    norm, _ = normalize(parse(gff), config=NormalizeConfig(taxid=1148))
    anti = [f for f in norm.features if f.type == "anticodon"]
    assert len(anti) == 1
    assert anti[0].attributes["amino_acid"] == ["glutamic acid"]
    # validator no longer flags the anticodon attribute (it was converted)
    codes = {d.code for d in validate(norm)}
    assert "noncanonical-special-case" not in codes


@pytest.mark.slow
def test_ecoli_transl_except_canonicalized_and_translates():
    import gzip
    p = ROOT / "examples" / "ecoli" / "GCF_000005845.2_ASM584v2_genomic.gff.gz"
    if not p.exists():
        pytest.skip(f"missing {p}")
    text = gzip.decompress(p.read_bytes()).decode("ascii", errors="replace")
    norm, _ = normalize(parse(text), config=NormalizeConfig(taxid=511145))
    # transl_except attributes converted to recoded_codon children
    recoded = [f for f in norm.features if f.type == "recoded_codon"]
    assert len(recoded) >= 3   # fdnG / fdoG / fdhF selenocysteine
    diags = validate(norm)
    # transl_except no longer flagged as noncanonical (attribute -> recoded_codon child)
    assert not any(d.code == "noncanonical-special-case" and "transl_except" in (d.message or "")
                   for d in diags)
    # canonical recoded_codon children are accepted by the validator (Task 5), not flagged
    assert not any("recoded_codon" in (d.message or "") for d in diags)
```

- [ ] **Step 2: 実行**
`docker exec ddbj-gff-dev uv run pytest tests/test_normalize_integration.py -v`（非slow の anticodon フィクスチャ pass）。
`docker exec ddbj-gff-dev uv run pytest tests/test_normalize_integration.py -m slow -v`（ecoli: recoded_codon 子が3件以上、transl_except 由来 noncanonical 消滅）。想定外なら実データ調査（パス実装は変えない）。

- [ ] **Step 3: 全体確認** — `docker exec ddbj-gff-dev uv run pytest -q`（slow 除外）→ 全 pass。

- [ ] **Step 4: Commit**
```bash
git add tests/test_normalize_integration.py
git commit -m "test(normalize): anticodon fixture + ecoli transl_except end-to-end (slow)"
```

---

## Self-Review

**1. Spec coverage**（spec §→タスク）:
- §3 pass_transl_except（recoded_codon/stop_codon 子）→ Task 1。§3 pass_anticodon → Task 2。ALL_PASSES/_APPLIED → Task 1,2。
- §4 vendoring → Task 3。§4 build_cds_feature 統合（属性＋子から例外収集・関数呼出）→ Task 4。
- §R-D7 aa マッピング集約（aa_names 共有）→ Task 1（3-B）＋ Task 4（Phase2 逆引き）。
- **3-B canonical 子feature を 3-A が受入**（ユーザー決定・spec 追補）→ Task 5（rule_feature_type accept-set）。
- §5 テスト（transl_except/anticodon/translate smoke/build_cds/validator accept/ecoli slow）→ Task 1-6。

ギャップ/意図的事項:
- no-start 判定: 例外 CDS は vendored 関数が M 強制のため発火しないが許容（spec §4 の「M 強制に整合」）。非例外 CDS は従来経路で不変（回帰なし）。
- 多exon跨ぎ codon（同ID 2行）は単一range前提のため未対応 → needs-manual（spec R-D4）。稀ケース、次以降で扱う。

**2. Placeholder scan**: 各ステップに実コード。"TBD"等なし。GFF フィクスチャは TAB（`\t`）＋ transl_except/anticodon 値の `,` は `%2C` エンコード（Phase1 が decode）。

**3. Type consistency**: `aa_names.full_name/to_abbrev/is_stop`（T1）、`_parse_pos_spec(spec)->dict|None`・`_attach_children`（T1, T2 で共有）、`pass_transl_except`/`pass_anticodon(doc,ctx)->list[Change]`、`Change("add-child-feature"/"needs-manual",…)`、`ALL_PASSES`（4件）、`_APPLIED`（+add-child-feature）、`_collect_transl_excepts(cds_feat)->list`（T4）、`translate_cds_with_transl_except(feature,parent_seq)`（T3）— 一致。`Feature(id,source,type,spans,attributes,parent_ids)`・`Span(seqid,start,end,strand,phase)` は既存に一致。子 attributes に ID/Parent を必ず設定（writer 要件）。
