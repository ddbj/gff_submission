# 設計書: INSDC 正規化器（フェーズ3 / サブプロジェクト 3B・共通ケース MVP）

- 日付: 2026-06-30
- 対象: 任意の GFF3（Phase1 オブジェクトモデル）を INSDC GFF3 形式へ近づける **GFF→GFF 正規化器** サブパッケージ `ddbj_gff.normalize`。本イテレーションは**共通ケース MVP**（ディレクティブ補完 ＋ SO-term 正規化）。特殊ケース canonical 化・環状座標は 3B-full（次イテレーション）。
- 前提: Phase1（`ddbj_gff` パーサ/モデル/ライター/`errors`）、Phase2（`ddbj_gff.mss`）、Phase3-A（`ddbj_gff.validate`）が `main` にある。

---

## 1. 背景・全体アーキテクチャ

Phase 3（INSDC GFF3 の正規化＋検証）の分離方針（ユーザー指針 — 生成と前処理を分けて実装を単純化）:
- **3A 検証（実装済み）**: detect-only チェッカ `ddbj_gff.validate`。
- **3B 正規化（本書）**: 任意 GFF → INSDC 整形済み GFF への GFF→GFF 変換。3A の検出と 1:1 対応する修正を行う。
- **Phase 2 生成（実装済み）**: 整形済み GFF → DDBJ MSS。

パイプライン: `任意GFF → [3B 正規化] → INSDC GFF → [3A 検証] → [Phase 2 生成]`。

3A 検出ルールと 3B 修正パスは 1:1 対応する。本 MVP は 3A の `missing-gff-version` / `missing-insdc-gff-version` / `missing-species-taxid` / `missing-sequence-region` / `cds-missing-transl-table` / `feature-type-not-insdc` に対応する修正を担う。

参照: INSDC GFF3 Specification、gff3tools `feature-mapping.tsv`（3A で同梱済み。列2=SO term、列4=INSDC Feature、列5–6=Qualifier）。

---

## 2. 確定した方針（意思決定ログ）

| # | 論点 | 決定 |
|---|---|---|
| N-D1 | スコープ | 共通ケース MVP: ディレクティブ補完 ＋ SO-term 正規化。特殊ケース canonical 化・環状座標は 3B-full へ |
| N-D2 | GFF に無い情報の取得 | **FASTA 任意 ＋ config**。`##sequence-region` は FASTA があれば真の長さ、無ければ feature 最大座標から近似（要警告）。`#!transl_table` は一貫した既存 CDS 属性を昇格、無ければ config 既定（1）。`##species` は config の taxid から（未指定なら補完せず報告） |
| N-D3 | 出力の形 | クリーンな INSDC 出力 GFF ＋ **構造化された変更レポート**（適用済み / 未解決を分離）。出力 GFF に provenance 属性は混ぜない |
| N-D4 | アーキテクチャ | 案A: パス・レジストリ（3A の `ALL_RULES` と対称な `ALL_PASSES`）。各パス `(doc, ctx) -> list[Change]`。Phase1 パーサ/ライター/モデルを再利用 |
| N-D5 | 変更モデル | `normalize()` は入力 doc を**変更しない**。作業コピーを変換し `(正規化doc, report)` を返す。冪等（正規化済みの再正規化は no-op） |
| N-D6 | qualifier 自動付与 | 具体値の qualifier のみ自動付与。プレースホルダ（`<…>`）・末尾 `*` を含む値は**捏造せず** report に「手動で値が必要」と記録 |
| N-D7 | exit code | 変換成功なら 0（未解決は report で通知。3B は検証器ではない）。parse 失敗・ファイル不在等のみ非0 |
| N-D8 | 範囲外 | 特殊ケース canonical 化、環状座標、プレースホルダ qualifier の自動補完、feature 削除/統合、深い生物学的補正、検証（3A）・MSS 生成（Phase2） |

依存: Phase1 `ddbj_gff`（`parse`/`GffDocument`/`Feature`/`Span`/writer/`errors`）、Phase3-A `ddbj_gff.validate`（`load_vocab`／round-trip 検証）。`normalize` 本体は stdlib のみ（FASTA 読込は CLI 側で Biopython）。開発・テストは amd64 コンテナ `ddbj-gff-dev` 内 `uv run pytest`。

---

## 3. アーキテクチャ・モジュール構成

```
src/ddbj_gff/normalize/
├── __init__.py     # normalize, NormalizationReport, NormalizeConfig 再エクスポート
├── normalize.py    # normalize(doc, *, seq_lengths=None, config=None) -> (GffDocument, NormalizationReport); ALL_PASSES
├── passes.py       # pass_directives, pass_so_terms（各 (doc, ctx) -> list[Change]）
├── report.py       # Change / NormalizationReport（dataclass）
├── config.py       # NormalizeConfig（taxid / transl_table / insdc_gff_version）＋ TOML ロード
├── cli.py
└── __main__.py
```

- データフロー: `parse(GFF) → normalize(doc, …) → (正規化doc, report) → writer で GFF 書き出し ＋ report 出力`。
- `ALL_PASSES` は 3A の `ALL_RULES` と 1:1 対応（本 MVP は `pass_directives`・`pass_so_terms` の2つ。3B-full で特殊ケース・環状パスを追加）。
- `NormalizeContext`（内部）: `vocab`（拡張 Vocab）・`seq_lengths`・`config` を各パスへ渡す。

**vocab 拡張（3-A の `vocab.py` を流用・拡張）:** 現状 `insdc_map` は列4 のみ。3B は列5–6（Qualifier）も要るので `Vocab` に `feature_qualifiers: dict[str, tuple[str, ...]]` を**追加**（既存 `insdc_map`/`feature_types`/`dbxref_dbtags` は不変）。`load_vocab()` の lru_cache を共有。
**重複 SO-term の解決:** 同一 SO-term が複数行ある場合（`LINE_element`/`SINE_element`/`gap`/`guide_RNA`/`mobile_genetic_element`）、**具体値（`<…>`・`*` を含まない）の行を優先**。残る矛盾は最初の行を採用し dedup 警告を記録。

---

## 4. パス1: ディレクティブ補完（`pass_directives`）

欠けている必須ディレクティブを補う。各々「既にあれば no-op（冪等）」。

| 補完するディレクティブ | 値の決め方 | 対応 3A コード |
|---|---|---|
| `##gff-version 3` | 固定（先頭に挿入） | missing-gff-version |
| `#!insdc-gff-version <V>` | config `insdc_gff_version`（既定 `1.0.0`） | missing-insdc-gff-version |
| `##species …?id=<taxid>` | config の `taxid`（任意）。未指定なら補完せず report に記録 | missing-species-taxid |
| `##sequence-region <seqid> 1 <len>` | feature があり region 無い seqid ごと。len は FASTA があれば真の長さ、無ければ feature 最大 end の近似＋警告 | missing-sequence-region |
| `#!transl_table primary:N` | file 先頭に既存無く transl_table を持たない CDS がある場合に追加。N=一貫した既存 CDS 値、無ければ config `transl_table`（既定 1） | cds-missing-transl-table |

- `##species` は外部情報（taxid）が要るため、config 未指定なら誤推測せず未補完として report に明記。
- `##sequence-region` の近似 fallback は、3A `feature-outside-region` を短い region で誤って消さないよう「近似」を必ず report に出す。

---

## 5. パス2: SO-term 正規化（`pass_so_terms`）

各 feature の `type` を `insdc_map` で引く:

1. **マップあり & 異なる**（`pseudogenic_CDS→CDS`, `coding_exon→exon`, `processed_pseudogene→gene`, `ncRNA_gene→ncRNA` 等）→ `f.type` を INSDC feature 名へ書き換え、qualifier を GFF 属性として付与:
   - `/pseudo`（フラグ）→ 属性 `pseudo`（値なし）
   - `/key="value"`（具体値）→ 属性 `key=value`（先頭 `/`・引用符を除去。例 `ncRNA_class=miRNA`, `pseudogene=processed`, `mod_base=ac4c`）
   - **プレースホルダ（`<…>`）・末尾 `*` を含む値は付与しない** → report に「feature X: `/operon` は手動で値が必要」を記録。
2. **マップあり & 同一**（`CDS→CDS` 等、既に INSDC 名）→ no-op。
3. **マップなし**（既知 SO-term でない type）→ 変更せず report に「未マップ type（手動対応）」。3A `feature-type-not-insdc`(WARN) と対応。

Phase 2 との整合: `miRNA→ncRNA`＋`ncRNA_class=miRNA` のように Phase 2 が期待する属性へ揃うため normalize→generate が滑らかに繋がる。

---

## 6. Change/レポート・API・CLI

**`report.py`:**
```python
@dataclass(frozen=True)
class Change:
    action: str    # "add-directive" | "rename-type" | "add-qualifier" | "unmapped-type"
                   # | "needs-manual" | "approx-region" | "no-taxid" | "dup-mapping"
    target: str    # feature id / seqid / directive 名
    message: str

@dataclass
class NormalizationReport:
    applied: list[Change]      # 実際に変換したもの
    unresolved: list[Change]   # 人手対応が要るもの
    # text レンダリング・件数サマリのヘルパ
```

**API（`normalize.py`、純粋・非破壊）:**
```python
def normalize(doc, *, seq_lengths=None, config=None) -> tuple[GffDocument, NormalizationReport]
```
- `seq_lengths`: `dict[str,int] | None`（seqid→長さ）。FASTA 読込は CLI 側（Biopython）で行い dict を渡す → `normalize` 本体は stdlib のみ。
- `config`: `NormalizeConfig(taxid=None, transl_table=1, insdc_gff_version="1.0.0")`。
- 作業コピーに `ALL_PASSES` を順次適用し `(正規化doc, report)` を返す。

**CLI（`python -m ddbj_gff.normalize`）:**
```
--gff IN.gff            （必須）
--fasta IN.fa           （任意; ##sequence-region の真の長さ源）
--config config.toml    （任意; [normalize] taxid / transl_table / insdc_gff_version）
--taxid / --transl-table / --insdc-gff-version   （任意; config を上書き）
--out OUT.gff           （省略時 stdout）
--report report.txt     （省略時 stderr）
```
- 流れ: parse → (FASTA→seq_lengths) → normalize → writer で正規化 GFF を `--out`、report を `--report` へ。
- exit code: 変換成功なら 0。parse 失敗・ファイル不在等のみ非0。
- config TOML は Phase 2 と同じ TOML スタイル（`[normalize]` セクション）。

---

## 7. テスト戦略・スコープ境界

**テスト（TDD・pytest・コンテナ内）:**
1. `test_normalize_vocab` — 拡張 vocab: `feature_qualifiers` ロード、重複 SO-term dedup（`LINE_element` は `"LINE"` を `"LINE*"` より優先）、既存 `insdc_map` 不変。
2. `test_normalize_passes` — パス単位（発火/非発火・冪等）。`pass_directives`: 5 ディレクティブの補完と冪等性。`pass_so_terms`: rename＋具体 qualifier 付与、プレースホルダの needs-manual 報告、未マップ type 不変＋報告、同名 no-op。
3. `test_normalize` — 集約: `(doc, report)`／入力 doc 不変／applied・unresolved／**ラウンドトリップ・オラクル: normalize → 3A `validate()` → 狙った ERROR が消える**。
4. `test_normalize_cli` — `--gff`→正規化 GFF＋report、`--fasta`で真の長さ、`--taxid`で species 追加、exit 0、ファイル不在で非0。
5. `test_normalize_integration`（slow）— 実 example（rice_cp 等）で normalize→validate により狙った ERROR が消え、特殊ケース系（`noncanonical-special-case`）は残ることを assert（スコープ境界の明文化）。

**フィクスチャ:** `tests/normalize_fixtures/` に「散らかった入力 GFF」＋小さな golden 正規化出力。構造アサート＋validate オラクルで二重検証。

**スコープ境界:**
- 内（本 MVP）: `pass_directives`／`pass_so_terms`（type rename＋具体 qualifier 付与＋未解決報告）／Change レポート／CLI／vocab qualifier 拡張・dedup。
- 外（3B-full 以降）: 特殊ケース canonical 化（trans-splicing join / transl_except→recoded_codon / anticodon→子feature）／環状座標正規化／プレースホルダ qualifier 自動補完／feature 削除・統合／深い生物学的補正。
- 非責務: 検証（3A）・MSS 生成（Phase2）はしない。
