# 設計書: INSDC プロファイル検証器（フェーズ3 / サブプロジェクト 3A）

- 日付: 2026-06-30
- 対象: GFF3（Phase1 オブジェクトモデル）を INSDC GFF3 プロファイルに対して検証し、診断を返す **detect-only** な検証器サブパッケージ `ddbj_gff.validate`。公式 SO-INSDC feature mapping 等の統制語彙を取り込む（「フル」）。
- 前提: Phase1（`ddbj_gff` パーサ/モデル/`errors`）と Phase2（`ddbj_gff.mss`）が `main` にある。

---

## 1. 背景・全体アーキテクチャ

Phase 3（INSDC GFF3 の正規化＋検証）を、独立したサブプロジェクトに分解する。本書はその最初 **3A: 検証器**。

**システム全体の関心事の分離**（ユーザー指針 — 実装を単純化するため生成と前処理を分ける）:
- **3B 正規化/前処理（次サブプロジェクト）**: 任意 GFF → INSDC 整形済み GFF への GFF→GFF 変換。例: `#!transl_table` 付与、`##sequence-region`/`#!insdc-gff-version` 補完、**INSDC 許可 SO-term への変換**（`feature-mapping.tsv` の親へ正規化）、特殊ケースの canonical 化（trans-splicing/transl_except/anticodon）、環状座標。
- **3A 検証（本書）**: 検出のみ・自動修正しない純粋なチェッカ。
- **Phase 2 登録ファイル生成（実装済み）**: 整形済み GFF → MSS。

パイプライン: `任意GFF → [3B 正規化] → INSDC GFF → [3A 検証] → [Phase 2 生成]`。検証器は 3B の出力検証にも使う。

検出ルールと 3B の修正は1対1対応する（3A が検出、3B が修正）。`feature-mapping.tsv` は 3A（列3許可判定）と 3B（SO→INSDC 親への変換）で共有するため、3A で同梱しておくことが 3B の土台になる。

参照: INSDC GFF3 Specification v0.5、reference 実装 `enasequence/gff3tools`（Apache-2.0、`src/main/resources/` に `feature-mapping.tsv`/`qualifier-mapping.tsv`/`default-rule-severities.properties` 等）。

---

## 2. 確定した方針（意思決定ログ）

| # | 論点 | 決定 |
|---|---|---|
| V-D1 | サブプロジェクト | 3A: INSDC プロファイル検証器（3B 正規化の前に） |
| V-D2 | 検査範囲 | **フル**（構造/構文＋生物学的＋統制語彙）。列3 SO-term は公式 `feature-mapping.tsv` で判定 |
| V-D3 | アーキテクチャ | 案A: ルールレジストリ（独立ルール関数群＋code＋severity） |
| V-D4 | 統制語彙の取得 | gff3tools `feature-mapping.tsv`（Apache-2.0・帰属）＋ INSDC dbxref 語彙を**スナップショット同梱**。`PROVENANCE.md`/`NOTICE` を併置。so.owl は不使用 |
| V-D5 | 修正の有無 | **detect-only**（自動修正なし）。修正は 3B の責務 |
| V-D6 | 重大度 | code→既定 severity 表（gff3tools `default-rule-severities`＋仕様の必須性）。`severity_overrides` で上書き（error/warning/info/off） |
| V-D7 | 特殊ケース形式 | 3A は `noncanonical-special-case`（INFO）検出のみ。canonical 変換は 3B |
| V-D8 | 範囲外 | 自動修正/正規化、特殊ケース canonical 変換、so.owl 完全探索、taxid 実在確認、MSS形式検証、深い翻訳妥当性 |

依存: Phase1 `ddbj_gff`（`parse`/`GffDocument`/`Feature`/`Span`/`errors`）。Biopython 不要。開発・テストは amd64 コンテナ `ddbj-gff-dev` 内 `uv run pytest`。

---

## 3. アーキテクチャ・モジュール構成・同梱データ

```
src/ddbj_gff/validate/
├── __init__.py        # validate, validate_cli の再エクスポート
├── vocab.py           # 同梱TSVロード: allowed_feature_types()/feature_insdc_map()/allowed_dbxref_dbtags()
├── severities.py      # rule code -> 既定 Severity 表
├── rules.py           # 個別ルール関数群（各 (doc, vocab) -> Iterable[Diagnostic]）
├── validate.py        # validate(doc, *, severity_overrides=None) -> list[Diagnostic]
├── cli.py             # CLI: parse -> validate -> レポート -> exit code
├── __main__.py
└── data/
    ├── feature-mapping.tsv      # gff3tools 由来スナップショット（列3許可SO term＋INSDC対応）
    ├── dbxref.tsv               # INSDC dbxref DBTAG 語彙スナップショット
    ├── PROVENANCE.md            # 取得元URL・版/コミット・ライセンス・手動リフレッシュ手順
    └── NOTICE                   # Apache-2.0 帰属
```

- `vocab.py`: data/ の TSV を遅延ロード・キャッシュ。`feature-mapping.tsv` の「SO term」列を許可集合に、「INSDC Feature」列を SO→INSDC 対応に。
- 診断は Phase1 `errors.Diagnostic`/`Severity` を再利用。`Diagnostic.code` に rule code、`message`/`line_no` に文脈。

---

## 4. 規則セット・既定重大度・特殊ケース境界

各ルールは安定 code を持つ。既定 severity は gff3tools `default-rule-severities`＋仕様の必須性に基づき、`severity_overrides` で上書き可。

**ヘッダ/ディレクティブ**: `missing-gff-version`(ERROR) / `missing-insdc-gff-version`(ERROR) / `missing-species-taxid`(ERROR) / `missing-sequence-region`(ERROR) / `duplicate-sequence-region`(ERROR)
**エンコード**: `non-ascii`(ERROR)
**座標/seqid**: `undefined-seqid`(ERROR) / `feature-outside-region`(ERROR、環状 is_circular の起点跨ぎは除外) / `start-gt-end-noncircular`(ERROR)
**列3 / SO term**: `feature-type-not-insdc`(WARN、`feature-mapping.tsv` の許可集合に無い)
**ID / Parent**: `missing-id-with-children`(ERROR) / `duplicate-id-different-type`(ERROR) / `multiple-parents`(ERROR、INSDCは単一) / `dangling-parent`(ERROR)
**CDS**: `cds-missing-transl-table`(ERROR、翻訳産物ありで transl_table も file 先頭 `#!transl_table` も無い) / `cds-invalid-phase`(ERROR、0/1/2でない)
**遺伝子注釈**: `gene-missing-locus-tag`(WARN)
**Dbxref**: `dbxref-unknown-dbtag`(WARN)

**特殊ケース形式（3A/3B 境界）**: 3A はプロファイル（構造＋統制語彙）に集中。trans-splicing/transl_except/anticodon の canonical 適合は 3B の責務。3A は非canonical表現（例: `exception=trans-splicing`＋`part=`、`transl_except=` 属性）を `noncanonical-special-case`(INFO) として検出するに留める。

各検出は 3B の修正に対応: `missing-insdc-gff-version`/`cds-missing-transl-table`→3B が付与、`feature-type-not-insdc`→3B が SO-term 変換、`noncanonical-special-case`→3B が canonical 化。

---

## 5. API・CLI・エラー処理

- `validate(doc: GffDocument, *, severity_overrides: dict[str, str] | None = None) -> list[Diagnostic]`：純粋関数。全ルール実行→`Diagnostic` を (line_no, code) 順に整列して返す。**自動修正なし**（入力 doc 不変）。
- `severity_overrides`: code→`error`/`warning`/`info`/`off`。`off` は当該ルールを無効化。
- CLI `python -m ddbj_gff.validate --gff X.gff [--severity CODE=LEVEL ...]`：
  - Phase1 `parse` で読込→ `doc.diagnostics`（パース診断）＋ `validate(doc)`（プロファイル診断）を統合してレポート。
  - 重大度別件数サマリ → 各診断（code/severity/位置/メッセージ）を stderr/stdout に。
  - **ERROR があれば exit 1**、なければ 0。
- 検証は例外を投げない（lenient）。

---

## 6. テスト戦略・スコープ境界

**テスト（TDD・pytest、コンテナ内）:**
1. `test_validate_vocab` — `feature-mapping.tsv` ロードで許可集合に CDS/mRNA/gene/exon 等が含まれ架空 term は含まれない、dbxref 語彙ロード。
2. `test_validate_rules` — 1ルール1テスト（全ルール）。発火/非発火の両ケースで code＋severity を検証。
3. `test_validate` — 集約・`severity_overrides`（off で消える/error 化）・整列・入力 doc 不変。
4. `test_validate_cli` — 有効 INSDC GFF→rc0、無効→rc1＋期待 code、`--severity` 上書き。
5. `test_validate_integration`（slow・存在時のみ）— 実 example（NCBI流）に具体的な期待診断（全例 `missing-insdc-gff-version` ERROR、rice_cp/ecoli は `noncanonical-special-case` INFO 等）。

**フィクスチャ**（`tests/validate_fixtures/`）: `valid_insdc.gff3`（全必須ディレクティブ＋locus_tag/transl_table/許可SO term の clean 例、ERROR ゼロ）＋ 各ルールを突く小さな `invalid_*.gff3`。

**同梱データ取得**: 計画の最初のタスクで `feature-mapping.tsv`（gff3tools raw, Apache-2.0）と INSDC dbxref 語彙スナップショットを取得し `data/` へ配置＋`PROVENANCE.md`/`NOTICE`。dbxref ルールは WARN のため語彙の不完全性は許容。

**スコープ境界:**
- 内（3A）: detect-only の INSDC プロファイル検証（§4 の全ルール）／severity 上書き／CLI／非canonical特殊ケースの INFO 検出。
- 外（3B 以降）: 自動修正・正規化／特殊ケース canonical 変換／so.owl 完全オントロジー探索／NCBI taxonomy の taxid 実在確認／DDBJ MSS 形式検証／CDS 翻訳妥当性の深い検査（Phase2 にあり）。
