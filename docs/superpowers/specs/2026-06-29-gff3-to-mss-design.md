# 設計書: GFF3 → DDBJ 登録形式 (MSS) 変換（フェーズ2 / MVP）

- 日付: 2026-06-29
- 対象: `ddbj_gff` ライブラリ上に、GFF3 + ゲノム FASTA を DDBJ MSS 形式（`.ann` 注釈 + 配列 FASTA）へ変換するサブパッケージ `ddbj_gff.mss` を実装する（minimal 相当: 遺伝子座あたり先頭 transcript のみ）。
- 位置づけ: フェーズ2。フェーズ1（GFF3 パーサ＋オブジェクトモデル）を土台に利用する。

---

## 1. 背景とゴール

### 1.1 全体の中での位置
`docs/project_goal.txt` の主目的の一つ「GFF から DDBJ への登録」。実務上の当面の主目的。
フェーズ1で作った `ddbj_gff`（`parse` / `GffDocument` / `Feature` / `to_biopython_location` / `errors`）を土台に、GFF3 → MSS 変換を載せる。

参照:
- DDBJ MSS（Mass Submission System）注釈形式: タブ区切り5列（Entry / Feature / Location / Qualifier_key / Qualifier_value）
- 既存の実験スクリプト `…/ddbj_submission_hifi/ddbj_gff/experimental/gff2mss_for_MP_minimum.py`（bcbio-gff 依存。本実装で置換）と、その出力 `SAMD00647143_marchantia_minimal.ann`、`common.metadata.tsv`

### 1.2 MVP のゴール
- 入力 GFF3 + ゲノム FASTA + `config.toml` + `common.metadata.tsv` から、MSS `.ann` と配列 FASTA を生成する。
- minimal 相当: gene → 代表 transcript（先頭 mRNA）の mRNA + CDS のみ。source、assembly_gap、locus_tag、partial、phase 由来 codon_start、翻訳検証（報告）、COMMON 逐語。
- 主対象データ: marchantia 核ゲノム（marpolbase GFF）。

---

## 2. 確定した方針（意思決定ログ）

| # | 論点 | 決定 |
|---|---|---|
| P2-D1 | MVP スコープ | **minimal 相当**（遺伝子座あたり先頭 transcript のみ）。複数 transcript / UTR / ncRNA 等は後続 |
| P2-D2 | 設定入力 | **`common.metadata.tsv`（COMMON 逐語）+ `config.toml`（source qualifier・locus_tag・transl_table・gap・ff_definition）** の2ファイル |
| P2-D3 | locus_tag | GFF の `locus_tag` 属性を優先。無ければ **`{PREFIX}_{NNNNNN}` 連番**（prefix/width/start/step は config）。元 ID は note に保存 |
| P2-D4 | 成功基準 | **ゴールドスナップショット（小フィクスチャの .ann/.fasta バイト一致）+ 構造アサーション + CDS 翻訳検証** |
| P2-D5 | 変換アーキテクチャ | **案A: 中間 MSS モデル + エミッタ**（converter: GffDocument→MssDocument、emitter: MssDocument→.ann） |
| P2-D6 | パッケージ配置 | フェーズ1 `ddbj_gff` 内の新サブパッケージ `ddbj_gff.mss` |
| P2-D7 | 複数 transcript | 先頭のみ採用し WARNING（残りは破棄） |
| P2-D8 | 翻訳検証 | 報告のみ（emit はブロックしない）。`strict=True` で ERROR 時に送出 |
| P2-D9 | feature 順序 | `source` → `assembly_gap`（位置順）→ gene の `mRNA`+`CDS`（位置順）。既存の受理実績ある構造に合わせる |
| P2-D10 | codon_start | フェーズ1の phase から算出（先頭セグメント phase+1）。従来の固定値 1 を改善 |

依存: フェーズ1 `ddbj_gff`、Biopython（`SeqIO` / `_insdc_location_string` / `Seq.translate` / `BeforePosition`・`AfterPosition`）、stdlib `tomllib`。開発・テストは amd64 Ubuntu コンテナ `ddbj-gff-dev` 内 `uv run pytest`。

---

## 3. アーキテクチャ & モジュール構成

```
src/ddbj_gff/mss/
├── __init__.py        # 公開API（convert, emit_ann, load_config 等）の再エクスポート
├── config.py          # config.toml 読込（型付き設定）, common.metadata.tsv 検証・逐語読込
├── model.py           # MSS データ構造（MssQualifier / MssFeature / MssEntry / MssDocument）
├── locus_tag.py       # locus_tag 解決（GFF属性優先→連番フォールバック）
├── gaps.py            # assembly_gap 検出（N連続→gap feature）
├── convert.py         # converter: GffDocument + seqs + config → MssDocument（+診断）
├── emit.py            # emitter: MssDocument → .ann テキスト / 配列 FASTA
└── cli.py             # CLI: gff/fasta/config/common/out → .ann + .fasta
```

**入力**: GFF3 ファイル / ゲノム FASTA / `config.toml` / `common.metadata.tsv`
**出力**: `<prefix>.ann` / `<prefix>.fasta`（エントリ毎に `//` 区切り）

**データフロー**: `ddbj_gff.parse(gff)` → `GffDocument`、`SeqIO` で FASTA をインメモリ辞書 `{seqid: Seq}` 化 → `convert(doc, seqs, config, common_rows)` → `(MssDocument, diagnostics)` → `emit_ann` / `emit_fasta`。

**CLI 例**: `python -m ddbj_gff.mss --gff X.gff --fasta X.fa --config config.toml --common common.metadata.tsv --out PREFIX`

---

## 4. 設定入力

### 4.1 `common.metadata.tsv`
DDBJ 標準の COMMON 行（DBLINK / SUBMITTER / REFERENCE / ST_COMMENT）。**逐語で読み込み `.ann` 先頭にそのまま出力**（パース・再生成しない）。検証は「存在・非空・先頭が `COMMON`」のみ。

### 4.2 `config.toml`
```toml
[source]                       # 各 source feature に付与する qualifier（記載順を保持）
organism        = "Marchantia polymorpha subsp. ruderalis"
sub_species     = "ruderalis"
strain          = "Tak-1"
mol_type        = "genomic DNA"
sex             = "male"
collection_date = "missing: lab stock"
country         = "missing: lab stock"

[source.chromosome]            # seqid → /chromosome 導出規則
pattern = "^chr(.+)$"          # 例: chr1 → chromosome "1"。不一致なら /submitter_seqid = seqid

[source.ff_definition]         # 一致可否でテンプレ選択（@@[...]@@ は MSS プレースホルダ）
chromosome = "@@[organism]@@ @@[strain]@@ DNA, chromosome: @@[chromosome]@@"
default    = "@@[organism]@@ @@[strain]@@ DNA, @@[entry]@@"

[locus_tag]
prefix = "MPTK1"
width  = 6                     # ゼロ詰め桁 → MPTK1_000010
start  = 10
step   = 10

[cds]
transl_table = 1               # 既定の翻訳表（NCBI table id）。CDS個別指定があれば優先

[assembly_gap]
min_length       = 10
gap_type         = "within scaffold"
linkage_evidence = "align genus"
estimated_length = "known"     # or "unknown"

[product]
default = "hypothetical protein"   # product 無し時。gene があれば "protein {gene}"
```
- `[source]` のキーは source qualifier として記載順に出力。`chromosome`/`ff_definition` は convert が末尾に追加。
- 読込は `tomllib`。未知キーは無視せず WARNING 診断。
- 欠落時の既定値: width=6, start=10, step=10, transl_table=1, min_length=10, estimated_length="known", product.default="hypothetical protein"。`[source]`・`[locus_tag].prefix`・`common.metadata.tsv` は必須。

---

## 5. MSS オブジェクトモデル + エミッタ

### 5.1 モデル（`mss/model.py`、データ構造のみ）
```
MssQualifier(key: str, value: str)
MssFeature(key: str, location: str, qualifiers: list[MssQualifier])   # location は INSDC文字列
MssEntry(name: str, features: list[MssFeature])                        # name = 配列名（例 chr1）
MssDocument(common_rows: list[str], entries: list[MssEntry])           # common_rows は逐語テキスト行
```

### 5.2 エミッタ（`mss/emit.py`）— 5列 .ann の桁規則
- **Entry名(列1)** はそのエントリの最初の行のみ、以降空。
- **Feature(列2)・Location(列3)** はその feature の最初の行のみ、以降の qualifier 行は空。
- qualifier は1つ1行。feature の先頭 qualifier が feature 行に同居。

```python
def emit_ann(doc: MssDocument) -> str:
    lines = list(doc.common_rows)
    for entry in doc.entries:
        first_row_of_entry = True
        for feat in entry.features:
            quals = feat.qualifiers or [MssQualifier("", "")]
            for i, q in enumerate(quals):
                col1 = entry.name if first_row_of_entry else ""
                col2 = feat.key if i == 0 else ""
                col3 = feat.location if i == 0 else ""
                lines.append("\t".join([col1, col2, col3, q.key, q.value]))
                first_row_of_entry = False
    return "\n".join(lines) + "\n"
```

### 5.3 配列 FASTA（`mss/emit.py`）
各エントリを FASTA（60桁折返し）で出力し、末尾に `//` 行を付す（MSS 配列ファイル形式、既存スクリプト準拠）。

---

## 6. コンバータ（`mss/convert.py`）

`convert(doc, seqs, config, common_rows) → (MssDocument, list[Diagnostic])`。配列(seqid)ごとに `MssEntry` を構築。

### 6.1 エントリ内順序
`source` → `assembly_gap`（位置順）→ 各 gene の `mRNA`+`CDS`（位置順）。

### 6.2 span の集約（重要）
GFF により CDS/exon は2通りの表現:
- marchantia: 各セグメントが別ID → 別 Feature（1 span ずつ）で mRNA の子
- NCBI流: 同一ID の 1 Feature に複数 span

どちらも `spans = [s for child in mrna.children if child.type == T for s in child.spans]`（T = "CDS" / "exon"）で **span を集めて union** することで統一的に扱う。strand に応じて整序（+鎖: start 昇順 / −鎖: start 降順）。

### 6.3 各要素
- **source**: location `1..{len(seq)}`。qualifier = config `[source]`（順序保持）+ chromosome（`pattern` 一致なら capture を `/chromosome`、不一致なら `/submitter_seqid = seqid`）+ `ff_definition`（テンプレ選択）。
- **assembly_gap**（`gaps.py`）: `seq.lower()` の `n{min_length,}` を検出。各 N連続 → location `{start}..{end}`（1-based）、qualifier `estimated_length`/`gap_type`/`linkage_evidence`。
- **gene→transcript（minimal）**: gene = root かつ type=="gene"、位置順。代表 transcript = 子 mRNA のうち ID が `.1` 終端のもの、無ければ最初の mRNA 子。複数 mRNA があれば `multi-transcript` WARNING を出し残りを破棄。mRNA が無ければ `no-transcript` WARNING で gene を skip。
- **mRNA**: location = exon span union の INSDC文字列。partial 判定: exon と CDS の端を比較し、exon端 == CDS端（UTR 無し）の端を partial に。strand を考慮して `<`/`>`。qualifier: `locus_tag`, `gene`(あれば), note(submitter id)。exon が無い場合は CDS span で代用し WARNING。
- **CDS**: location = CDS span union の INSDC文字列。CDS が無ければ `no-cds` WARNING で CDS skip。
  - **codon_start**: 先頭（strand順）CDS span の phase から `phase + 1`。phase 不明なら 1。
  - **partial**: codon_start 補正後に配列抽出し、先頭コドン∉開始コドン集合 → 5'partial、末尾コドン∉終止コドン集合 → 3'partial。strand を考慮して `<`/`>`。
  - qualifier: `locus_tag`, `transl_table`(CDS feature の `transl_table` 属性優先・無ければ config 既定), `codon_start`, `product`(GFF の product / 無ければ config 既定。gene 名ありで product 既定なら `protein {gene}`), `gene`(あれば), `inference`(あれば), note(submitter id)。
- **locus_tag**（`locus_tag.py`）: gene の `locus_tag` 属性があればそれ。無ければ `{prefix}_{n:0{width}d}`（n = start から step 刻み、gene 毎に増分）。
- **note**: `submitter_gene_id: {gene.id}, submitter_transcript_id: {mrna.id}`。

### 6.4 partial の表現
Biopython `BeforePosition`/`AfterPosition` を span に適用し、`_insdc_location_string(compound_location, len(seq))` で `<`/`>` 付き INSDC 文字列を生成。`to_biopython_location()` で得た location を基に、strand と partial 端に応じて先頭/末尾 span を置換（既存 `fix_partial_location` 同等のロジック）。

### 6.5 qualifier 整理
GFF の構造用属性（`ID`/`Parent`/`Name`/`part`/`gene_biotype`/`gbkey` 等）は出力しない。キー変換: `Dbxref→db_xref`, `Note→note`。`product`/`inference`/`gene` は保持。

### 6.6 翻訳検証（診断）
codon_start 補正後の CDS 配列を `transl_table` で翻訳し、(a) 長さが3の倍数、(b) 内部終止コドン無し（末尾の終止を除く）、(c) 非 5'partial なら開始コドン（M）で始まる、を確認。違反は `Diagnostic`（WARNING/ERROR）に収集。emit はブロックしない。

---

## 7. エラー処理 / 診断

- フェーズ1 `ddbj_gff.errors`（`Severity` / `Diagnostic` / `GffParseError`）を再利用。
- converter は診断を集約し `(MssDocument, list[Diagnostic])` を返す。既定 lenient。`strict=True` で最初の ERROR にて `GffParseError` を送出。
- 診断コード: `unknown-config-key`, `no-transcript`, `no-cds`, `multi-transcript`, `no-exon`(CDSで代用), `missing-sequence`(GFF の seqid が FASTA に無→ERROR・entry skip), `translation-internal-stop`, `translation-not-multiple-of-3`, `translation-no-start`。
- CLI は重大度別件数を stderr に要約。任意で診断レポートファイル出力。

---

## 8. テスト戦略（TDD・pytest）

1. `test_mss_config.py` — config.toml 読込（必須欠落エラー・既定値補完・未知キー警告）、common.tsv 検証。
2. `test_mss_model.py` — MSS dataclass の基本。
3. `test_mss_emit.py` — 5列桁規則の**ゴールドスナップショット**（手組み MssDocument → 期待 .ann 文字列）、FASTA の `//` と60桁折返し。
4. `test_mss_locus_tag.py` — GFF 属性優先／フォールバック連番（prefix/width/start/step、複数 gene で増分）。
5. `test_mss_gaps.py` — N連続検出（min_length 境界・複数・小文字 n・gap 無し・座標 1-based）。
6. `test_mss_convert.py` — 構造アサーション: source qualifier + chromosome/ff_definition、mRNA/CDS の location 文字列（join / complement）、**span union を両表現で**（marchantia式 複数Feature CDS ＋ NCBI式 複数span CDS）、partial（`<`/`>`：UTR 欠落・開始/終止欠落）、phase 由来 codon_start、product 既定（`protein {gene}`）、note、翻訳検証診断、診断コード（no-transcript / missing-sequence 等）。
7. `test_mss_snapshot.py`（E2E 回帰）— 小さな自作 GFF + FASTA → `.ann`/`.fasta` を生成し、コミット済み期待ファイルと**バイト一致**（一度生成・レビュー・固定）。
8. `test_mss_integration.py`（slow・存在時のみ）— 実 marchantia marpolbase GFF + ゲノム FASTA（gitignore 対象）で変換し、ERROR 診断ゼロ・既知遺伝子の location/locus_tag をスポット確認・性能スモーク。

**フィクスチャ**（`tests/mss_fixtures/`）:
- `mini.gff3` — +鎖 多exon UTR付き完全CDS / −鎖 UTR欠落 partial / NCBI式 単一ID複数span CDS / N連続を含む contig / `locus_tag` 属性ありの gene
- `mini.fa` — 整合する小 FASTA（開始 ATG・終止コドンを正しく配置し CDS が綺麗に翻訳されるよう手組み）
- `config.toml` / `common.metadata.tsv` — 小さな設定
- `expected.ann` / `expected.fasta` — ゴールドスナップショット（converter 完成後に生成・レビュー・コミット）

---

## 9. スコープ境界

### 9.1 内（フェーズ2 MVP）
GFF3 + FASTA + config + common → MSS `.ann` + 配列 FASTA。source / assembly_gap / gene→先頭 mRNA+CDS / locus_tag（属性・連番フォールバック）/ partial 判定 / phase 由来 codon_start / 翻訳検証（報告）/ COMMON 逐語 / CLI。

### 9.2 外（後続フェーズ / 後日）
- 複数 transcript・重複 CDS 戦略（nonredundant / redundant_as_misc / full）
- UTR feature、ncRNA / tRNA / rRNA / miRNA / precursor_RNA
- オルガネラ特殊ケース（trans_splicing / transl_except / 環状座標）
- 厳密な INSDC バリデーション（フェーズ3）
- AGAT 前処理（外部・手動。フェーズ2 の入力前処理として推奨だが本ツールには組み込まない。trans-splicing 等を壊し得るため一律適用しない）
- DDBJ MSS バリデータ連携
- pseudogene、partial の `start_range`/`end_range` 表現

---

## 10. 既存資産（参考）
`…/ddbj_submission_hifi/ddbj_gff/experimental/gff2mss_for_MP_minimum.py` が minimal 変換の実装参考。bcbio-gff 依存・codon_start 固定値1・source/locus_tag ハードコードを、本実装ではフェーズ1ライブラリ + config 外部化 + phase 由来 codon_start で作り直す。出力 `SAMD00647143_marchantia_minimal.ann` が MSS 形式の実例。
