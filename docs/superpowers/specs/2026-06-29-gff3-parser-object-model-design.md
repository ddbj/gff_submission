# 設計書: GFF3 パーサ＋オブジェクトモデル（フェーズ1）

- 日付: 2026-06-29
- 対象: 自作 GFF3 パーサ／オブジェクトモデル／ライター（INSDC GFF3 プロファイル対応）
- 位置づけ: DDBJ GFF 変換ツール群の**土台（フェーズ1）**。後続の全変換方向が依存する。

---

## 1. 背景とゴール

### 1.1 プロジェクト全体（`docs/project_goal.txt`）
GFF から DDBJ への登録、および DDBJ Flat File から GFF への変換ツール群を開発する。
INSDC GFF3 仕様は **交換形式**であり、登録用 GFF とは別物。将来的に DDBJ flat file →
交換用 GFF3、DDBJ DB → 交換用 GFF の出力も想定される。

参照仕様:
- INSDC GFF3 Specification v0.5（`docs/INSDC GFF3 Specification - v0.5.docx`）
- Sequence Ontology GFF3（Lincoln Stein 1.26 をベースとする）

### 1.2 フェーズ分割（全体像）
各フェーズを「仕様→計画→実装」サイクルで進める。一気に最終目標までは行わない。

| フェーズ | 内容 | 主テストデータ |
|---|---|---|
| **1（本書）** | 自作 GFF3 パーサ＋オブジェクトモデル＋ライター。変換意味論は持たず、読込→モデル→書出しの**意味的ラウンドトリップ**を TDD で保証 | 全 example |
| 2 | GFF3 → DDBJ 登録形式（MSS）。実務上の主目的 | marchantia |
| 3 | INSDC GFF3 の正規化＋バリデーション（trans-splicing / transl_except / 環状座標の canonical 化） | rice_cp, ecoli |
| 4 | DDBJ Flat File → GFF3（逆方向） | — |

### 1.3 フェーズ1のゴール
- INSDC/SO GFF3 を読み込み、自作オブジェクトモデルへ変換する。
- オブジェクトモデルを GFF3 テキストへ書き出す。
- **意味的ラウンドトリップ**（後述 §5.3）を品質基準とし TDD で担保する。
- 特殊ケース（trans-splicing、transl_except、環状座標、order）を**データとして忠実に保持**する（解釈・変換はしない）。

---

## 2. 確定した方針（意思決定ログ）

| # | 論点 | 決定 |
|---|---|---|
| D1 | 開始サブプロジェクト | フェーズ1（パーサ＋オブジェクトモデル）から |
| D2 | フェーズ1の範囲 / 品質基準 | パーサ＋モデル＋ライター。**意味的ラウンドトリップ**（バイト一致は狙わない） |
| D3 | CDS phase の扱い | **各 span に phase を保持**し、`codon_start` は派生プロパティ |
| D4 | スケール / 性能 | **インメモリ全体モデル**（arabidopsis 121MB を現実的に処理。定数メモリ・ストリーミングは将来課題） |
| D5 | 依存・環境管理 | 開発は **amd64 ベースの Ubuntu コンテナ**。**uv + pyproject.toml** |
| D6 | オブジェクトモデル構成 | **案A: Feature 中心の集約モデル**（ID 単位の単一 Feature＋複数 span＋親子グラフ） |
| D7 | AGAT による標準化 | **フェーズ2の推奨・任意の外部前処理**。Python 依存に含めない。特殊ケースを壊し得るため一律適用しない（caveat 記載） |
| D8 | `part` と `number` の区別 | `part`=位置断片の順序（構造・Span 単位）／`number`,`exon_number`=生物学的序数（修飾子・Feature 単位）。混同しない |
| D9 | ライターの出力順 | 既定は **初出順**（ラウンドトリップが安定）。INSDC 推奨の正規ソートは**オプション**（正規化は Phase3 の責務） |
| D10 | パーサ厳格度 | 既定 **lenient ＋診断収集**。`strict=True` で例外送出（厳格検証は Phase3） |

その他の既定: Python 3.11+、テストは pytest、`src/` レイアウト、依存は **Biopython のみ**
（bcbio-gff は使用しない＝`project_goal.txt` 準拠）。

---

## 3. アーキテクチャ & プロジェクト構成

### 3.1 モジュール分割（単一責務・独立テスト可能）
- `model.py` — データ構造のみ（`GffDocument`, `Feature`, `Span`, `Directive`）。パース/出力ロジックを持たない
- `parser.py` — GFF3 テキスト → `GffDocument`
- `writer.py` — `GffDocument` → GFF3 テキスト
- `attributes.py` — 列9の属性パース/エスケープ（`%XX` デコード、`,` 区切り、順序保持）
- `errors.py` — 例外（`GffParseError`）と診断（`Diagnostic`、重大度 enum）

### 3.2 ディレクトリ構成
```
gff_submission/                 # 既存リポジトリ root（git 未初期化 → 初期化予定）
├── pyproject.toml              # uv 管理, package = ddbj_gff, dep = biopython のみ
├── uv.lock
├── Dockerfile                  # linux/amd64 ubuntu + uv + Python 3.11
├── .devcontainer/devcontainer.json
├── .gitignore                  # examples の巨大ファイル(fasta / 121MB gff)等を除外
├── src/ddbj_gff/
│   ├── __init__.py
│   ├── model.py
│   ├── parser.py
│   ├── writer.py
│   ├── attributes.py
│   └── errors.py
├── tests/
│   ├── fixtures/               # 小さな自作 GFF 断片（難所を凝縮）
│   ├── test_attributes.py
│   ├── test_model.py
│   ├── test_parser.py
│   ├── test_writer.py
│   ├── test_roundtrip.py
│   └── test_integration.py     # 実ファイル（巨大ファイルは存在時のみ・slow marker）
├── examples/                   # 既存（巨大バイナリは git 管理外）
└── docs/
```

### 3.3 開発環境
- 開発は `linux/amd64` の Ubuntu コンテナ内で行う。`Dockerfile`（ubuntu base + uv + Python 3.11）と `.devcontainer/devcontainer.json` を用意。
- macOS(arm64) 上でも `--platform linux/amd64` で再現可能にする。
- 開発ループ: コンテナ内 `uv run pytest`（既定は高速テスト、`-m slow` で統合/性能）。

### 3.4 巨大ファイル方針
arabidopsis GFF(121MB)・各 fasta(数百MB) は git に入れず `.gitignore`。テストは `tests/fixtures/`
の小さな自作 GFF を主とし、実ファイル統合テストは「存在すれば実行（slow マーカー）」とする。

---

## 4. オブジェクトモデル（案A: Feature 中心の集約モデル）

GFF の「行」固有の情報（座標・strand・phase・score・part）と「feature」固有の情報
（ID・type・source・属性）を分離する。同一 ID の複数行は1つの `Feature` に集約し、各行は `Span` になる。

```
GffDocument
├── directives    : list[Directive]      # 出現順保持
│     便利アクセサ: gff_version / insdc_gff_version / species(taxid)
│                  / sequence_regions{seqid:(start,end)} / transl_table_map
├── features      : list[Feature]        # 出現順（全 feature）
├── feature_index : dict[id, Feature]    # ID 索引
├── roots         : list[Feature]        # Parent を持たない最上位 feature
├── fasta         : dict[id, Seq] | None # ##FASTA 以降のペプチド（INSDC: ペプチドのみ）
├── sequences     : dict[seqid, Seq]|None# 任意で外部核酸 FASTA を attach（Phase1 のパースには不要）
└── diagnostics   : list[Diagnostic]     # パース時の診断

Feature
├── id          : str | None             # ID 属性
├── source      : str                    # 列2（feature 内で一定を期待・不一致は診断）
├── type        : str                    # 列3（SO term）
├── spans       : list[Span]             # 1個以上（multi-exon / discontinuous）
├── attributes  : dict[str, list[str]]   # 列9（順序保持・%XX デコード済・未知タグも保持）
├── parent_ids  : list[str]              # Parent 属性（INSDC=単一 / SO=複数 を一般化）
├── children    : list[Feature]          # グラフ構築後に解決
├── parents     : list[Feature]          # 逆参照
└── 便利プロパティ: name / dbxref / note / locus_tag / gene / product / protein_id
                  / transl_table / is_circular / is_trans_spliced / is_ordered
                  / number / exon_number / codon_start (CDS のみ; 先頭 span の phase から派生)
   メソッド: to_biopython_location() → FeatureLocation / CompoundLocation を生成
            （part 順 or strand 順に整序。座標切り出し・翻訳は後フェーズで使用）

Span                                      # GFF 1行ぶんの位置情報
├── seqid  : str                          # 通常 feature 内で一定。trans-splicing では異なり得る
├── start  : int                          # 1-based inclusive（GFF 流をそのまま保持）
├── end    : int
├── strand : str                          # '+' '-' '.' '?'
├── phase  : int | None                   # CDS のみ 0/1/2
├── score  : float | None                 # 列6
└── part   : int | None                   # part= 属性（順序の明示: 環状 / trans-splicing / order）

Directive
├── raw    : str                          # 元の行
├── kind   : str                          # gff-version / sequence-region / species
│                                         #  / insdc-gff-version / transl_table / FASTA / unknown ...
└── value  : (kind ごとに構造化)
```

### 4.1 集約ルール（discontinuous feature）
- 同一 `ID` を持つ複数行 → 1つの `Feature`（spans に集約）。全行は同一 type を期待（不一致は診断）。
- ID 無し行 → それぞれ独立 `Feature`（集約しない）。
- span 順序: `part=` があれば part 順、無ければファイル出現順を保持。
  `to_biopython_location()` 側で「+鎖は列4昇順／−鎖は列5降順、ただし part・trans_splicing・環状が
  ある場合はそれを優先」という SO 規則で整序する。

### 4.2 親子グラフ
全 feature 収集後の2パス目で `parent_ids` を解決し `children`/`parents` をリンク（前方参照に対応）。
`roots` = Parent 無し feature。

### 4.3 属性の3層分類（D8）
全属性は `attributes: dict[str, list[str]]` に**逐語保持**され自動的にラウンドトリップする。
そのうえで認識する一群に便利アクセサを付け、構造属性と序数修飾子を**混同しない**。

- **構造属性（位置組立てに影響）**: `part`（Span 単位・断片順序）、`is_ordered`（Feature 単位・order/join）
- **認識する修飾子（Feature 単位・逐語保持＋アクセサ）**: `number`, `exon_number`, `Name`, `Dbxref`,
  `Note`, `locus_tag`, `gene`, `product`, `protein_id`, `transl_table`, `exception`, `Is_circular` …
- **その他/未知（逐語保持）**: `gbkey`, `gene_biotype`, `sub-species` 等

原則: **`part` だけを構造的順序として扱い、`number` 系は順序付けに使わない**
（同一 feature の断片が同じ `number` を持ち得るため。例: `chloroplast.gff3` の trans-spliced intron-1）。

### 4.4 Biopython 境界
`Span`/`Feature`（1-based・phase/part 付き）が真実の源。座標演算が要る時だけ
`to_biopython_location()` で Biopython の 0-based location へ橋渡しする。配列 I/O は `SeqIO`、配列は `Seq`。
phase は独自に保持し、`codon_start` は派生（codon_start = 先頭 span の phase + 1）。

---

## 5. パース処理フロー・ライター・ラウンドトリップ

### 5.1 パース処理フロー（`parser.py`、2パス・インメモリ）

**Pass 0 — 行分類（行単位ストリーム読み）**
- 空行 / 素の `#` コメント → ドロップ（件数のみ診断）
- `##`・`#!` ディレクティブ → `Directive` 化して順序保持
- feature 行（9列）→ 生レコードへトークナイズ
- `###`（前方参照解決境界）→ Phase1 は全体2パスのため意味的には無視、ディレクティブとして保持
- `##FASTA` → 以降を FASTA モードへ
- ASCII 想定（INSDC）。非 ASCII は WARNING を記録しつつ受理

**Pass 1 — レコード→Feature 集約**
- 列9を `attributes.py` で解析（`;` 分割 → `key=value` → `,` で多値分割 → `%XX` デコード、順序保持）
- 列1,4,5,6,7,8 ＋ `part` から `Span` を生成
- 集約: ID 既出 → 既存 Feature に `Span` 追加（type/source/属性の整合を検査、相違は診断）／
  ID 新規 → Feature 生成＋Span／ID 無し → 無名の単独 Feature
- `feature_index[id]` と `features`（初出順）を更新

**Pass 2 — グラフ解決**
- `parent_ids` を `feature_index` で解決し `children`/`parents` をリンク（前方参照対応）。
  親不在は WARNING（`parent_ids` は保持）
- `roots` = Parent 無し feature。`part` があれば各 feature の span を part 順に整序

**FASTA**: `##FASTA` 以降を `SeqIO` で読み `document.fasta`（INSDC はペプチドのみ）。
核酸配列は別ファイルを `attach_sequences(path)` で後付けできるフック（Phase1 のパース自体には不要）。

### 5.2 ライター（`writer.py`）
- ディレクティブを保持順に出力
- feature は `document.features`（初出順）で走査 → ラウンドトリップが自然に安定。
  INSDC 推奨の正規ソート（seqid→start、親→子）出力は**オプション**として用意し、既定にはしない。
- 1 Span = 1行。列1-8は Span から、source/type は Feature から、列9は属性を保持順に再構成
  （`%XX` 再エンコード、`part` は Span から注入、phase は Span の値）。
- 末尾に `##FASTA` ＋ペプチド配列（あれば）。

### 5.3 意味的ラウンドトリップ（テストの判定基準）
`parse(write(parse(text)))` が `parse(text)` と**モデル等価**であること。等価の定義:
- ディレクティブ集合（kind ＋ value）が一致
- feature 集合が、`id / type / source / spans(seqid,start,end,strand,phase,score,part をソート)
  / attributes(key→値の多重集合) / parent_ids / 親子関係(id ベース)` で一致
- **feature 出現順・属性の字面順・空白・エンコード詳細は一致不問**
- `GffDocument.semantically_equals()` を実装しテストのオラクルにする

### 5.4 エラー処理 / パーサ厳格度（D10）
- 既定は **lenient ＋診断収集**。`document.diagnostics` に `Diagnostic(severity, line_no, code, message)` を蓄積。
- 重大度: ERROR（行解析不能 → スキップ記録）／WARNING（親不在・同 ID type 不一致・非 ASCII・
  start>end 等）／INFO（コメント drop 件数・未知ディレクティブ・FASTA に核酸 等）。
- `strict=True` オプションで最初の ERROR（任意で WARNING）で `GffParseError` を送出。
- 代表対応:
  - 列数不正 / 座標非整数 → ERROR、行スキップ
  - start>end → WARNING（環状起点跨ぎの可能性、トポロジ検証は Phase1 ではしない）
  - 同 ID 異 type → マージせず別扱い ＋ ERROR
  - 親不在（dangling Parent）→ WARNING（`parent_ids` 保持）
  - 未知ディレクティブ → `kind='unknown'` で保持、INFO

---

## 6. テスト戦略（TDD、pytest）

仕様と example から落ちるテストを先に書き、実装で通す。

1. **`test_attributes.py`** — `%XX` エンコード/デコード往復、未エスケープ `,` での多値分割、順序保持、
   `Note` の `%3B` や literal `~`
2. **`test_model.py`** — 複数 span、codon_start 派生、`to_biopython_location()`
   （+鎖昇順 / −鎖降順 / part 順 / trans-splicing remote seqid）、便利プロパティ、`semantically_equals()` の真偽
3. **`test_parser.py`** — 各ディレクティブ、同 ID 集約、ID 無し単独、親子リンク（前方参照含む）、
   診断（親不在 / 列数不正 / 非 ASCII / 同 ID 異 type）、part 順・strand `?`
4. **`test_writer.py`** — 1 span=1行、part 注入、span 別 phase、属性再エンコード、ディレクティブ/FASTA 出力
5. **`test_roundtrip.py`** — 各フィクスチャで `parse(write(parse(t))).semantically_equals(parse(t))`
6. **`test_integration.py`**（slow マーカー・ファイル無ければ skip）— 実ファイル
   （rice_cp / ecoli[gunzip] / chloroplast / marchantia / arabidopsis）の往復＋不変条件
   （例: rps12 の CDS が 3 span、dnaX の ribosomal slippage、transl_except 存在）と
   arabidopsis 121MB のパース性能スモーク

**キュレートされた小フィクスチャ**（難所を1つずつ凝縮し高速・意図明示）:
- `canonical_gene.gff3` — INSDC 仕様の FABP5 例（gene/mRNA/exon/UTR/CDS multi-span）
- `discontinuous_cds.gff3` — multi-exon CDS が1 ID を共有（ycf3 3 セグメント）
- `trans_splicing.gff3` — rps12（exception・part・number・is_ordered・strand `?`・discontinuous intron）
- `transl_except.gff3` — ecoli Sec（transl_except 属性形）＋ ribosomal slippage（dnaX）
- `circular.gff3` — 起点跨ぎ feature（end>length）＋ Is_circular
- `attributes_escaping.gff3` — `%XX`・多値 Dbxref・長い Note
- `forward_parent_ref.gff3` / `dangling_parent.gff3` / `malformed_lines.gff3` — 診断
- `peptide_fasta.gff3` — `##FASTA` ペプチドブロック

**判定基準**: 数値カバレッジ目標ではなく「文書化した各エッジケースにテストがある」。

---

## 7. スコープ境界

### 7.1 内（フェーズ1）
- GFF3 ⇄ オブジェクトモデルの読み書きと意味的ラウンドトリップ
- ディレクティブ・属性・集約・親子グラフ・span/phase/part・診断
- `to_biopython_location()` 橋渡し（後フェーズの抽出・翻訳のため）
- 特殊ケース（trans-splicing・transl_except 属性形・環状座標・order）を**データとして忠実に保持**

### 7.2 外（後フェーズ）
- INSDC プロファイル検証 ＝ Phase3
- NCBI 流 ⇄ INSDC canonical の**正規化変換**
  （transl_except 属性 ⇄ recoded_codon 子 feature、exception+part ⇄ location=join 等）＝ Phase3
- 配列切り出し・翻訳・タンパク検証 ＝ Phase2 以降
- DDBJ MSS 出力 ＝ Phase2
- DDBJ flat file → GFF3 ＝ Phase4
- 定数メモリ・ストリーミング ＝ 将来
- AGAT 連携 ＝ Phase2（外部・任意）
- 環状トポロジ検証 ＝ Phase3

### 7.3 重要な原則
Phase1 は特殊ケースを**忠実に「表現」するが「解釈・変換」はしない**。
例: 環状の end>配列長 は `Span(end>len)` としてそのまま往復させ、座標の巻き戻しや配列抽出は後フェーズ。
trans-splicing の remote seqid は Span に保持し、連結産物の組立ては後フェーズ。

---

## 8. 検証データ（`examples/`）

| データ | 役割 | 含まれる難所 |
|---|---|---|
| marchantia/MpTak1_v7.1.marpolbase.gff | 実際の登録元入力（marpolbase） | gene-mRNA-exon/CDS/UTR、miRNA、複数 transcript 由来の重複 CDS |
| marchantia/chloroplast.gff3 | DDBJ 作成の INSDC 流 GFF（AP025455.1） | trans-splicing＋discontinuous intron、order(is_ordered)、strand `?`、number/exon_number、sequence_feature |
| rice_cp/rice_cp.gff3 | NCBI RefSeq 葉緑体・環状 | trans-splicing(rps12, exception+part)、multi-segment CDS、Is_circular |
| ecoli/*.gff.gz | NCBI RefSeq・環状 | transl_except(Sec)、ribosomal slippage、pseudogene |
| arabidopsis/AT_chr1.gff3 | 巨大（121MB） | 選択的スプライシング、性能 |
| yeast/*.gff.gz | NCBI RefSeq R64 | 一般ケース |

注: example の NCBI/DDBJ 製 GFF は INSDC canonical 形式とは書き方が異なる（例: trans-splicing は
`exception=trans-splicing`＋`part=N` であり `location=join(...)` ではない）。フェーズ1はこれらを忠実に
保持し、canonical 形式への正規化は Phase3 の責務とする。

---

## 9. 既存資産（参考）

`/Users/tanizawa/projects/marchantia/assembly_hifi/annotation/ddbj_submission_hifi/ddbj_gff`
に GFF → DDBJ MSS 変換の実験スクリプト（`gff2mss_for_MP*.py`）がある。bcbio-gff 依存で、
複数 transcript 由来の重複 CDS の扱いに 4 バリエーション（minimum / nonredundant /
redundant_as_misc / full）を作って試行錯誤した経緯がある。フェーズ2 のロジック参考にするが、
土台は本書のとおり自作する。

---

## 10. 未決事項 / 将来検討
- INSDC SO-INSDC feature mapping（列3で許容される SO term サブセット）の取り込みは Phase3。
- `transl_table` のファイルレベル（`#!transl_table`）と CDS レベルの優先関係の解釈は Phase2/3。
- 正規ソート出力（INSDC export 用）の詳細仕様は Phase3 で確定。
- 定数メモリ・ストリーミングパーサ（ヒトゲノム規模）は将来。
