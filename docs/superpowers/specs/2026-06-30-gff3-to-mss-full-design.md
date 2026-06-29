# 設計書: GFF3 → DDBJ MSS 変換 — Phase 2 残り機能（複数 transcript / 非コードRNA）

- 日付: 2026-06-30
- 対象: フェーズ2 minimal MVP（`ddbj_gff.mss`、merged）を拡張し、複数 transcript（重複CDSの扱いをモード切替）と非コードRNA（miRNA/pre_miRNA ほか）に対応する。
- 前提: `docs/superpowers/specs/2026-06-29-gff3-to-mss-design.md`（minimal MVP）と、その実装（`ddbj_gff.mss` 全モジュール）が `main` に存在する。

---

## 1. 背景とゴール

minimal MVP は「遺伝子座あたり先頭 transcript の mRNA+CDS」のみを出力した（spec §9.2 で複数 transcript・UTR・各種RNA を後回し）。本拡張で marchantia 核ゲノム登録に必要な残りを実装する:
1. 複数 transcript（重複CDSの扱い）
2. 非コードRNA（miRNA / pre_miRNA、汎用 tRNA/rRNA/ncRNA マッピング）

UTR feature は出力しない（mRNA の join 位置で UTR 領域は暗黙表現される。決定済み）。
オルガネラ特殊ケース（trans_splicing / transl_except / 環状）は引き続き **Phase 3**（INSDC 正規化・検証）の範囲。

marchantia GFF の規模: 18,200 gene / 22,172 mRNA（約4,000が複数 transcript、最大9）/ 28,537 five_prime_UTR / 24,791 three_prime_UTR / 315 miRNA / 256 pre_miRNA。tRNA/rRNA/ncRNA は核ゲノムには無い。

---

## 2. 確定した方針（意思決定ログ）

| # | 論点 | 決定 |
|---|---|---|
| F-D1 | 複数 transcript の重複CDS | **モード切替** `minimal` / `nonredundant` / `full` を選択可能に |
| F-D2 | モード指定 | `config.toml [transcript] mode` ＋ CLI `--mode`（CLI 優先）。不正値は `GffParseError` |
| F-D3 | 既定モード | **`nonredundant`** |
| F-D4 | UTR feature | 出力しない（mRNA 位置で暗黙表現） |
| F-D5 | 非コードRNA | 出力。型→MSS feature マッピング（miRNA→ncRNA[ncRNA_class]、pre_miRNA→precursor_RNA、tRNA/rRNA/ncRNA…） |
| F-D6 | CDS 重複排除 | `nonredundant` で**遺伝子座内・location 文字列キー**で排除。note に共有 transcript 群を列挙 |
| F-D7 | アーキテクチャ | 案A: `convert` の per-gene 処理を `build_gene_features(gene, mode, …)` に置換。既存ビルダ再利用 |
| F-D8 | 範囲外 | オルガネラ特殊ケース・INSDC 正規化/検証・AGAT・DDBJバリデータ（= Phase 3 以降） |

---

## 3. モード設定（config / CLI）

- `MssConfig` に `transcript_mode: str = "nonredundant"` を追加。
- `config.toml`:
  ```toml
  [transcript]
  mode = "nonredundant"   # minimal | nonredundant | full
  ```
- `load_config` は `[transcript].mode` を読み、`{"minimal","nonredundant","full"}` 以外なら `GffParseError`（code `invalid-mode`）。未指定なら既定 `nonredundant`。
- CLI に `--mode {minimal,nonredundant,full}` を追加。指定時は config 値を上書き（`main` が `cfg.transcript_mode` を差し替えてから `convert` を呼ぶ）。
- `convert(doc, seqs, cfg, common_rows, *, strict=False)` は `cfg.transcript_mode` を参照（シグネチャは不変）。

出力方針（per 遺伝子座）:

| モード | mRNA | CDS |
|---|---|---|
| `minimal` | 先頭 transcript のみ | 先頭 transcript のCDSのみ |
| `nonredundant` | 全 transcript | location 同一CDSは1回（note に由来 transcript 群） |
| `full` | 全 transcript | 全 transcript のCDS（重複も出力） |

---

## 4. per-gene 構築ロジック（`convert.py`）

現在の per-gene ブロック（`_representative_mrna` → mRNA+CDS）を新関数に置換:

```
build_gene_features(gene, mode, locus_tag, genome_seq, cfg, diagnostics) -> list[MssFeature]
```

locus_tag は1遺伝子座につき1つ（ゲノム全体連番アサイナ、既存どおり。`convert` 側で `assigner.assign(gene)` を呼び渡す）。

### 4.1 タンパク質コード遺伝子（mRNA 子あり）
- transcripts = type=="mRNA" の子。id 順（`.1`,`.2`,…）で安定ソート。
- `minimal`: 代表 transcript（`.1` 優先・無ければ先頭）のみ。複数あれば `multi-transcript` WARNING。→ mRNA + CDS。
- `full`: 各 transcript → `build_mrna_feature` ＋ `build_cds_feature`（重複CDSもそのまま append）。
- `nonredundant`:
  - 各 transcript → `build_mrna_feature` を append。
  - 各 transcript の `build_cds_feature` を構築し、**CDS の location 文字列をキー**に遺伝子座内 dict で集約。初出時は (cds, [tid])、既出時は tid を追記。
  - 全 mRNA を出した後、ユニークCDSを初出順に append。共有が複数 transcript の CDS は、その `note` 値を `submitter_gene_id: G, submitter_transcript_id: t1, t2, …` に書き換え（ヘルパで note qualifier を置換）。
- 出力順（遺伝子座内）: `full`/`minimal` は mRNA→CDS の対を transcript 順。`nonredundant` は 全mRNA→ユニークCDS群。

### 4.2 非コード遺伝子（mRNA 子なし）
`build_noncoding_features(gene, locus_tag, genome_seq, cfg, diagnostics) -> list[MssFeature]`:
- 型→MSS feature マッピング:
  `pre_miRNA→precursor_RNA` / `miRNA→ncRNA`(+`ncRNA_class=miRNA`) / `tRNA→tRNA` / `rRNA→rRNA` / `ncRNA,snRNA,snoRNA→ncRNA`(+`ncRNA_class=<型>`) / `tmRNA→tmRNA` / その他認識外→`misc_RNA`。
- 認識可能な RNA 子が無ければ `no-rna` WARNING で空を返す。
- 各 RNA 子について:
  - location = `collect_spans(rna, "exon")` があればそれ、無ければ `rna.spans`（marchantia の miRNA/pre_miRNA は単一 span）。`build_insdc_location` を使用。**partial 判定はしない**（非コードに開始/終止コドンの概念なし）。
  - qualifier: `locus_tag` / ncRNA 系なら `ncRNA_class` / `gene`（あれば）/ GFF `Note`→`note` / submitter note（`submitter_gene_id`, `submitter_transcript_id: <rna.id>`）。
- marchantia の miRNA 遺伝子（gene＋pre_miRNA＋miRNA×N）→ `precursor_RNA` 1個＋`ncRNA` N個、同一 locus_tag。モードに依存しない（CDS が無く重複問題が起きない）。

### 4.3 既存ビルダの再利用
`build_mrna_feature` / `build_cds_feature` / `build_insdc_location` / `collect_spans` / `_submitter_note` はそのまま使用。CDS dedup の note 書き換えのみ新ヘルパ。

---

## 5. 診断 / エラー処理
- `multi-transcript` WARNING は **`minimal` モードのみ**（他モードは全 transcript を保持し破棄しないため）。
- 非コード遺伝子に認識可能 RNA 子が無い → `no-rna` WARNING。
- `[transcript].mode` または `--mode` が不正値 → `GffParseError`（`invalid-mode`）。
- 既存の診断（`missing-sequence` ERROR、`no-cds`、`source-missing-qualifier`、翻訳検証 等）は不変。

---

## 6. テスト戦略（TDD・pytest）
1. `test_mss_config`（拡張）: `[transcript].mode` 読込（既定 `nonredundant`、不正値で `GffParseError`）。
2. `test_mss_convert`（拡張）:
   - 2 transcript・CDS同一/UTR差の遺伝子座: `minimal`→mRNA1/CDS1、`full`→mRNA2/CDS2、`nonredundant`→mRNA2/CDS1（CDS note に両 transcript id）。
   - CDS が異なる 2 transcript → `nonredundant` でも CDS2（排除されない）。
   - miRNA 遺伝子（gene＋pre_miRNA＋miRNA×2）→ `precursor_RNA`＋`ncRNA`×2（`ncRNA_class=miRNA`、同一 locus_tag）。
   - `no-rna`（RNA 子なし非コード遺伝子）WARNING。
3. `test_mss_cli`（拡張）: `--mode full`/`--mode nonredundant` が反映、不正 `--mode` でエラー終了。
4. `test_mss_snapshot`（拡張）: mini フィクスチャに「複数 transcript 遺伝子」と「miRNA 遺伝子」を追加し、既定 `nonredundant` の期待 .ann を再生成・固定。単一 transcript 部分は全モード同一なので既存ゴールド該当行は不変。必要なら `full` の期待 .ann も別途固定。
5. `test_mss_integration`（更新）: 実 marchantia を既定 `nonredundant` で変換 → ERROR 0、mRNA は全 transcript 分（>22,000）、ユニークCDS < mRNA 数（重複排除が起きている）、`precursor_RNA`/`ncRNA` feature が存在、各 CDS が必須 qualifier を持つ。

---

## 7. スコープ境界
### 内
transcript 3モード（minimal/nonredundant/full）／CDS 遺伝子座内 location 重複排除＋note集約／非コードRNA（miRNA・pre_miRNA＋汎用 tRNA/rRNA/ncRNA マッピング）／モード config・CLI。

### 外（Phase 3 以降）
オルガネラ特殊ケース（trans_splicing / transl_except / 環状座標）／厳密 INSDC バリデーション／NCBI流⇄INSDC正規化／AGAT 前処理／DDBJ MSS バリデータ連携／pseudogene／partial の start_range/end_range 表記。

---

## 8. 既存資産（参考）
`…/ddbj_submission_hifi/ddbj_gff/experimental/` の 4 スクリプト（gff2mss_for_MP{,_minimum,_nonredundant,_redundant_as_misc}.py）が各モードの実装参考。`full` は「重複が多すぎて断念」、`nonredundant` は「同一CDSは1つ・note に由来記載・mRNA:CDS が 1:1 にならない」と README にある（本実装の `nonredundant` と一致）。`create_other_RNA_features` の型マッピングが §4.2 の参考。
