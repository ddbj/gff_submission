# Heterosigma akashiwo ゲノム DDBJ 登録パイプライン 設計書

- 日付: 2026-07-02
- 対象: `dev/heterosigma/` のヘテロシグマ (Heterosigma akashiwo) ゲノムアノテーションを DDBJ MSS 形式へ変換・登録する
- スコープ: **heterosigma 登録ファイル (.ann/.fasta) の生成完了まで**。gff_submission → ddbj_mss_tools の本格移植は別 spec の後続プロジェクトとする(本設計は移植を見据えたモジュール構成にする)

---

## 1. 背景・目的

`dev/heterosigma/braker_with_ncRNA_mtcp.gff3`(以下 mtcp GFF)を入力に、DDBJ MSS 登録ファイルを生成する。mtcp GFF には**構造の異なる 3 系統**のアノテーションが混在している:

| 区分 | seqid | 由来 | GFF 構造 | col-9 product | transl_table |
|---|---|---|---|---|---|
| 核 | scaffold_1..43 | BRAKER4 → AGAT 標準化済 | gene → **transcript** → exon/CDS/UTR/intron/start・stop_codon | なし | 記載なし |
| 葉緑体 | CP | LiftOff | gene → **CDS**(mRNA 層なし) | あり | 11 |
| ミトコンドリア | MT | LiftOff | gene → CDS(タンパク) と gene → tRNA → exon | あり | 記載なし |

核にはさらに、**親 gene を持たない独立 feature** として ncRNA 918(Infernal)・tRNA 712(tRNAscan-SE)・rRNA 37(pybarrnap)、および pseudogene 1356(主に tRNAscan-SE の tRNA 擬遺伝子・子なし)が含まれる。

### 既存 `mss` 変換器との不整合(要対応の核心)

現状の `src/ddbj_gff/mss/convert.py` は **gene → mRNA → CDS** を前提とするため、このままでは:

- 核の `transcript`(SO 上 `misc_RNA` 相当)が `mRNA` にマップされず、**核の全遺伝子が脱落**
- organelle の gene → CDS(mRNA 層なし)も**脱落**
- 親 gene のない ncRNA/tRNA/rRNA 約 1,700 feature も**脱落**
- `_product()` に仕様外の `"protein {gene名}"` フォールバックがあり、id→product 表の入口もない

→ 本質的な課題は「AGAT で綺麗にする」ことではなく、**構造を gene → mRNA →(exon/CDS)に統一**し、**mss 変換器の product / feature 出力ルールを仕様に合わせる**ことである。

---

## 2. 入力データ(`dev/heterosigma/`)

- `Haka_JPv1.fa.gz` — ゲノム FASTA(gzip, 約 347 MB)。scaffold_* と MT・CP を含む
- `braker_with_ncRNA_mtcp.gff3` — mtcp GFF(約 53 MB)
- `annotation_kaas.tsv` — KAAS アノテーション(**CRLF 改行**)。ヘッダ `protein_acc<TAB>Description_2`。キーは `anno1.` を除くと GFF の transcript ID と 13,575 件 100% 一致(例: `anno1.g1.t1`→`g1.t1`、`anno1.both_agree_g273.t2`→`both_agree_g273.t2`)。記載あり 4,446 件 / 空 9,129 件

---

## 3. 全体アーキテクチャ

登録は **2 run** に分ける(核 WGS ↔ organelle):

| run | 対象 seqid | submission_category | mol_type / topology | source 特記 | transl_table 既定 | product 主ソース |
|---|---|---|---|---|---|---|
| **nuclear** | scaffold_1..43 | `WGS` | genomic DNA / linear(+assembly_gap) | submitter_seqid(WGS) | 1 | annotation_kaas.tsv → 無ければ hypothetical protein |
| **organelle** | MT, CP | `GNM` | genomic DNA / **circular** | organelle=(mitochondrion / plastid:chloroplast) | 1(CP は CDS 属性 11 が優先) | col-9 product |

organelle を 1 run にまとめられるのは、egapx2mss 方式の source ビルダが seq ごとに `/organelle`・topology・`ff_definition` を出し分けられるため。transl_table は CDS 属性優先で CP=11(属性)/MT=1(既定)を自動処理する。

### パイプライン(データフロー)

```
braker_with_ncRNA_mtcp.gff3 + Haka_JPv1.fa.gz + annotation_kaas.tsv
        │
 [Step1] AGAT 全体標準化 ──▶ 検証ゲート
        │    NG時 → organelle限定AGAT → コード処理 にフォールバック
        ▼   standardized.agat.gff3                      (dev/heterosigma)
 [Step2] heterosigma 固有 (dev/heterosigma/scripts)
        │    split_by_compartment.py … seqid で nuclear/organelle に分割・pseudogene 除外
        │    build_product_map.py    … kaas → product_map.tsv
        ▼   nuclear.gff3, organelle.gff3, product_map.tsv
 [Step3] normalize (共有, 小改修)  … coding transcript→mRNA 矯正 + directive 付与
        ▼   nuclear.normalized.gff3, organelle.normalized.gff3 (+ validate レポート)
 [Step4] feature = 本プロジェクト mss / COMMON+source = ddbj_mss_tools.common
        │    make_ann.py (adapter, dev/heterosigma/scripts)
        ▼   submission/nuclear.ann/.fasta, submission/organelle.ann/.fasta
```

---

## 4. ステップ詳細

### Step 1 — AGAT 標準化 + 検証ゲート(`dev/heterosigma/`)

目的: gene → mRNA → exon/CDS への構造統一、organelle への mRNA・exon 補完、**親なし RNA への gene 親付与**。

- 実行(amd64 コンテナ): `agat_convert_sp_gxf2gxf.pl -g braker_with_ncRNA_mtcp.gff3 -o standardized.agat.gff3`
- **検証ゲート**(`scripts/verify_agat.py`)— 以下を検査し、NG なら方針をフォールバック:
  1. 階層: すべての mRNA/exon/CDS/RNA が gene まで解決する Parent 鎖を持つ(孤児ゼロ)
  2. organelle: CP/MT の protein 遺伝子に mRNA・exon が補完されている
  3. 親なし RNA: 核 ncRNA/tRNA/rRNA に gene 親が付与されている
  4. ID 健全性: 空 ID / 重複 ID / marchantia 空 Note 事案(645 孤児)のパターンが出ていない
  5. 件数: gene/mRNA/CDS 数が入力から想定外に増減していない
- **フォールバック順**: 戦略3(全体 AGAT)→ 戦略2(organelle 限定 AGAT + 核はそのまま)→ 戦略1(AGAT 不使用・normalize/Step2 コードで構造統一)
- 注意点(ゲートで重点確認): 括弧付き tRNA 名 `trnA(tgc)`、分割 rpl2 遺伝子、`transcript`→`mRNA` 改名の有無、属性(`gene_biotype`/`product`/`Dbxref`/`transl_table`)の保持

### Step 2 — heterosigma 固有スクリプト(`dev/heterosigma/scripts/`)

- **`build_product_map.py`**: `annotation_kaas.tsv` を CRLF 除去 → `anno1.` 除去 → 末尾 ` [EC:...]` 除去(EC 除去のみ。`/EC_number` 化や `---` 正規化はしない)→ `product_map.tsv` を出力。
  - キーは **transcript_id と gene_id の両方**を同一 product で出力(`g1.t1` と `g1` の両方)。これにより mss の代表 transcript 選択に依存せず引ける。
  - 空 description の行は出力しない(= hypothetical protein に落ちる)。
  - EC 除去後に空文字になる場合も出力しない。
- **`split_by_compartment.py`**: seqid で `nuclear.gff3`(scaffold_*)と `organelle.gff3`(MT, CP)に分割。
  - **pseudogene 除外**もここで実施: `type == "pseudogene"` または属性 `gene_biotype == "pseudogene"` の feature とその子孫を出力しない(除外件数をレポート)。

### Step 3 — normalize(共有モジュール, 小改修)

各分割 GFF に対し実行。

- **新規パス**: `pass_coerce_transcript_to_mrna` — CDS を子に持つ `transcript` を `mRNA` に矯正(config gate)。Step1 AGAT が既に改名済みなら no-op。
- 既存の directive 付与(`##species` taxid=2829 / `#!transl_table` / `##sequence-region` は FASTA 長から。FASTA は gzip 対応)。
- 出力後に `validate` を実行しレポートを保存(登録前レビュー用)。

### Step 4 — feature(本プロジェクト mss)+ COMMON/source(ddbj_mss_tools.common)

**feature 生成**は本プロジェクト `mss` が担当、**COMMON/source/gap** は `ddbj_mss_tools` の共有 `common` パッケージを**実コードとして再利用**する。両者を `dev/heterosigma/scripts/make_ann.py`(adapter)が結合する。

adapter の流れ(egapx2mss の `write_ddbj_ann` を踏襲):
1. 本プロジェクト mss の **feature-only API** で、entry(seqid)ごとの feature 行(source を除く)を得る
2. `common.models.load_common_json(common_*.json)` で `CommonModel` を得る
3. `common.source_builder.load_sequence_roles(sequence_roles.tsv)` で per-entry の role を得る
4. `common.common_builder.create_common(...)` で COMMON 行を生成
5. entry ごとに: `is_circular` なら `TOPOLOGY circular` 行 → `source_qualifier()` + `ff_definition()` による source 行 → feature 行 → `common.gap_annotator.annotate_gaps()` によるアセンブリギャップ行、を連結
6. `.ann` と cleaned `.fasta`(gzip 入力を解凍)を出力

`is_wgs` 判定は `common.source_builder` と同じ「全 entry が unplaced か」に従う(nuclear=全 unplaced→WGS、organelle=organelle 型→非 WGS で per-entry source + TOPOLOGY)。

---

## 5. 共有モジュール改修(`src/ddbj_gff/`)

移植を見据え、**feature 生成**と**ドキュメント組み立て(source/COMMON/gap/emit)**を分離する。

### `mss/convert.py`

1. **product ルール是正**(`_product`)
   - `(1)` 外部 `product_map` を `mRNA.id` → 無ければ `gene.id` で参照(非空なら採用)→ `(2)` col-9 の `product` → `(3)` `product_default`("hypothetical protein")
   - **`"protein {gene名}"` フォールバックを削除**
2. **feature-only API**: `build_entry_features(doc, seqs, cfg, product_map, diagnostics) -> dict[seqid, list[MssFeature]]`(source を含まない)を新設。既存 `convert()` はこれを呼び、standalone 用途では source を前置し gap を後置する(後方互換)。adapter は `build_entry_features` を直接使う。
3. **親 gene のない ncRNA/tRNA/rRNA のトップレベル出力**: gene ルート処理後、`doc.roots` の type ∈ {ncRNA, tRNA, rRNA(+SO 変種)} を単独 feature として emit(locus_tag 採番、product/`ncRNA_class`/`Dbxref`→db_xref/note を付与)。AGAT で gene 親が付いた分は従来の gene→RNA 経路で処理され、**両輪**で取りこぼしを防ぐ。
4. **pseudogene 除外**: 念のため変換器側でも `type=="pseudogene"` / `gene_biotype=="pseudogene"` をスキップ(Step2 で除去済みでも二重防御)。
5. **内部 stop → `misc_feature`**: CDS を翻訳し body に内部 `*` があれば、CDS/mRNA を出さずにその位置へ `misc_feature` を 1 つ出力(translate/product なし、`note` に概要: 遺伝子名・内部 stop 検出の旨)。末端の開始/終止 codon 欠如のみの partial は従来どおり `<`,`>` 付き CDS。
6. **transl_table**: CDS 属性優先(既存)。CP=11(属性)/MT・核=config 既定 1。

### `mss/config.py`

- `[product] map = "path"` を追加し product_map を読み込む(feature-side 設定)。
- 本 heterosigma 経路では `[source]` / `[assembly_gap]` / `[chromosome]` / `[ff_definition]` は使わない(COMMON/source/gap は `common` パッケージが担当)。使用するのは `[locus_tag]`・`[cds]`・`[product]` セクション。

### `mss/emit.py`

- feature 1 件を 5 列行に変換する `feature_rows(entry_or_blank, feat)` を切り出し、adapter と共有できるようにする。

### `normalize/`

- `pass_coerce_transcript_to_mrna` 追加(§Step3)。gzip FASTA 対応。

### FASTA gzip 対応(`mss/cli.py`, `normalize/cli.py`, adapter)

- パスが `.gz` の場合は `gzip.open` で読む。

---

## 6. COMMON / source 詳細(egapx2mss 踏襲)

`ddbj_mss_tools` の共有 `common` パッケージを再利用する(移植の第一歩)。

- **入力**:
  - `common_nuclear.json` / `common_organelle.json`(`CommonModel`): `DATATYPE` / `KEYWORD` / `DBLINK`(project/biosample/SRA)/ `SUBMITTER` / `REFERENCE` / `DATE`(hold_date)/ `SOURCE`(organism="Heterosigma akashiwo", strain=…, mol_type="genomic DNA")/ `SOURCE_IDENTIFIER`(例 "strain")/ `ASSEMBLY_GAP`(nuclear のみ)
  - `sequence_roles.tsv`(seq_id / type / seq_name / status / topology):
    ```
    # nuclear は記載しない(未記載 = unplaced = WGS 判定)
    MT	organelle	mitochondrion	complete	circular
    CP	organelle	plastid:chloroplast	complete	circular
    ```
    `seq_name` が `/organelle` の値になり、`ff_definition` では `_organelle_code` で `plastid:chloroplast`→`chloroplast`、`mitochondrion`→`mitochondrial` に変換される。
- **category**: nuclear=`WGS`(DATATYPE=WGS, GNM 継承)、organelle=`GNM`(完全ゲノム)。
- **依存の閉じ込め**: adapter が `sys.path` に `ddbj_mss_tools/src` を追加して `common` を import する。依存は `dev/heterosigma` に閉じ、gff_submission 公開レポジトリは汚さない。

---

## 7. 設定ファイル(`dev/heterosigma/`)

- `nuclear.mss.toml` / `organelle.mss.toml`(feature-side): `[locus_tag]` prefix=仮値(例 `HAKA`)・width/start/step、`[cds] transl_table`(1 / 1)、`[product] default="hypothetical protein"`・nuclear は `map="product_map.tsv"`
- `common_nuclear.json` / `common_organelle.json`(§6)
- `sequence_roles.tsv`(§6)
- locus_tag prefix は仮値。BioSample 登録後の正式 prefix に config で差し替える。

---

## 8. エラー処理・検証

- Step1: 検証ゲート(§Step1)。NG はフォールバック。
- Step3: `validate` レポートを保存し登録前にレビュー。
- Step4: mss diagnostics(no-cds / internal-stop→misc_feature / translation 系)を集計。`common` の category バリデーション(必須 source/DBLINK/ST_COMMENT 欠落警告)も活かす。
- 共通: FASTA seqid(scaffold_*/MT/CP)と GFF seqid の一致を検証。

---

## 9. テスト方針(TDD)

共有コードは `tests/` に小 fixture で先にテスト:

- product 3 ルール(map hit / col-9 / hypothetical、`"protein {gene}"` 廃止の確認)
- 親なし ncRNA/tRNA/rRNA のトップレベル出力
- pseudogene 除外
- 内部 stop → misc_feature(translate/product なし・note あり)
- 末端 partial の `<`,`>` 表記(既存挙動維持)
- transl_table 属性優先(CP=11)
- gzip FASTA 読込
- feature-only API(source を含まないこと)

heterosigma スクリプトは単体テスト:

- `build_product_map`: CRLF / `anno1.` 除去 / 末尾 EC 除去 / 空行スキップ / transcript+gene 両キー
- `split_by_compartment`: seqid 分割 / pseudogene 除外

統合(slow マーカー): 数 scaffold + MT + CP の小サブセットで end-to-end(Step1→4)。最後に amd64 コンテナで全ゲノム実走。

adapter は `ddbj_mss_tools.common` を import するため、テストは `ddbj_mss_tools/src` を `sys.path` に追加して実行(依存が無い CI では skip マーカー)。

---

## 10. リスク・未決事項

- **AGAT の挙動**: 括弧付き tRNA 名・分割 rpl2・`transcript`→`mRNA` 改名の有無・属性保持 → 検証ゲートで判定、NG はフォールバック。
- **性能・メモリ**: 全体パース + 約 1 GB 相当 FASTA。amd64 コンテナで実走。
- **COMMON 実値**: BioProject/BioSample/submitter/reference/hold_date は登録直前にユーザが実値提供(それまでプレースホルダ)。
- **locus_tag prefix**: 仮値。正式 prefix に差し替え。
- **MT の遺伝コード**: table 1 を採用(要 CDS 翻訳での妥当性確認: 内部 stop が出れば §内部stop→misc_feature へ)。
- **organelle の BioSample/hold**: MT/CP を 1 run(共通 BioProject/BioSample 前提)。別々が必要なら run を分割(3 run へ)。

---

## 11. 成果物

- `dev/heterosigma/standardized.agat.gff3`(Step1)
- `dev/heterosigma/scripts/{verify_agat.py, build_product_map.py, split_by_compartment.py, make_ann.py}`
- `dev/heterosigma/{nuclear.gff3, organelle.gff3, product_map.tsv}`(Step2)
- `dev/heterosigma/{nuclear.normalized.gff3, organelle.normalized.gff3}` + validate レポート(Step3)
- `dev/heterosigma/{nuclear.mss.toml, organelle.mss.toml, common_nuclear.json, common_organelle.json, sequence_roles.tsv}`
- `dev/heterosigma/submission/{nuclear.ann, nuclear.fasta, organelle.ann, organelle.fasta}`(Step4)
- `src/ddbj_gff/` 改修(mss: product ルール / feature-only API / 親なし RNA / pseudogene 除外 / misc_feature / gzip、normalize: transcript→mRNA)+ テスト

## 12. 将来(別 spec)

gff_submission の feature 生成(normalize/mss/validate)を `ddbj_mss_tools` のサブツールとして移植し、COMMON/source/gap/fasta/submission_category を共有 `common` に一本化する。本設計の feature/document 分離と `common` 再利用は、その移植の前提を満たすように作る。
