# 設計書: 正規化器 feature-type 修正（フェーズ3-B / Finding A・B 是正）

- 日付: 2026-06-30
- 対象: `ddbj_gff.normalize.passes.pass_so_terms` の局所修正。SO-term の collapse 対象を**コア型 whitelist** に絞り、INSDC feature-key への書換をやめ、SO→INSDC-feature 変換は Phase 2 生成に一本化する。新規モジュールなし。
- 前提: Phase 3-B（共通ケース MVP）が `main` にある。本修正は 3-B 実走（ゼニゴケ核ゲノム）で露呈した2件の合成バグ（Finding A/B）の是正。

---

## 1. 背景: ゼニゴケ実走で露呈した2つの合成バグ

2026-06-30、実ゼニゴケ核ゲノムで `3B normalize → 3A validate → Phase2 MSS` を通した結果、パイプラインは exit 0 で完走したが2件の問題が判明:

- **Finding A（col3 語彙の不整合）**: `pass_so_terms` が col3 を INSDC feature 名（`feature-mapping.tsv` 列4）へ書換（`five_prime_UTR→5'UTR` 等）。だが 3-A の `feature-type-not-insdc` は SO-term 列（列2）との照合のため、正規化後に **53,328 件の新規 WARNING** が発生（生GFFでは 0）。3-B の「修正」が 3-A の警告を生んでいた。col4 が SO term でもある型（exon/CDS/gene/ncRNA）は無害だが、SO term でない col4（5'UTR/mobile_element/misc_RNA…）で齟齬。
- **Finding B（pre_miRNA の劣化）**: 3-B が `miRNA/pre_miRNA→ncRNA` に書換えたため、Phase 2 が必要とする元 SO type が破壊された。Phase 2 は `_RNA_MAP`（`pre_miRNA→precursor_RNA`, `miRNA→ncRNA`）と「ncRNA_class を type から導出」する実装を持つため、3-B の書換後は type=`ncRNA` を見て `ncRNA_class=other` に落ちた（3-B が付けた属性は無視）。

**根本原因**: Phase 2 は既に SO→INSDC-feature の変換を自前で持つ（`_RNA_MAP`/`_STRUCTURAL`/UTR は mRNA partial 表現）。3-B が col3 を INSDC 名へ書換えるのは責務違反かつ情報破壊。3-A は列2の全 SO term を許可するため、SO term を保てば warning も出ない。

---

## 2. 確定した方針（意思決定ログ）

| # | 論点 | 決定 |
|---|---|---|
| F-D1 | スコープ | feature-type 修正のみを単独サブプロジェクトで先行（特殊ケース canonical化/環状/プレースホルダ qualifier は後続の 3-B-full-2 へ） |
| F-D2 | collapse 基準 | **コア型 whitelist** `_COLLAPSE_TARGETS = {gene, mRNA, CDS, exon, intron}`。書換先(col4)がこの集合に含まれる時のみ collapse。col4 が SO term でも whitelist 外（ncRNA/5'UTR 等）なら放置 |
| F-D3 | SO→INSDC feature 変換 | 3-B では行わない。Phase 2 生成の責務（既存 `_RNA_MAP` 等が担う） |
| F-D4 | qualifier 付与 | 既存ロジック維持（具体値→属性、プレースホルダ→needs-manual、既存キー非上書き）。whitelisted target の qualifier は実データ上すべて具体値なので needs-manual は実質不発火だが防御的に残す |
| F-D5 | 放置型の report | whitelist 外で放置する型は report に出さない（53k ノイズ回避）。未マップ型のみ `unmapped-type` を従来通り report |

依存・不変部分: `vocab.py`（`feature_qualifiers`/dedup 含む）、`pass_directives`、`normalize()`、CLI、`report.py` は変更なし。

---

## 3. 変更内容（`pass_so_terms`）

`_COLLAPSE_TARGETS = {"gene", "mRNA", "CDS", "exon", "intron"}` を追加。各 feature について:

```
target = vocab.insdc_map.get(f.type)
if target is None:                              # 既知 SO term でない
    → Change("unmapped-type", …)（unresolved）, 不変
elif target == f.type:                          # 既にコア名
    → no-op
elif target in _COLLAPSE_TARGETS:               # コア型へ collapse
    f.type = target
    → Change("rename-type", …)
    各 qualifier: 具体値→属性付与+Change("add-qualifier"); プレースホルダ→Change("needs-manual"); 既存キーは skip
else:                                            # target はあるが whitelist 外（ncRNA/5'UTR/mobile_element 等）
    → 放置（書換えない・Change を出さない）       ★新規分岐
```

**型ごとの挙動:**

| 入力 type | col4 | 動作 |
|---|---|---|
| `coding_exon` / `noncoding_exon` / … | exon | →exon |
| `pseudogenic_CDS` | CDS | →CDS ＋ `pseudo` |
| `processed_pseudogene` 等 pseudogene | gene | →gene ＋ `pseudogene=…` |
| `spliceosomal_intron` 等 | intron | →intron |
| `five_prime_UTR` / `three_prime_UTR` | 5'UTR/3'UTR | 放置（Phase2 が partial 表現） |
| `miRNA` | ncRNA | 放置（Phase2: ncRNA[ncRNA_class=miRNA]） |
| `pre_miRNA` | ncRNA | 放置（Phase2: precursor_RNA） |
| `ncRNA_gene` / `mobile_genetic_element` / `binding_site` 等 | ncRNA/mobile_element/misc_binding | 放置（whitelist 外） |
| 未マップ（既知 SO term でない） | — | `unmapped-type` report、不変 |

---

## 4. テスト

**`tests/test_normalize_pass_so_terms.py` 更新:**
- 維持: `pseudogenic_CDS→CDS+/pseudo`、`processed_pseudogene→gene+/pseudogene`、`CDS→CDS` no-op、未マップ→`unmapped-type`。
- 差し替え（WL外で放置に変わる）: `binding_site`・`mobile_genetic_element` のテストを「type 不変・rename/qualifier Change なし」を確認する形に。
- 新規: `miRNA`/`pre_miRNA`/`five_prime_UTR` が**不変**（rename Change なし）、`coding_exon→exon`、`spliceosomal_intron→intron`。

**`tests/test_normalize_integration.py`:**
- 既存 `test_messy_fixture_normalizes_and_validates`（`coding_exon`/`pseudogenic_CDS` 使用 = WL内）は green のまま。
- 受け入れテスト追加: `gene/mRNA/coding_exon/pseudogenic_CDS/miRNA/pre_miRNA/five_prime_UTR` を含む小フィクスチャを normalize → (a) `coding_exon→exon`・`pseudogenic_CDS→CDS` のみ書換、miRNA/pre_miRNA/UTR は不変、(b) `validate(normalized)` の code 集合に `feature-type-not-insdc` が**出ない**。

**受け入れ基準:**
1. 上記ユニット＋統合テスト green、全体スイート回帰なし（特に `test_normalize_vocab.py` の dedup テストは vocab 不変につき不変）。
2. マージ後にゼニゴケ実走を再実行し、`feature-type-not-insdc` 0 件・MSS で `pre_miRNA→precursor_RNA`・`miRNA→ncRNA[ncRNA_class=miRNA]` を確認。

**スコープ境界:**
- 内: `pass_so_terms` の collapse 基準を whitelist に変更＋テスト。
- 外（後続 3-B-full-2）: 特殊ケース canonical化（trans-splicing/transl_except/anticodon）、環状座標、プレースホルダ qualifier の能動的処理、`ncRNA_gene` 等 gene-level RNA の扱い精緻化。
