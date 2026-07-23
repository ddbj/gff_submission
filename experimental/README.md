# Data conversion scripts from GFF3 to MSS format
These script files were used to submit the genome sequences for Marchantia polymorpha.  
The script file is specific to this data, and cannot be applied to other data without modification.

## Prerequisites
- The script files are tested with Python 3.11
- Requires bcbiogff (tested with v.0.6.6)
```
conda install -c bioconda bcbiogff
```

Due to recent updates of bcbiogff and Biopython, the script files may not work with the latest versions.

## Download sample data
```
curl -LO https://hifi.marchantia.info/data/experimental/MpTak1_v7.1.ddbj.fa
curl -LO https://hifi.marchantia.info/data/experimental/MpTak1_v7.1.ddbj.gff
```

## Run
```
python gff2mss_for_MP.py

or either of
python gff2mss_for_MP_minimum.py
python gff2mss_for_MP_nonredundant.py
python gff2mss_for_MP_redundant_as_misc.py
```

## 使用上の注意点
### GFFの構成
gene-ｍRNA-CDS,exon　という階層構造になっていて、ID-parentで親子関係が記載されていること。

### gene, product, locus_tag の各qualifierの扱い
本サンプルデータでは、CDSおよびmRNAを記載する際に必要なこれらの情報は GFFファイルのmRNA行にあらかじめ埋め込まれていることを前提とする。
同じ遺伝子座かRa複数のmRNA, CDSが生じる場合、それらに対して共通する gene, locus_tag を記載する必要がある。

### mRNA, CDS, UTR
mRNA は GFF中の exon 行の位置情報を連結し join を使って位置を表す。  
CDS は GFF中の CDS 行の位置情報を連結し join を使って位置を表す。  
5'UTR および 3'UTR は GFF中の three_prime_UTR または five_prime_UTR 行の位置情報を連結し join を使って位置を表す。  

- partial CDS, mRNAの扱い  
UTRに該当する領域が存在しない場合、ｍRNA は partial の扱いとし、 `<` または `>` を使った partial location で位置を表す。  
start codon および stop codon が存在しない場合、CDS の位置は同様に partial locationで表す。　　


### script fileについて
- gff2mss_for_MP.py
    GFF3に含まれる全ての情報 (mRNA, CDS, UTR, miRNA) をMSS形式に変換。  
    同じ遺伝子座から複数のtranscriptが生じる場合に、重複するCDSやUTRが多すぎたため、使用は断念 (UTR だけが異なり CDS が共通する場合や、途中のexonの有無の違いでUTRが共通する場合が該当)  

- gff2mss_for_MP_nonredundant.py
    mRNA, CDS, miRNA を変換  
    同一の CDS は１つのみ記載する (note にどのtranscriptに由来するかを記載)  
    mRNAとCDS が 1:1 にならないのが欠点  

- gff2mss_for_MP_redundant_as_misc.py
    mRNA, CDS, miRNA を変換  
    同一の CDS が存在する場合最初の１つのみ記載し、それ以外は misc_feature として記載 (misc_featureにはどのCDSと同一かを記載)  
    misc_featureも含めれば、mRNAとCDS との関係は 1:1 になる  

- gff2mss_for_MP_minimum.py
    mRNA, CDS, miRNA を変換  
    mRNAとCDSをそれぞれの遺伝子座について１件のみ記載した最小構成。