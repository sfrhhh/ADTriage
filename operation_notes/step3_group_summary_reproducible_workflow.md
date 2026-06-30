# Step 3: ADT Group Summary 可复现操作说明

## 1. 目标

Step 3 的目标是按用户指定的 Seurat metadata 分组列，计算每个 group level 中的 feature 平均信号，用于后续按细胞群或 cluster 进行 ADT/RNA/CD45 isoform pseudobulk 级别复核。

输出目录：

```text
ADTriage_output/03_group_summary/
```

核心输出：

```text
adt_group_summary.tsv
rna_group_summary.tsv
cd45isoform_group_summary.tsv
```

## 2. 输入文件

```text
tests/data/test_seurat_qc.rds
ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv
```

本次指定的分组列：

```text
Adult_cHSPCs_Ultima_majority_voting
Standardpcaclusters
```

输入 Seurat 对象检查结果：

```text
cells = 13289
assays = RNA, HTO, ADT, CD45isoforms
Adult_cHSPCs_Ultima_majority_voting levels = 9
Standardpcaclusters levels = 15
```

## 3. 使用脚本

```text
scripts/build_adt_group_summary.R
```

脚本功能：

```text
1. 读取 Seurat RDS 和 Step 1 llm_review mapping 表。
2. 检查 group_by metadata 列和必要 assay。
3. 在内存中重新计算 ADT/CD45isoforms/HTO 的 CLR data。
4. 对每个 group_by、group_level、feature 计算平均信号。
5. 对缺失 feature 不中断，写入 status/message。
```

## 4. 归一化规则

脚本读取 RDS 后在内存中重新归一化，不覆盖原始 RDS。

```r
bm <- NormalizeData(bm, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = FALSE)
bm <- NormalizeData(bm, assay = "CD45isoforms", normalization.method = "CLR", verbose = FALSE)
bm <- NormalizeData(bm, assay = "HTO", normalization.method = "CLR", verbose = FALSE)
```

RNA 使用对象中既有 `RNA` assay 的 `data` layer，即 log-normalized signal。

## 5. 复现命令

在项目根目录运行：

```bash
cd /data/project/yanhuixu/frsui/sc_project/ADTriage
```

### 5.1 R 脚本语法检查

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript -e 'parse(file="scripts/build_adt_group_summary.R"); cat("R parse OK\n")'
```

### 5.2 生成 group summary

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript scripts/build_adt_group_summary.R \
  --seurat-rds tests/data/test_seurat_qc.rds \
  --mapping ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv \
  --outdir ADTriage_output/03_group_summary \
  --group-by Adult_cHSPCs_Ultima_majority_voting,Standardpcaclusters
```

## 6. 输出文件说明

### 6.1 adt_group_summary.tsv

每个 ADT/HTO feature 在每个 group level 的平均 CLR 信号。

关键字段：

```text
group_by
group_level
n_cells
feature_id
feature_name
target_class
assay
feature
mean_CLR
mean_ADT_CLR
mean_HTO_CLR
status
message
```

说明：

```text
target_class = hashtag 时 assay = HTO。
其他 feature 默认 assay = ADT。
```

### 6.2 rna_group_summary.tsv

普通 gene-level ADT 的 paired RNA 验证信号。

关键字段：

```text
group_by
group_level
n_cells
feature_id
feature_name
target_class
validation_assay
validation_feature
mean_ADT_CLR
mean_RNA_lognormalized
status
message
```

### 6.3 cd45isoform_group_summary.tsv

CD45 isoform ADT 的 paired CD45isoforms 验证信号。

关键字段：

```text
group_by
group_level
n_cells
feature_id
feature_name
target_class
validation_assay
validation_feature
mean_ADT_CLR
mean_CD45isoform_signal
status
message
```

## 7. 本次运行统计

```text
group_by columns = 2
total group levels = 24
adt_group_summary.tsv rows = 4200
rna_group_summary.tsv rows = 4032
cd45isoform_group_summary.tsv rows = 24
```

状态统计：

```text
ADT summary:
ok = 4152
missing_feature = 48

RNA summary:
ok = 3984
missing_validation_feature = 48

CD45isoform summary:
ok = 24
```

`n_cells` 加总检查：

```text
Adult_cHSPCs_Ultima_majority_voting = 13289
Standardpcaclusters = 13289
```

## 8. 本次已确认的缺失项

ADT/HTO summary 中缺失：

```text
mAOC433 / HTO-433
mAOC434 / HTO-434
```

RNA summary 中缺失：

```text
FADT0170 / Tag-ALFA -> RNA:Tag-ALFA
mAOC415 / mAOC-CD158b1 -> RNA:KIR2DL2
```

用户已确认这些缺失符合预期。

## 9. 需要手动检查的部分

每次重新运行 Step 3 后，建议检查：

```text
1. 每个 group_by 的 n_cells 加总是否等于 Seurat 对象总细胞数。
2. status/message 中的缺失项是否符合预期。
3. `mAOC-CD16` 是否使用 `RNA:FCGR3A`，`mAOC-HLA-ABC` 是否使用 `RNA:HLA-A`。
4. group_level 顺序是否满足后续绘图或统计需求。
5. 如果更换 mapping 表，应优先使用 `adt_feature_mapping.llm_review.tsv`，避免多基因 validation_feature 直接进入汇总。
```

## 10. 注意事项

```text
1. Step 3 不修改原始 Seurat RDS。
2. ADT 与 HTO 输出在同一个 `adt_group_summary.tsv` 中，通过 `assay` 字段区分。
3. RNA 和 CD45isoforms 验证表只包含 mapping 中对应 validation_assay 的条目。
4. 当前输出为长表格式，便于后续按 feature_id、group_by 和 group_level 连接其他 QC 结果。
```
