# Step 4: Signal Consistency 可复现操作说明

## 1. 目标

Step 4 的目标是基于 Step 3 的 group-level pseudobulk summary，计算 ADT 信号与验证特征信号在不同 group level 间的一致性。

普通 gene-level ADT：

```text
ADT CLR group mean
vs
RNA log-normalized group mean
```

CD45 isoform ADT：

```text
ADT CLR group mean
vs
CD45isoforms group mean
```

输出目录：

```text
ADTriage_output/04_signal_consistency/
```

核心输出：

```text
rna_adt_correlation_summary.tsv
cd45isoform_adt_correlation_summary.tsv
rna_adt_<group_by>_consistency_scatter_plots.pdf
cd45isoform_adt_<group_by>_consistency_scatter_plots.pdf
```

## 2. 输入文件

来自 Step 3：

```text
ADTriage_output/03_group_summary/rna_group_summary.tsv
ADTriage_output/03_group_summary/cd45isoform_group_summary.tsv
```

本次输入包含两个 `group_by`：

```text
Adult_cHSPCs_Ultima_majority_voting
Standardpcaclusters
```

## 3. 使用脚本

```text
scripts/calculate_signal_consistency.R
```

脚本功能：

```text
1. 读取 Step 3 的 RNA 和 CD45isoform group summary。
2. 对每个 group_by + feature 取 status = ok 且信号值完整的 group level。
3. 计算 Pearson correlation、Spearman correlation、linear model slope、linear model R²。
4. 记录 n_groups 和 min_cells_per_group。
5. 按 group_by 分别输出 scatter plot PDF。
```

## 4. 统计规则

每个 feature 在每个 `group_by` 下独立计算。

使用的统计量：

```text
pearson_correlation
spearman_correlation
linear_model_slope
linear_model_r_squared
n_groups
min_cells_per_group
```

默认要求：

```text
min_groups = 3
```

若完整 group 数不足 3：

```text
status = insufficient_groups
相关统计量写为 NA
```

线性模型：

```r
lm(mean_ADT_CLR ~ validation_group_mean)
```

其中 scatter plot：

```text
x axis = validation feature group mean
y axis = ADT CLR group mean
point = group level
point size = n_cells
orange line = least-squares linear fit
subtitle = Spearman r, Pearson r, R², n_groups
```

## 5. 复现命令

在项目根目录运行：

```bash
cd /data/project/yanhuixu/frsui/sc_project/ADTriage
```

### 5.1 R 脚本语法检查

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript -e 'parse(file="scripts/calculate_signal_consistency.R"); cat("R parse OK\n")'
```

### 5.2 计算信号一致性

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript scripts/calculate_signal_consistency.R \
  --group-summary-dir ADTriage_output/03_group_summary \
  --outdir ADTriage_output/04_signal_consistency \
  --min-groups 3
```

## 6. 输出文件说明

### 6.1 rna_adt_correlation_summary.tsv

普通 gene-level ADT 与 RNA 验证特征的 group-level 一致性统计。

关键字段：

```text
group_by
feature_id
feature_name
target_class
validation_assay
validation_feature
validation_signal
pearson_correlation
spearman_correlation
linear_model_slope
linear_model_r_squared
n_groups
min_cells_per_group
status
message
```

### 6.2 cd45isoform_adt_correlation_summary.tsv

CD45 isoform ADT 与 CD45isoforms assay 验证特征的 group-level 一致性统计。

字段同 `rna_adt_correlation_summary.tsv`。

### 6.3 scatter plot PDF

由于 Step 3 可包含多个 metadata 分组列，Step 4 每个 `group_by` 单独输出 PDF。

本次输出：

```text
rna_adt_Adult_cHSPCs_Ultima_majority_voting_consistency_scatter_plots.pdf
rna_adt_Standardpcaclusters_consistency_scatter_plots.pdf
cd45isoform_adt_Adult_cHSPCs_Ultima_majority_voting_consistency_scatter_plots.pdf
cd45isoform_adt_Standardpcaclusters_consistency_scatter_plots.pdf
```

页面数：

```text
RNA PDF: 每个 group_by 168 页
CD45isoform PDF: 每个 group_by 1 页
```

## 7. 本次运行统计

输出表：

```text
rna_adt_correlation_summary.tsv rows = 336 + header
cd45isoform_adt_correlation_summary.tsv rows = 2 + header
```

状态统计：

```text
RNA:
ok = 332
insufficient_groups = 4

CD45isoform:
ok = 2
```

`insufficient_groups` 条目：

```text
Adult_cHSPCs_Ultima_majority_voting / FADT0170 / Tag-ALFA / RNA:Tag-ALFA
Adult_cHSPCs_Ultima_majority_voting / mAOC415 / mAOC-CD158b1 / RNA:KIR2DL2
Standardpcaclusters / FADT0170 / Tag-ALFA / RNA:Tag-ALFA
Standardpcaclusters / mAOC415 / mAOC-CD158b1 / RNA:KIR2DL2
```

这些条目与 Step 3 中的缺失验证特征一致。

CD45 isoform 结果：

```text
Adult_cHSPCs_Ultima_majority_voting / mAOC-CD45RA / PTPRC-RA: ok
Standardpcaclusters / mAOC-CD45RA / PTPRC-RA: ok
```

## 8. 验证记录

已完成：

```text
1. R parse 检查通过。
2. Step 4 完整运行成功。
3. 输出 2 个 summary TSV 和 4 个 PDF。
4. summary 行数与 Step 3 feature 数和 group_by 数一致。
5. PDF 文件均已生成且大小非 0。
```

说明：

```text
系统没有 pdfinfo 命令，因此未用 pdfinfo 做页数检查。
R 脚本运行日志显示 RNA PDF 每个 168 页，CD45isoform PDF 每个 1 页。
```

## 9. 需要手动检查的部分

每次重新运行 Step 4 后，建议检查：

```text
1. summary 表中 status != ok 的条目是否可由 Step 3 的缺失 feature 解释。
2. 随机打开每个 group_by 的 RNA PDF，确认 x/y 轴、拟合线和统计量标注正确。
3. 检查 CD45isoform PDF 中 CD45RA 与 PTPRC-RA 的趋势是否符合预期。
4. 如果后续新增 group_by，确认是否生成对应 group_by 文件名的独立 PDF。
5. 关注 n_groups 很少或 min_cells_per_group 很低的 feature，相关性解释应谨慎。
```

## 10. 注意事项

```text
1. Step 4 不重新读取 Seurat 对象，只依赖 Step 3 的 group summary。
2. 相关性只反映 group-level 平均信号一致性，不代表 single-cell 层面的相关性。
3. RNA 低表达或 validation feature 缺失会导致 insufficient_groups 或不稳定相关性。
4. 如果 Step 3 的 group_by、mapping 或归一化规则改变，必须重新运行 Step 4。
```
