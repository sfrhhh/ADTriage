# Step 2: ADT Validation Feature Plots 可复现操作说明

## 1. 目标

Step 2 的目标是基于 Seurat RDS 对象和 Step 1 的 feature mapping 表，绘制 ADT feature 与对应验证 feature 的 UMAP FeatureDimPlot。

当前版本只输出 `plot_feature`，不输出 density 图。

输出目录：

```text
ADTriage_output/02_feature_plots/
```

核心输出：

```text
adt_validation_plot_feature_top175_page001.png
adt_validation_plot_feature_top175_page002.png
adt_validation_plot_feature_top175_page003.png
adt_validation_plot_feature_top175_page004.png
adt_validation_plot_feature_top175_page005.png
adt_validation_plot_feature_top175_page006.png
adt_validation_plot_feature_top175_page007.png
adt_validation_plot_feature_top175_page008.png
adt_validation_plot_layout_index_top175.tsv
adt_validation_plot_status_top175.tsv
panel_png_top175/plot_feature/
```

## 2. 输入文件

```text
tests/data/test_seurat_qc.rds
ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv
```

输入 Seurat 对象要求：

```text
RNA assay
ADT assay
CD45isoforms assay
可选 HTO assay
StandardpcaUMAP2D reduction
```

脚本会自动查找 reduction，优先匹配：

```text
standardpcaumap2d
standardumap2d
harmony.umap
cca.umap
umap
```

本次实际使用：

```text
StandardpcaUMAP2D
```

## 3. 使用脚本

```text
scripts/plot_adt_validation_features.R
```

该脚本使用 `scop::FeatureDimPlot()` 绘图，并通过 ImageMagick 拼接 PNG：

```text
/data/project/yanhuixu/frsui/software/bin/magick
```

## 4. 运行环境

R 使用 `scop` conda 环境：

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript
```

本步骤遵循 `scop-r` skill：使用本地 `scop` 包，并对 reduction 名称使用精确 regex。

## 5. 绘图逻辑

### 5.1 assay 归一化

脚本读取 RDS 后在内存中重新归一化，不覆盖原始 RDS。

```r
bm <- NormalizeData(bm, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = FALSE)
bm <- NormalizeData(bm, assay = "CD45isoforms", normalization.method = "CLR", verbose = FALSE)
bm <- NormalizeData(bm, assay = "HTO", normalization.method = "CLR", verbose = FALSE)
```

注意：

```text
ADT 使用 CLR margin=2。
CD45isoforms 使用 CLR 默认参数。
HTO 使用 CLR 默认参数，不使用 margin=2。
```

### 5.2 成对 panel

每个 mapping row 输出两个 panel：

```text
左侧：validation feature
右侧：ADT 或 HTO feature
```

普通 gene-level ADT：

```text
左侧 RNA:<validation_feature>
右侧 ADT:<feature_name>
```

CD45 isoform：

```text
左侧 CD45isoforms:<validation_feature>
右侧 ADT:<feature_name>
```

control / spike / hashtag：

```text
左侧 No validation feature 占位图
右侧 ADT 或 HTO feature
```

如果绘图 feature 不存在，脚本会在原位置输出红字错误占位图，并在 status 表记录原因。

### 5.3 关键绘图参数

本次最终配置：

```text
plot_type = plot_feature
panel_px = 800
panel_dpi = 200
pt.size = 0.5
raster = FALSE
pt.alpha = 0.8
pairs_per_row = 2
max_page_height_px = 9600
legend text/title size = 10
```

`pt.size` 必须显式固定。`scop::FeatureDimPlot()` 默认 `pt.size = NULL` 时会按细胞数自动计算点大小，本项目测试中会导致 top10 与 top175 批量绘图的单 panel 点大小不一致。

## 6. 复现命令

在项目根目录运行：

```bash
cd /data/project/yanhuixu/frsui/sc_project/ADTriage
```

### 6.1 R 脚本语法检查

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript -e 'parse(file="scripts/plot_adt_validation_features.R"); cat("R parse OK\n")'
```

### 6.2 全量绘图

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript scripts/plot_adt_validation_features.R \
  --seurat-rds tests/data/test_seurat_qc.rds \
  --mapping ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv \
  --outdir ADTriage_output/02_feature_plots \
  --n-features 175 \
  --panel-px 800 \
  --panel-dpi 200 \
  --pt-size 0.5 \
  --max-page-height-px 9600 \
  --magick /data/project/yanhuixu/frsui/software/bin/magick
```

## 7. 输出文件说明

### 7.1 分页大图

全量 175 个 feature 会生成 350 个 800x800 panel。单张总图高度会超过 ImageMagick 默认限制，因此脚本按页输出。

本次输出 8 页：

```text
ADTriage_output/02_feature_plots/adt_validation_plot_feature_top175_page001.png
...
ADTriage_output/02_feature_plots/adt_validation_plot_feature_top175_page008.png
```

页面尺寸：

```text
page001-page007: 3200x9600
page008: 3200x3200
```

### 7.2 单 panel PNG

目录：

```text
ADTriage_output/02_feature_plots/panel_png_top175/plot_feature/
```

文件命名格式：

```text
<pair_index>_<side>_<feature_id>_<assay_feature>.png
```

示例：

```text
009_R_FADT0042_ADT_mAOC-CD21.png
```

### 7.3 layout index 表

```text
ADTriage_output/02_feature_plots/adt_validation_plot_layout_index_top175.tsv
```

用途：

```text
记录每个 panel 在分页大图中的顺序、页码、行列位置、feature_id、feature_name、assay、feature、状态和对应单 panel 文件。
```

通过搜索该表可以定位某个 feature 在大图中的位置。

### 7.4 status 表

```text
ADTriage_output/02_feature_plots/adt_validation_plot_status_top175.tsv
```

用途：

```text
快速查看每个 panel 是否成功绘制，以及错误原因。
```

## 8. 本次运行统计

```text
feature rows = 175
panel PNG = 350
layout index rows = 350 + header
status rows = 350 + header
page PNG = 8
ok = 340
placeholder = 6
error = 4
```

本次仍为 error 的 panel：

```text
mAOC433 / HTO-433      Feature not found in HTO: HTO-433
mAOC434 / HTO-434      Feature not found in HTO: HTO-434
FADT0170 / Tag-ALFA    Feature not found in RNA: Tag-ALFA
mAOC415 / mAOC-CD158b1 Feature not found in RNA: KIR2DL2
```

已通过 Step 1 多基因收敛修复的条目：

```text
mAOC0301 / mAOC-CD16     RNA:FCGR3A, ADT:mAOC-CD16, status = ok
mAOC0453 / mAOC-HLA-ABC  RNA:HLA-A, ADT:mAOC-HLA-ABC, status = ok
```

## 9. 需要手动检查的部分

每次重新运行 Step 2 后，建议检查：

```text
1. page001 和 page008 是否能正常打开，排版是否为每行两对图。
2. `adt_validation_plot_layout_index_top175.tsv` 中目标 feature 的页码和行列是否可检索。
3. `adt_validation_plot_status_top175.tsv` 中 error 项是否符合预期。
4. CD45 isoform 是否使用 CD45isoforms assay，而不是 RNA:PTPRC。
5. HTO 缺失项是否只是 Seurat 对象缺 feature，而不是脚本中 assay 选择错误。
6. 多基因 symbol 是否已经在 Step 1 的 llm_review 表中变成单一 RNA validation_feature。
7. 关键视觉回归 panel `009_R_FADT0042_ADT_mAOC-CD21.png` 的点大小是否保持当前效果。
```

## 10. 注意事项

```text
1. 不要使用 `raster = TRUE`，当前图中点会变得难以观察。
2. 不再输出 density 图；如果后续恢复 density，需要重新评估绘图效果。
3. 单张超高总图会触发 ImageMagick 限制，因此保留分页输出。
4. top10 测试文件只用于调参，正式 Step 2 输出不保留 top10 文件。
5. 绘图输入应使用 `adt_feature_mapping.llm_review.tsv`，不要直接用 initial 表。
```
