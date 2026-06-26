# ADTriage

**ADTriage: an ADT feature-level triage workflow for CITE-seq quality assessment**

ADTriage 是一个面向 CITE-seq / Feature Barcode 数据的 **ADT 特征级质量评估工作流**。它以质控后的 Seurat 对象和 Cell Ranger Feature Barcode reference 文件为输入，逐个评估 Antibody-Derived Tag（ADT）特征是否具有合理的生物学信号、是否与 RNA 或 CD45 isoform 信息一致、是否符合细胞类型先验表达模式，并输出可人工审查的 ADT 质量标签。

---

## 1. 项目目标

CITE-seq 数据中的 ADT 信号受抗体质量、panel 设计、染色条件、背景噪音、ambient antibody、isotype control、feature reference 配置和细胞类型组成影响。常规单细胞分析流程通常关注细胞层面的质控，而缺少对 **每一个 ADT feature 是否可信** 的系统审查。

ADTriage 的目标是建立一个可复现、可审计、可人工复核的工作流，用于回答以下问题：

1. 某个 ADT feature 是否能映射到明确的人类基因符号或特定 CD45 isoform？
2. 该 ADT 的 CLR 信号是否与对应 RNA 表达或 isoform assay 信号一致？
3. 该 ADT 的阳性细胞类型是否符合已知免疫学或造血系统 marker 先验知识？
4. 该 ADT 的信号分布是否具备可分离的阳性/阴性群体？
5. 该 ADT 应被标记为 `correct`、`suspicious`、`wrong`、`no_signal`、`control` 还是 `manual_review`？

---

## 2. 输入文件

ADTriage 需要两个核心输入。

### 2.1 质控后的 Seurat 对象

输入对象必须为 `.rds` 格式的 Seurat object，并至少包含以下三个 assay：

```text
RNA
ADT
CD45-isoform
```

其中：

* `RNA` assay 用于验证 ADT 对应基因的转录本表达。
* `ADT` assay 用于评估抗体标签信号。
* `CD45-isoform` assay 用于验证 CD45RA、CD45RO 等无法直接映射为独立 gene symbol 的 CD45 剪接变体信号。

Seurat 对象中还需要包含细胞聚类信息或细胞类型注释信息。若对象中尚未包含细胞类型注释，ADTriage 可调用 CellTypist 进行 cluster-level 和 single-cell-level 注释。

### 2.2 Feature Barcode reference 文件

第二个输入为 Cell Ranger 上游分析时使用的 `FB_ref.csv` 文件。该文件应符合 10x Genomics Feature Barcode reference 格式，至少包含以下列：

```text
id
name
read
pattern
sequence
feature_type
```

ADTriage 会读取 `id` 和 `name` 字段，构建 ADT feature 与 gene symbol、CD45 isoform 或 control feature 的映射关系。

---

## 3. 输出结果

ADTriage 的主要输出为一组逐步递进的表格和图形文件。每一步均生成新文件，不直接覆盖上一步结果。

### 3.1 ADT feature 映射表

```text
adt_feature_mapping.tsv
```

该表记录每个 ADT feature 的基础注释信息。

推荐字段如下：

```text
feature_id
feature_name
feature_type
target_class
human_gene_symbol
validation_assay
validation_feature
mapping_method
mapping_confidence
mapping_evidence
needs_manual_review
```

其中 `target_class` 可取以下值：

```text
gene
cd45_isoform
control
hashtag
unknown
ambiguous
```

### 3.2 RNA–ADT 相关性结果表

```text
rna_adt_correlation_summary.tsv
```

该表记录每个 ADT 与其对应 RNA gene 或 CD45 isoform feature 的信号一致性。

推荐字段如下：

```text
feature_id
feature_name
validation_feature
annotation_level
n_celltypes
min_cells_per_type
pearson_r
spearman_r
linear_model_slope
linear_model_intercept
r_squared
correlation_status
```

其中 `annotation_level` 包括：

```text
cluster_level
single_cell_level
```

### 3.3 ADT 分布指标表

```text
adt_distribution_metrics.tsv
```

该表记录每个 ADT 的 CLR 分布、双峰拟合、阳性率和信噪比等指标。

推荐字段如下：

```text
feature_id
feature_name
n_components_selected
mu_negative
mu_positive
sd_negative
sd_positive
threshold
positive_rate
separation_index
overlap_area
snr_raw
tnr_raw
dynamic_range_clr
positive_cv
distribution_status
```

### 3.4 Marker 先验表达审查表

```text
marker_prior_review.tsv
```

该表记录每个 ADT 理论上应表达的细胞类型及其与观测结果的一致性。

推荐字段如下：

```text
feature_id
feature_name
expected_positive_celltypes
expected_negative_celltypes
observed_positive_celltypes
observed_negative_celltypes
prior_confidence
prior_evidence_level
prior_consistency_score
needs_manual_review
review_reason
```

### 3.5 最终 ADT 质量标签表

```text
final_adt_qc_labels.tsv
```

该表为 ADTriage 的核心输出，整合映射、相关性、分布和先验表达信息，为每个 ADT feature 分配最终标签。

推荐字段如下：

```text
feature_id
feature_name
target_class
human_gene_symbol
validation_assay
validation_feature
mapping_confidence
cluster_spearman_r
single_cell_spearman_r
separation_index
snr_raw
positive_rate
prior_consistency_score
final_label
final_confidence
review_reason
```

`final_label` 可取以下值：

```text
correct
suspicious
wrong
no_signal
control
manual_review
```

---

## 4. 工作流概览

ADTriage 的分析流程分为 9 个步骤。

### 4.1 输入检查

检查 Seurat 对象和 `FB_ref.csv` 是否满足运行要求。

检查内容包括：

```text
Seurat object 是否可读取
RNA assay 是否存在
ADT assay 是否存在
CD45-isoform assay 是否存在
ADT feature 是否能与 FB_ref.csv 对应
counts/data slot 是否存在
metadata 中是否已有 cluster 或 cell type 信息
```

### 4.2 构建 ADT 与验证特征映射表

读取 `FB_ref.csv`，根据 ADT 名称构建 feature-level 映射表。

映射规则包括：

1. 普通 ADT marker 映射到对应 human gene symbol。
2. CD45 映射到 `PTPRC`。
3. CD45RA、CD45RO、CD45RB、CD45RC 等 CD45 isoform 不映射为独立 gene symbol，而是映射到 `CD45-isoform` assay 中的对应 feature。
4. isotype control、IgG control、Fc control、spike-in 等特征标记为 `control`。
5. 无法自动判断的特征标记为 `ambiguous` 或 `unknown`，进入人工审查。

### 4.3 运行标准化处理

使用 `scop` R 包中的 `standard_scop` 函数对 Seurat 对象进行标准化处理，生成用于后续分析的新 Seurat 对象。

### 4.4 运行 CellTypist 注释

调用 `RunCellTypist`，基于指定 `.pkl` 模型文件进行细胞类型注释。

ADTriage 同时保存两类注释结果：

```text
cluster-level annotation
single-cell-level annotation
```

两类注释结果分别用于后续 RNA–ADT 相关性分析和 marker 先验表达一致性分析。

### 4.5 绘制配对 FeaturePlot

对每个 ADT feature 绘制配对 FeaturePlot。

常见图形组合包括：

```text
ADT CLR signal
RNA log-normalized signal
CD45 isoform signal
CellTypist cluster-level annotation
CellTypist single-cell-level annotation
```

每个 ADT 输出独立图形，所有图形合并保存为 PDF。

### 4.6 计算 RNA–ADT 或 isoform–ADT 相关性

基于细胞类型进行 pseudobulk 聚合。

每个散点代表一种细胞类型：

```text
x = mean RNA log-normalized expression 或 mean CD45 isoform signal
y = mean ADT CLR signal
```

分别基于 cluster-level 和 single-cell-level 注释计算：

```text
Pearson correlation
Spearman correlation
linear model slope
R²
n_celltypes
```

### 4.7 构建 marker 先验表达审查结果

根据每个 ADT 的 target 和实际出现的细胞类型，判断该 marker 理论上应表达于哪些细胞类型。

该步骤输出：

```text
expected_positive_celltypes
expected_negative_celltypes
observed_positive_celltypes
prior_consistency_score
prior_confidence
manual_review flag
```

对于置信度低于阈值的判断，标记为人工审查。

默认阈值：

```text
prior_confidence < 0.80 -> manual_review
```

### 4.8 评估 ADT CLR 分布指标

对每个 ADT 的 CLR 值分布进行建模，评估是否存在可分离的阴性群体和阳性群体。

推荐指标包括：

```text
positive_rate
separation_index
overlap_area
snr_raw
tnr_raw
dynamic_range_clr
positive_cv
```

如果 ADT 分布不具备可分离信号，则标记为 `no_signal` 或 `suspicious`。

### 4.9 综合打分与最终标签

整合以下信息：

```text
feature mapping confidence
RNA–ADT correlation
isoform–ADT correlation
ADT distribution metrics
cell type prior consistency
control feature status
manual review flag
```

最终为每个 ADT feature 分配标签：

```text
correct
suspicious
wrong
no_signal
control
manual_review
```

---

## 5. 推荐项目结构

```text
ADTriage/
├── README.md
├── config/
│   ├── config.yaml
│   ├── marker_prior.yaml
│   ├── cd45_isoform_map.yaml
│   └── scoring_rules.yaml
├── R/
│   ├── 00_validate_input.R
│   ├── 01_standard_scop.R
│   ├── 02_run_celltypist.R
│   ├── 03_featureplot_pairs.R
│   ├── 04_rna_adt_correlation.R
│   ├── 05_adt_distribution_metrics.R
│   └── 06_final_scoring.R
├── python/
│   ├── parse_fb_ref.py
│   ├── hgnc_query.py
│   └── marker_prior_review.py
├── prompts/
│   ├── adt_mapping_review.md
│   └── marker_prior_review.md
├── schemas/
│   ├── adt_mapping.schema.json
│   ├── marker_prior.schema.json
│   └── final_qc.schema.json
├── outputs/
│   ├── tables/
│   ├── figures/
│   ├── reports/
│   └── objects/
└── SKILL.md
```

---

## 6. 配置文件示例

推荐使用 `config/config.yaml` 统一管理输入、输出、assay 名称、阈值和工具路径。

```yaml
project:
  name: "ADTriage"
  species: "human"
  genome: "GRCh38"

input:
  seurat_rds: "data/qc_seurat.rds"
  fb_ref_csv: "data/FB_ref.csv"
  celltypist_model_pkl: "ref/celltypist_model.pkl"

assays:
  rna_assay: "RNA"
  adt_assay: "ADT"
  cd45_isoform_assay: "CD45-isoform"

metadata:
  cluster_col: "seurat_clusters"
  celltypist_cluster_col: "celltypist_cluster_label"
  celltypist_single_cell_col: "celltypist_cell_label"

tools:
  Rscript: "/path/to/Rscript"
  python: "/path/to/python"
  celltypist: "/path/to/celltypist"

reference:
  hgnc_tsv: "ref/hgnc_complete_set.tsv"
  marker_prior_yaml: "config/marker_prior.yaml"
  cd45_isoform_map_yaml: "config/cd45_isoform_map.yaml"
  scoring_rules_yaml: "config/scoring_rules.yaml"

output:
  outdir: "outputs"
  table_dir: "outputs/tables"
  figure_dir: "outputs/figures"
  report_dir: "outputs/reports"
  object_dir: "outputs/objects"
  cache_dir: "outputs/cache"

mapping:
  confidence_cutoff_manual_review: 0.80
  use_hgnc: true
  use_llm_for_ambiguous_only: true

cell_annotation:
  run_celltypist: true
  min_cells_per_type: 30

correlation:
  aggregation: "celltype_mean"
  methods:
    - "pearson"
    - "spearman"
  min_celltypes_for_correlation: 5

distribution:
  method: "gaussian_mixture"
  max_components: 3
  min_positive_rate: 0.005

final_scoring:
  manual_review_cutoff: 0.80
  labels:
    - "correct"
    - "suspicious"
    - "wrong"
    - "no_signal"
    - "control"
    - "manual_review"
```

---

## 7. 依赖工具

### 7.1 R 依赖

```r
Seurat
Matrix
data.table
dplyr
ggplot2
patchwork
cowplot
future
mclust
jsonlite
yaml
scop
```

可选依赖：

```r
dsb
ThresholdR
ADTnorm
CITESeQC
mixtools
flexmix
pROC
rmarkdown
```

### 7.2 Python 依赖

```python
pandas
numpy
scipy
anndata
scanpy
celltypist
pyyaml
requests
```

可选依赖：

```python
mygene
scvi-tools
```

---

## 8. 运行方式

ADTriage 推荐按步骤运行，而不是一次性黑箱执行。

### 8.1 检查输入

```bash
Rscript R/00_validate_input.R \
  --config config/config.yaml
```

预期输出：

```text
outputs/tables/input_validation_summary.tsv
```

### 8.2 构建 ADT feature 映射表

```bash
python python/parse_fb_ref.py \
  --config config/config.yaml
```

预期输出：

```text
outputs/tables/adt_feature_mapping.tsv
```

### 8.3 运行标准化处理

```bash
Rscript R/01_standard_scop.R \
  --config config/config.yaml
```

预期输出：

```text
outputs/objects/seurat_standard_scop.rds
```

### 8.4 运行 CellTypist 注释

```bash
Rscript R/02_run_celltypist.R \
  --config config/config.yaml
```

预期输出：

```text
outputs/objects/seurat_celltypist_annotated.rds
outputs/tables/celltypist_annotation_summary.tsv
```

### 8.5 绘制配对 FeaturePlot

```bash
Rscript R/03_featureplot_pairs.R \
  --config config/config.yaml
```

预期输出：

```text
outputs/figures/adt_featureplot_pairs.pdf
```

### 8.6 计算 RNA–ADT 相关性

```bash
Rscript R/04_rna_adt_correlation.R \
  --config config/config.yaml
```

预期输出：

```text
outputs/tables/rna_adt_correlation_summary.tsv
outputs/figures/rna_adt_correlation_scatter.pdf
```

### 8.7 计算 ADT 分布指标

```bash
Rscript R/05_adt_distribution_metrics.R \
  --config config/config.yaml
```

预期输出：

```text
outputs/tables/adt_distribution_metrics.tsv
outputs/figures/adt_distribution_density.pdf
```

### 8.8 综合打分

```bash
Rscript R/06_final_scoring.R \
  --config config/config.yaml
```

预期输出：

```text
outputs/tables/final_adt_qc_labels.tsv
```

---

## 9. 标签定义

### 9.1 correct

该 ADT feature 满足以下主要条件：

```text
映射关系明确
分布存在可识别阳性群体
阳性细胞类型符合 marker 先验表达
RNA–ADT 或 isoform–ADT 信号具有一致性
无明显背景异常
```

### 9.2 suspicious

该 ADT feature 存在部分异常，但不足以直接判定为错误。

常见情况包括：

```text
分布可分离，但阳性细胞类型不完全符合预期
RNA–ADT 相关性低，但 marker 本身可能存在 RNA/protein discordance
positive_rate 异常偏高或偏低
不同注释层级下结论不一致
```

### 9.3 wrong

该 ADT feature 的观测信号与预期明显冲突。

常见情况包括：

```text
阳性细胞类型与已知 marker 生物学背景冲突
对应 RNA 或 isoform 信号完全不支持 ADT 分布
feature 名称或 reference 映射存在高风险错误
```

### 9.4 no_signal

该 ADT feature 缺乏可检测信号。

常见情况包括：

```text
CLR 分布接近单峰背景
positive_rate 低于阈值
dynamic range 低
SNR 低
无明确阳性细胞群
```

### 9.5 control

该 feature 为实验对照或技术对照。

常见情况包括：

```text
isotype control
IgG control
Fc control
spike-in
hashtag
non-targeting control
```

### 9.6 manual_review

该 feature 需要人工审查。

触发条件包括：

```text
mapping_confidence < 0.80
prior_confidence < 0.80
feature_name 无法自动解析
marker 生物学背景存在歧义
自动规则之间结论冲突
```

---

## 10. 设计原则

ADTriage 遵循以下原则：

1. **每一步生成新表格，不覆盖上一步结果。**
2. **所有标签必须可追溯到具体指标和规则。**
3. **不使用 RNA–ADT correlation 单独判定 ADT 错误。**
4. **CD45 isoform 不应强制映射为独立 human gene symbol。**
5. **LLM 只用于模糊命名解析和 marker 先验表达审查，不读取原始表达矩阵。**
6. **所有低置信度结果进入人工审查，而不是自动删除。**
7. **最终结果用于 triage，不直接等同于统计学显著性检验结论。**

---

## 11. 适用场景

ADTriage 适用于以下数据类型：

```text
CITE-seq
Feature Barcode
TotalSeq-A/B/C panel
含 ADT assay 的 Seurat object
含 CD45 isoform antibody panel 的单细胞数据
多抗体 panel 混样设计
```

尤其适用于以下问题：

```text
判断 ADT panel 中哪些抗体信号可信
检查 FB_ref.csv 是否存在 feature 命名或映射问题
评估 CD45RA/CD45RO 等 isoform 抗体是否合理
识别无信号或背景异常的 ADT feature
为下游 ADT gating、WNN、marker-based annotation 提供质量依据
```

---

## 12. 局限性

ADTriage 的结果依赖于输入数据质量和参考信息完整性。

主要限制包括：

1. **RNA–protein discordance**
   RNA 表达与表面蛋白表达并非总是一致。低相关性不等同于 ADT 错误。

2. **细胞类型注释误差**
   CellTypist 或 cluster annotation 错误会影响 marker prior consistency 评估。

3. **稀有细胞类型不足**
   当某个 marker 的理论阳性细胞类型在数据中细胞数不足时，无法可靠评估该 ADT 是否合理。

4. **control feature 依赖实验设计**
   spike-in、IgG、Fc block、isotype 等对照需要结合具体实验设计解释。

5. **缺少训练数据时无法建立监督分类器**
   第一版 ADTriage 使用规则和指标打分进行 triage，最终标签仍需人工审查确认。

---

## 13. 推荐审查流程

完成自动分析后，建议按以下顺序人工审查：

1. 先检查 `manual_review` 和 `ambiguous` feature。
2. 再检查 `wrong` 和 `suspicious` feature。
3. 对 `no_signal` feature 查看 CLR density 和 FeaturePlot。
4. 对 CD45 isoform feature 单独检查 `CD45-isoform` assay。
5. 最后检查 `correct` feature 是否存在明显假阳性。

人工审查后建议新增一列：

```text
manual_label
```

并保留自动标签：

```text
final_label
```

不要直接覆盖自动结果。这样可以保留自动规则与人工判断之间的差异。

---

## 14. 最小可用版本

ADTriage 的最小可用版本应至少输出以下 5 个文件：

```text
outputs/tables/adt_feature_mapping.tsv
outputs/tables/rna_adt_correlation_summary.tsv
outputs/tables/adt_distribution_metrics.tsv
outputs/tables/marker_prior_review.tsv
outputs/tables/final_adt_qc_labels.tsv
```

这 5 个文件构成 ADT feature-level triage 的核心结果。后续版本可继续增加 PDF report、HTML report、交互式审查界面和监督学习模型。
