# ADTriage

**ADTriage: an ADT feature-level triage workflow for CITE-seq quality assessment**

ADTriage 是一个面向 CITE-seq / Feature Barcode 数据的 **ADT feature-level quality triage workflow**。
其目标不是重新完成单细胞上游分析、细胞注释或整合分析，而是针对每一个 Antibody-Derived Tag（ADT）特征，整合以下信息后生成可人工审查的质量评估结果：

```text
1. FB_ref.csv 中的 feature barcode 信息
2. mAOC_gene_symbol_map.csv 中的 ADT id → gene symbol 映射
3. Seurat 对象中的 RNA assay、ADT assay、CD45isoforms assay
4. Seurat metadata 中用户指定的 >=1 个分组变量
5. ADT CLR 信号与 RNA log-normalized 信号的一致性
6. ADT CLR 分布形态、信噪比、动态范围、阳性比例等指标
7. CD45 剪接变体特异性信号
8. 规则判断、大模型辅助判断与人工审查结果
```

ADTriage 输出的核心结果是每个 ADT feature 的质量标签，例如：

```text
correct
suspicious
wrong
no_signal
control
manual_review
```

---

## 1. 项目定位

ADTriage 关注的是：

```text
一个 ADT feature 本身是否合理？
它的信号是否可解释？
它是否与 RNA、CD45 isoform 或细胞分组结果一致？
它是否需要人工审查？
```

ADTriage 不负责以下任务：

```text
1. 不运行 standard_scop
2. 不运行 CellTypist
3. 不执行细胞注释
4. 不执行 Seurat 对象标准化主流程
5. 不执行整合分析
6. 不决定最终生物学结论
```

用户需要在输入的 Seurat 对象中提前准备好用于分组的 metadata 列。
ADTriage 只读取这些列，并基于这些列绘制分组图、计算分组统计量和生成审查表。

---

## 2. 输入文件

### 2.1 质控后的 Seurat 对象

输入 Seurat 对象必须至少包含以下 assay：

```text
RNA
ADT
CD45isoforms
```

其中：

```text
RNA           : 基因表达矩阵
ADT           : ADT feature barcode 矩阵
CD45isoforms  : CD45 剪接变体特征矩阵
```

`CD45isoforms` assay 必须包含以下四个 feature：

```text
PTPRC-RA
PTPRC-RB
PTPRC-RC
PTPRC-RO
```

Seurat metadata 中必须包含至少 **1** 个用户指定的分组列，例如：

```text
sample
seurat_clusters
manual_annotation
treatment
disease_status
celltype_major
celltype_minor
```

用户指定多少个分组列，ADTriage 就会针对多少种分组方式绘图和计算统计结果。

示例配置：

```yaml
input:
  seurat_rds: "data/qc_seurat.rds"

assays:
  rna_assay: "RNA"
  adt_assay: "ADT"
  cd45_isoform_assay: "CD45isoforms"

group_by:
  - "seurat_clusters"
  - "manual_annotation"
  - "sample"
```

---

### 2.2 FB_ref.csv

`FB_ref.csv` 是运行 Cell Ranger 上游分析时使用的 feature reference 文件。
该文件至少需要包含以下列：

```text
id
name
read
pattern
sequence
feature_type
```

示例：

```csv
id,name,read,pattern,sequence,feature_type
FADT0038,mAOC-CD8a,R2,5PNNNNNNNNNN,ACGT...,Antibody Capture
FADT0039,mAOC-CD3,R2,5PNNNNNNNNNN,TGCA...,Antibody Capture
FADT0040,mAOC-CD45RA,R2,5PNNNNNNNNNN,GGCC...,Antibody Capture
```

其中 `id` 是 ADTriage 中最重要的关联键。
后续 `mAOC_gene_symbol_map.csv` 会通过 `id` 与 `FB_ref.csv` 进行关联。

---

### 2.3 mAOC_gene_symbol_map.csv

`mAOC_gene_symbol_map.csv` 是用户额外提供的 ADT id 到 gene symbol 的映射表。
该表第一列 `id` 必须与 `FB_ref.csv` 中的 `id` 对应。

必须包含以下列：

```text
id
genesymbol
genesymbol_source
notes
last_seen_name
last_seen_reference
updated_at
```

示例：

```csv
id,genesymbol,genesymbol_source,notes,last_seen_name,last_seen_reference,updated_at
FADT0038,CD8A,reference:对应基因名,,mAOC-CD8a,抗体信息对照表260611-FY.xlsx,
FADT0039,CD3E,reference:对应基因名,,mAOC-CD3,抗体信息对照表260611-FY.xlsx,
FADT0040,PTPRC,reference:对应基因名,,mAOC-CD45RA,抗体信息对照表260611-FY.xlsx,
FADT0041,CD80,reference:对应基因名,,mAOC-CD80,抗体信息对照表260611-FY.xlsx,
FADT0042,CR2,reference:对应基因名,,mAOC-CD21,抗体信息对照表260611-FY.xlsx,2026-06-24 11:22:42
```

该文件的作用是：

```text
1. 优先通过 id 关联 ADT feature 与 gene symbol
2. 减少大模型逐行判断 ADT 名称的 token 消耗
3. 保留 gene symbol 来源、历史名称和人工备注
4. 为后续人工审查提供可追溯证据
```

---

## 3. ADT feature 映射规则

ADTriage 的第一步是构建 ADT feature 映射表。

映射优先级如下：

```text
1. 通过 FB_ref.csv$id 与 mAOC_gene_symbol_map.csv$id 关联 gene symbol
2. 对 CD45 / CD45 isoform 相关特征应用规则判断
3. 对空字段、无法识别字段、对照抗体字段应用规则判断
4. 仅对 unresolved / ambiguous 条目调用大模型
5. 大模型判断结果进入人工审查表
```

---

### 3.1 普通 ADT feature

对于可以明确映射到 human gene symbol 的 ADT，例如：

```text
CD3   -> CD3E / CD3D / CD3G，依据抗体实际 target 决定
CD4   -> CD4
CD8a  -> CD8A
CD14  -> CD14
CD19  -> CD19
CD34  -> CD34
CD38  -> CD38
CD80  -> CD80
CD21  -> CR2
```

此类 feature 的验证方式为：

```text
ADT assay 中的对应 ADT feature
vs
RNA assay 中的对应 gene symbol
```

输出字段中：

```text
target_class = gene
validation_assay = RNA
validation_feature = 对应 human gene symbol
```

---

### 3.2 CD45 总表达

如果 ADT feature 是 CD45 总表达，而不是 CD45RA、CD45RB、CD45RC、CD45RO 等剪接变体，则对应 gene symbol 为：

```text
PTPRC
```

此类 feature 可以使用 RNA assay 中的 `PTPRC` 作为验证特征。

输出字段中：

```text
target_class = gene
human_gene_symbol = PTPRC
validation_assay = RNA
validation_feature = PTPRC
```

---

### 3.3 CD45 剪接变体

如果 ADT feature 是 CD45 剪接变体，例如：

```text
CD45RA
CD45RB
CD45RC
CD45RO
```

则不应将其简单视为普通 `PTPRC` gene-level RNA 表达验证。

此类 feature 的 gene symbol 可以记录为：

```text
PTPRC
```

但其主要验证特征应来自 `CD45isoforms` assay。

映射规则：

```text
CD45RA -> CD45isoforms: PTPRC-RA
CD45RB -> CD45isoforms: PTPRC-RB
CD45RC -> CD45isoforms: PTPRC-RC
CD45RO -> CD45isoforms: PTPRC-RO
```

输出字段中：

```text
target_class = cd45_isoform
human_gene_symbol = PTPRC
validation_assay = CD45isoforms
validation_feature = PTPRC-RA / PTPRC-RB / PTPRC-RC / PTPRC-RO
```

注意：
CD45 isoform feature 的 RNA–ADT correlation 不应使用 `RNA:PTPRC` 作为唯一判定依据。
`RNA:PTPRC` 代表 PTPRC gene-level 总表达，不能区分 RA、RB、RC、RO isoform epitope。

---

### 3.4 Spike-in、IgFc、isotype 和其他对照

对于以下类型：

```text
SPIKE
Spike
IgG
IgFc
isotype
control
negative control
hashing control
background control
```

此类 feature 通常不应强制映射到 human gene symbol。

输出字段中：

```text
target_class = control
human_gene_symbol = NA
validation_assay = none
validation_feature = NA
needs_manual_review = TRUE 或 FALSE，取决于用户规则
```

如果用户提供了 `genesymbol`，例如：

```text
SPIKE-RBD1-S2E12
SPIKE-NTD-4A8
SPIKE-FP-76E1
```

则 ADTriage 保留该字段，但不会将其视为 human gene symbol。

建议输出：

```text
target_class = spike_control
human_gene_symbol = NA
external_target_name = SPIKE-RBD1-S2E12
validation_assay = none
```

---

### 3.5 无法判断或歧义 feature

如果 feature 无法通过以下方式确定：

```text
1. mAOC_gene_symbol_map.csv
2. CD45 isoform 规则
3. control / spike / isotype 规则
4. 本地 marker dictionary
```

则进入大模型辅助判断。

大模型只接收压缩后的信息：

```text
feature_id
feature_name
last_seen_name
genesymbol
genesymbol_source
notes
feature_type
```

大模型输出必须为结构化 JSON，例如：

```json
{
  "feature_id": "FADT0000",
  "feature_name": "mAOC-unknown",
  "target_class": "ambiguous",
  "human_gene_symbol": null,
  "validation_assay": "manual",
  "validation_feature": null,
  "confidence": 0.42,
  "reason": "Feature name does not provide enough information to infer a unique human gene symbol.",
  "needs_manual_review": true
}
```

当 `confidence < 0.8` 时，必须进入人工审查。

---

## 4. 输出目录结构

推荐输出目录：

```text
ADTriage_output/
├── 00_input_validation/
│   ├── input_validation_report.tsv
│   └── assay_feature_summary.tsv
├── 01_feature_mapping/
│   ├── adt_feature_mapping.initial.tsv
│   ├── adt_feature_mapping.llm_review.tsv
│   └── adt_feature_mapping.review_ready.tsv
├── 02_feature_plots/
│   ├── group_by_seurat_clusters/
│   ├── group_by_manual_annotation/
│   └── group_by_sample/
├── 03_group_summary/
│   ├── adt_group_summary.tsv
│   ├── rna_group_summary.tsv
│   └── cd45isoform_group_summary.tsv
├── 04_signal_consistency/
│   ├── rna_adt_correlation_summary.tsv
│   ├── cd45isoform_adt_correlation_summary.tsv
│   └── consistency_scatter_plots/
├── 05_distribution_metrics/
│   ├── adt_distribution_metrics.tsv
│   ├── adt_gmm_metrics.tsv
│   └── adt_density_plots/
├── 06_prior_review/
│   ├── marker_prior_review.tsv
│   └── marker_prior_review.llm.jsonl
├── 07_final_triage/
│   ├── final_adt_triage_table.tsv
│   ├── final_adt_triage_table.xlsx
│   └── manual_review_table.tsv
└── report/
    ├── ADTriage_report.html
    └── ADTriage_report.pdf
```

---

## 5. 核心输出表

### 5.1 adt_feature_mapping.initial.tsv

该表由 `FB_ref.csv`、`mAOC_gene_symbol_map.csv` 和规则映射生成。

建议字段：

```text
feature_id
feature_name
feature_type
genesymbol_from_map
genesymbol_source
last_seen_name
last_seen_reference
target_class
human_gene_symbol
external_target_name
validation_assay
validation_feature
mapping_method
mapping_confidence
mapping_evidence
needs_llm_review
needs_manual_review
```

---

### 5.2 adt_feature_mapping.review_ready.tsv

该表在大模型辅助判断后生成，但仍不是最终结果。

建议字段：

```text
feature_id
feature_name
target_class
human_gene_symbol
external_target_name
validation_assay
validation_feature
mapping_method
mapping_confidence
llm_suggested_target_class
llm_suggested_gene_symbol
llm_confidence
llm_reason
needs_manual_review
manual_curated_target_class
manual_curated_gene_symbol
manual_notes
```

---

### 5.3 rna_adt_correlation_summary.tsv

针对普通 gene-level ADT feature，计算 ADT CLR 信号与 RNA log-normalized 信号在不同分组变量下的一致性。

每个散点代表一个分组水平，例如一个 cluster、一个 cell type 或一个 sample。

建议字段：

```text
feature_id
feature_name
human_gene_symbol
group_by
n_groups
min_cells_per_group
pearson_r
pearson_p
spearman_r
spearman_p
linear_model_slope
linear_model_intercept
linear_model_r2
correlation_status
```

推荐规则：

```text
n_groups < 5:
  correlation_status = insufficient_groups

spearman_r >= 0.6:
  correlation_status = concordant

0.3 <= spearman_r < 0.6:
  correlation_status = weakly_concordant

spearman_r < 0.3:
  correlation_status = discordant
```

阈值可以在配置文件中调整。

---

### 5.4 cd45isoform_adt_correlation_summary.tsv

针对 CD45 isoform ADT feature，计算 ADT CLR 信号与 `CD45isoforms` assay 中对应 feature 的一致性。

建议字段：

```text
feature_id
feature_name
cd45_isoform
validation_feature
group_by
n_groups
pearson_r
pearson_p
spearman_r
spearman_p
linear_model_slope
linear_model_r2
correlation_status
```

其中：

```text
CD45RA -> PTPRC-RA
CD45RB -> PTPRC-RB
CD45RC -> PTPRC-RC
CD45RO -> PTPRC-RO
```

---

### 5.5 adt_distribution_metrics.tsv

该表记录每个 ADT feature 的 CLR 分布指标。

建议字段：

```text
feature_id
feature_name
target_class
n_cells
mean_clr
median_clr
sd_clr
q05_clr
q25_clr
q75_clr
q95_clr
dynamic_range_clr
positive_rate
negative_component_mean
positive_component_mean
negative_component_sd
positive_component_sd
gmm_n_components
gmm_delta_mu
gmm_separation_index
gmm_overlap_area
snr_raw
tnr_raw
pos_cv
distribution_status
```

推荐解释：

```text
dynamic_range_clr:
  ADT CLR 信号的动态范围

positive_rate:
  根据阈值判定的阳性细胞比例

gmm_delta_mu:
  阳性成分均值 - 阴性成分均值

gmm_separation_index:
  双峰分离程度

snr_raw:
  raw ADT count 层面的 signal-to-noise ratio

tnr_raw:
  true negative ratio 或 background ratio，具体定义需在配置文件中固定
```

---

### 5.6 marker_prior_review.tsv

该表用于记录每个 ADT feature 根据先验知识预期表达的细胞群。

注意：
ADTriage 不负责产生细胞注释。
该表中的细胞类型或分组名称来自用户指定的 Seurat metadata 分组列。

建议字段：

```text
feature_id
feature_name
target_class
human_gene_symbol
group_by
observed_positive_groups
expected_positive_groups
expected_negative_groups
prior_confidence
prior_source
prior_reason
prior_consistency_score
needs_manual_review
```

如果大模型无法以 `confidence >= 0.8` 判断预期表达群，则：

```text
needs_manual_review = TRUE
```

---

### 5.7 final_adt_triage_table.tsv

最终汇总表。

建议字段：

```text
feature_id
feature_name
feature_type
target_class
human_gene_symbol
external_target_name
validation_assay
validation_feature
mapping_confidence
best_group_by
rna_adt_spearman_r
cd45isoform_adt_spearman_r
gmm_separation_index
snr_raw
dynamic_range_clr
positive_rate
prior_consistency_score
final_score
final_label
final_confidence
review_reason
needs_manual_review
manual_final_label
manual_notes
```

---

## 6. 主要分析流程

### Step 0. 输入检查

检查内容：

```text
1. Seurat 对象是否可读取
2. RNA assay 是否存在
3. ADT assay 是否存在
4. CD45isoforms assay 是否存在
5. CD45isoforms assay 是否包含 PTPRC-RA、PTPRC-RB、PTPRC-RC、PTPRC-RO
6. 用户指定的 group_by metadata 列是否存在
7. FB_ref.csv 是否包含必要列
8. mAOC_gene_symbol_map.csv 是否包含必要列
9. FB_ref.csv$id 是否能与 ADT assay feature 对应
10. mAOC_gene_symbol_map.csv$id 是否能与 FB_ref.csv$id 对应
```

失败条件：

```text
RNA assay 缺失
ADT assay 缺失
CD45isoforms assay 缺失
用户指定 group_by 列全部缺失
FB_ref.csv 缺少 id 或 name 列
mAOC_gene_symbol_map.csv 缺少 id 或 genesymbol 列
```

---

### Step 1. 构建 ADT feature 映射表

输入：

```text
FB_ref.csv
mAOC_gene_symbol_map.csv
规则字典
可选：大模型辅助判断
```

输出：

```text
01_feature_mapping/adt_feature_mapping.initial.tsv
01_feature_mapping/adt_feature_mapping.review_ready.tsv
```

核心逻辑：

```text
1. 以 FB_ref.csv 为主表
2. 使用 id 左连接 mAOC_gene_symbol_map.csv
3. 对 CD45 / CD45RA / CD45RB / CD45RC / CD45RO 应用规则
4. 对 spike、IgFc、isotype、control 应用对照规则
5. 对无法判断的条目调用大模型
6. 所有 confidence < 0.8 的条目进入人工审查
```

---

### Step 2. 绘制 ADT 与验证特征图

对于每个 ADT feature，根据 `validation_assay` 绘制配对图。

普通 gene-level ADT：

```text
ADT feature CLR signal
RNA gene log-normalized signal
```

CD45 isoform ADT：

```text
ADT feature CLR signal
CD45isoforms assay 中对应 isoform signal
```

control / spike feature：

```text
仅绘制 ADT feature signal
不绘制 RNA 验证特征
```

每个 `group_by` 列单独输出一套图。

示例：

```text
02_feature_plots/group_by_seurat_clusters/FADT0038_CD8A.pdf
02_feature_plots/group_by_manual_annotation/FADT0038_CD8A.pdf
02_feature_plots/group_by_sample/FADT0038_CD8A.pdf
```

---

### Step 3. 按分组变量计算 pseudobulk summary

对每个用户指定的 `group_by` 列，计算每个分组水平中的平均信号。

示例：

```text
group_by = seurat_clusters
group_level = 0, 1, 2, 3, ...
```

普通 ADT：

```text
mean_ADT_CLR
mean_RNA_lognormalized
n_cells
```

CD45 isoform ADT：

```text
mean_ADT_CLR
mean_CD45isoform_signal
n_cells
```

输出：

```text
03_group_summary/adt_group_summary.tsv
03_group_summary/rna_group_summary.tsv
03_group_summary/cd45isoform_group_summary.tsv
```

---

### Step 4. 计算信号一致性

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

推荐统计量：

```text
Pearson correlation
Spearman correlation
linear model slope
linear model R²
n_groups
min_cells_per_group
```

散点图要求：

```text
1. 每个点代表一个 group level
2. x 轴为验证特征信号
3. y 轴为 ADT CLR 信号
4. 使用最小二乘法拟合直线
5. 标注 Spearman r、Pearson r、R²、n_groups
```

输出：

```text
04_signal_consistency/rna_adt_correlation_summary.tsv
04_signal_consistency/cd45isoform_adt_correlation_summary.tsv
04_signal_consistency/consistency_scatter_plots/
```

---

### Step 5. 先验表达审查

对于每个非 control ADT feature，ADTriage 根据以下信息判断其理论预期表达群：

```text
feature_name
human_gene_symbol
target_class
用户指定 group_by 列中的 group levels
观察到的 ADT 阳性 group
观察到的 RNA 或 CD45isoform 高表达 group
```

优先级：

```text
1. 本地 marker prior 配置文件
2. 项目内已有人工审查记录
3. 大模型辅助判断
4. 人工审查
```

大模型只处理缺失或不确定项。
当大模型置信度 `< 0.8` 时，必须进入人工审查。

输出：

```text
05_prior_review/marker_prior_review.tsv
05_prior_review/marker_prior_review.llm.jsonl
```

---

### Step 6. 评估 ADT CLR 分布指标

对每个 ADT feature 计算分布指标。

基础指标：

```text
mean
median
standard deviation
quantile 5%
quantile 25%
quantile 75%
quantile 95%
dynamic range
positive rate
```

模型指标：

```text
Gaussian mixture model components
negative component mean
positive component mean
component separation
overlap area
SNR
TNR
positive CV
```

输出：

```text
06_distribution_metrics/adt_distribution_metrics.tsv
06_distribution_metrics/adt_gmm_metrics.tsv
06_distribution_metrics/adt_density_plots/
```

---

### Step 7. 最终 triage

最终标签由规则门控和加权评分共同决定。

强制规则示例：

```text
1. target_class = control 或 spike_control
   -> final_label = control

2. mapping_confidence < 0.8
   -> final_label = manual_review

3. validation_assay = RNA 且 RNA feature 不存在
   -> final_label = manual_review

4. validation_assay = CD45isoforms 且 isoform feature 不存在
   -> final_label = manual_review

5. ADT 分布无可检测阳性峰，positive_rate 低于阈值
   -> final_label = no_signal

6. ADT 阳性 group 与先验预期 group 冲突
   -> final_label = suspicious 或 wrong

7. RNA–ADT correlation 低，但 ADT 属于 CD45 isoform
   -> 不作为单独负面证据
```

推荐标签定义：

```text
correct:
  映射明确，ADT 信号存在，分布可分离，观察阳性群与先验一致。

suspicious:
  映射基本明确，但至少一个核心证据不一致，例如 RNA–ADT 相关性低、分布分离弱或阳性群偏离预期。

wrong:
  映射明确但观察信号与预期方向冲突，或疑似 feature 名称、id、抗体 target 错配。

no_signal:
  ADT feature 无可检测阳性信号，或阳性比例低于设定阈值。

control:
  spike、IgFc、isotype、negative control 或其他实验对照。

manual_review:
  映射不确定、证据不足、大模型置信度不足或规则冲突。
```

输出：

```text
07_final_triage/final_adt_triage_table.tsv
07_final_triage/manual_review_table.tsv
```

---

## 7. 配置文件示例

推荐使用 `config/ADTriage.yaml` 管理参数。

```yaml
project:
  name: "ADTriage"
  description: "ADT feature-level triage workflow for CITE-seq quality assessment"
  species: "human"

input:
  seurat_rds: "data/qc_seurat.rds"
  fb_ref_csv: "data/FB_ref.csv"
  maoc_gene_symbol_map_csv: "data/mAOC_gene_symbol_map.csv"

assays:
  rna_assay: "RNA"
  adt_assay: "ADT"
  cd45_isoform_assay: "CD45isoforms"

cd45_isoforms:
  CD45RA: "PTPRC-RA"
  CD45RB: "PTPRC-RB"
  CD45RC: "PTPRC-RC"
  CD45RO: "PTPRC-RO"

group_by:
  - "seurat_clusters"
  - "manual_annotation"
  - "sample"

mapping:
  confidence_cutoff_manual_review: 0.8
  use_maoc_gene_symbol_map: true
  use_rule_based_mapping: true
  use_llm_for_unresolved_only: true

correlation:
  min_cells_per_group: 20
  min_groups_for_correlation: 5
  methods:
    - "pearson"
    - "spearman"

distribution:
  use_gmm: true
  max_gmm_components: 3
  min_positive_rate: 0.005
  min_gmm_separation_index: 1.0

llm:
  enabled: true
  confidence_cutoff: 0.8
  cache_enabled: true
  cache_dir: "ADTriage_output/cache/llm"

output:
  outdir: "ADTriage_output"
  overwrite: false
```

---

## 8. 推荐运行方式

推荐入口：

```bash
Rscript scripts/run_ADTriage.R \
  --config config/ADTriage.yaml
```

推荐分步运行：

```bash
Rscript scripts/00_validate_input.R \
  --config config/ADTriage.yaml

Rscript scripts/01_build_feature_mapping.R \
  --config config/ADTriage.yaml

Rscript scripts/02_plot_features.R \
  --config config/ADTriage.yaml

Rscript scripts/03_group_summary.R \
  --config config/ADTriage.yaml

Rscript scripts/04_signal_consistency.R \
  --config config/ADTriage.yaml

Rscript scripts/05_distribution_metrics.R \
  --config config/ADTriage.yaml

Rscript scripts/06_prior_review.R \
  --config config/ADTriage.yaml

Rscript scripts/07_final_triage.R \
  --config config/ADTriage.yaml
```

---

## 9. 软件依赖

### 9.1 R 依赖

推荐 R 版本：

```text
R >= 4.3.0
```

核心 R 包：

```r
Seurat
SeuratObject
Matrix
data.table
dplyr
tidyr
stringr
ggplot2
patchwork
cowplot
mclust
jsonlite
yaml
readr
openxlsx
```

可选 R 包：

```r
mixtools
flexmix
pROC
rmarkdown
knitr
```

---

### 9.2 Python 依赖

如果启用大模型辅助判断或外部 gene symbol 查询，可配置 Python 环境。

推荐 Python 版本：

```text
Python >= 3.10
```

核心 Python 包：

```text
pandas
numpy
pyyaml
jsonschema
requests
```

---

## 10. 大模型调用原则

ADTriage 中的大模型只用于低频、结构化、可审查的判断任务。

允许大模型处理：

```text
1. 无法通过规则识别的 ADT feature 名称
2. 无法确定的 gene symbol 映射
3. marker prior 缺失时的预期表达群判断
4. manual_review 行的解释性总结
```

禁止大模型处理：

```text
1. 原始表达矩阵
2. 完整 Seurat 对象
3. 逐细胞 metadata
4. 大规模 FeaturePlot 图像
5. 可由 R/Python 确定性计算完成的统计任务
```

大模型输出必须是 JSON 或 JSONL，并保留：

```text
model_name
prompt_version
input_hash
confidence
reason
created_at
```

推荐缓存策略：

```text
cache_key = sha256(feature_id + feature_name + last_seen_name + prompt_version)
```

---

## 11. 人工审查机制

以下情况必须进入人工审查：

```text
1. mapping_confidence < 0.8
2. llm_confidence < 0.8
3. human_gene_symbol 为空且 target_class 不是 control
4. validation_assay 指向的 feature 不存在
5. RNA–ADT 或 CD45isoform–ADT 结果与先验表达冲突
6. GMM 无法稳定拟合
7. final_label = suspicious
8. final_label = wrong
9. final_label = manual_review
```

人工审查表：

```text
07_final_triage/manual_review_table.tsv
```

建议人工审查后填写：

```text
manual_final_label
manual_curated_gene_symbol
manual_curated_target_class
manual_notes
reviewer
reviewed_at
```

---

## 12. 推荐项目结构

```text
ADTriage/
├── README.md
├── config/
│   ├── ADTriage.yaml
│   ├── marker_prior.yaml
│   ├── cd45_isoform_rules.yaml
│   └── scoring_rules.yaml
├── scripts/
│   ├── run_ADTriage.R
│   ├── 00_validate_input.R
│   ├── 01_build_feature_mapping.R
│   ├── 02_plot_features.R
│   ├── 03_group_summary.R
│   ├── 04_signal_consistency.R
│   ├── 05_distribution_metrics.R
│   ├── 06_prior_review.R
│   └── 07_final_triage.R
├── R/
│   ├── io.R
│   ├── validation.R
│   ├── mapping.R
│   ├── plotting.R
│   ├── correlation.R
│   ├── distribution_metrics.R
│   ├── scoring.R
│   └── report.R
├── python/
│   ├── llm_mapping_review.py
│   ├── llm_marker_prior_review.py
│   └── validate_llm_json.py
├── prompts/
│   ├── adt_mapping_review.md
│   └── marker_prior_review.md
├── schemas/
│   ├── adt_mapping_review.schema.json
│   ├── marker_prior_review.schema.json
│   └── final_triage.schema.json
├── tests/
│   ├── test_feature_mapping.R
│   ├── test_cd45isoform_rules.R
│   └── test_distribution_metrics.R
└── example/
    ├── FB_ref.example.csv
    ├── mAOC_gene_symbol_map.example.csv
    └── ADTriage.example.yaml
```

---

## 13. 设计原则

ADTriage 遵循以下原则：

```text
1. 以 feature_id 为主键，保证 FB_ref.csv、mAOC_gene_symbol_map.csv、Seurat ADT feature 可追溯。
2. 每一步生成新表，不直接覆盖上一阶段结果。
3. 统计计算由 R/Python 完成，不交给大模型。
4. 大模型只处理 unresolved / ambiguous / prior-missing 项。
5. CD45 isoform 与 PTPRC gene-level 表达分开评估。
6. correlation 是诊断证据，不是单独判定标准。
7. 最终标签必须保留 review_reason。
8. 所有低置信度结果进入人工审查。
```

---

## 14. 最小可用输出

ADTriage 的 MVP 至少应生成以下 5 个文件：

```text
01_feature_mapping/adt_feature_mapping.review_ready.tsv
04_signal_consistency/rna_adt_correlation_summary.tsv
04_signal_consistency/cd45isoform_adt_correlation_summary.tsv
05_distribution_metrics/adt_distribution_metrics.tsv
07_final_triage/final_adt_triage_table.tsv
```

其中最终主表为：

```text
07_final_triage/final_adt_triage_table.tsv
```

该表是 ADTriage 的核心结果，用于后续人工审查、报告生成和 ADT panel 质量评估。
