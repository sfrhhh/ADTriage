# ADTriage

ADTriage 是一个面向 CITE-seq / Feature Barcode 数据的 **ADT feature-level quality triage workflow**。

它的目标不是重新做单细胞上游分析、细胞注释或整合分析，而是对每一个 ADT feature 汇总 mapping、RNA/CD45isoform 验证、group-level 信号一致性、marker prior、ADT 分布指标和 target enrichment，最终生成可追溯的 ADT triage 标签。

当前主结果：

```text
ADTriage_output/07_final_triage/final_adt_triage_table.llm_resolved.tsv
```

当前最终标签：

```text
correct
suspicious
wrong
no_signal
control
```

---

## 1. 项目定位

ADTriage 回答的问题是：

```text
一个 ADT feature 的 target mapping 是否明确？
ADT 信号是否存在且可分离？
ADT 阳性群体是否符合 marker prior 或 RNA/CD45isoform 验证方向？
该 feature 最终应作为 correct、suspicious、wrong、no_signal 还是 control？
```

ADTriage 不负责：

```text
1. 运行 standard_scop
2. 运行 CellTypist
3. 生成细胞注释
4. 生成 Seurat 对象主分析结果
5. 做整合分析
6. 给出最终生物学结论
```

输入 Seurat 对象应已经包含可用的 assay、reduction 和 metadata。ADTriage 只在此基础上做 feature-level triage。

---

## 2. 当前项目结构

```text
ADTriage/
├── AGENTS.md
├── README.md
├── scripts/
│   ├── build_adt_feature_mapping.py
│   ├── apply_llm_feature_mapping_review.py
│   ├── plot_adt_validation_features.R
│   ├── build_adt_group_summary.R
│   ├── calculate_signal_consistency.R
│   ├── build_marker_prior_review.py
│   ├── apply_marker_prior_llm_review.py
│   ├── build_adt_distribution_metrics.R
│   └── build_final_adt_triage.py
├── priors/
│   └── adult_cHSPC_marker_prior.strict.tsv
├── operation_notes/
│   ├── step1_feature_mapping_reproducible_workflow.md
│   ├── step2_adt_validation_feature_plots_reproducible_workflow.md
│   ├── step3_group_summary_reproducible_workflow.md
│   ├── step4_signal_consistency_reproducible_workflow.md
│   ├── step5_prior_review_reproducible_workflow.md
│   ├── step6_distribution_metrics_reproducible_workflow.md
│   └── step7_final_triage_reproducible_workflow.md
├── tests/data/
│   ├── test_FB_ref.csv
│   ├── mAOC_gene_symbol_map.csv
│   └── test_seurat_qc.rds
└── ADTriage_output/
    ├── 01_feature_mapping/
    ├── 02_feature_plots/
    ├── 03_group_summary/
    ├── 04_signal_consistency/
    ├── 05_prior_review/
    ├── 06_distribution_metrics/
    └── 07_final_triage/
```

每个 Step 的完整复现细节以 `operation_notes/` 中的 Markdown 为准。

---

## 3. 输入数据

### 3.1 Feature reference

本项目测试输入：

```text
tests/data/test_FB_ref.csv
```

核心字段：

```text
id
name
read
pattern
sequence
feature_type
```

`id` 是 Step 1 mapping 的主关联键。

### 3.2 ADT id 到 gene symbol 映射表

本项目测试输入：

```text
tests/data/mAOC_gene_symbol_map.csv
```

核心字段：

```text
id
genesymbol
genesymbol_source
notes
last_seen_name
last_seen_reference
updated_at
```

注意：该文件可能带 BOM，Python 读取时使用 `utf-8-sig`。

### 3.3 Seurat RDS

本项目测试输入：

```text
tests/data/test_seurat_qc.rds
```

当前 workflow 依赖：

```text
assay:
  RNA
  ADT
  CD45isoforms

metadata:
  Adult_cHSPCs_Ultima_majority_voting
  Standardpcaclusters

reduction:
  standardpcaumap2d
```

`CD45isoforms` assay 至少需要：

```text
PTPRC-RA
PTPRC-RB
PTPRC-RC
PTPRC-RO
```

---

## 4. 环境

Python 脚本使用 base conda 环境：

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B <script.py>
```

R 脚本使用 `scop` conda 环境：

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript <script.R>
```

Step 2 需要将 PDF/图像拼接结果转为 PNG 时使用：

```text
/data/project/yanhuixu/frsui/software/bin/magick
```

---

## 5. Workflow 总览

当前 workflow 共 7 步：

```text
Step 1  ADT feature mapping
Step 2  ADT 与 RNA/CD45isoform validation feature plot
Step 3  group-level pseudobulk summary
Step 4  ADT 与 validation feature 的 signal consistency
Step 5  marker prior review
Step 6  ADT distribution metrics
Step 7  final ADT triage
```

---

## 6. Step 1: ADT Feature Mapping

脚本：

```text
scripts/build_adt_feature_mapping.py
scripts/apply_llm_feature_mapping_review.py
```

输入：

```text
tests/data/test_FB_ref.csv
tests/data/mAOC_gene_symbol_map.csv
```

输出：

```text
ADTriage_output/01_feature_mapping/adt_feature_mapping.initial.tsv
ADTriage_output/01_feature_mapping/adt_feature_mapping.review_ready.tsv
ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv
```

核心逻辑：

```text
1. 以 FB_ref.csv 为主表。
2. 使用 id 左连接 mAOC_gene_symbol_map.csv。
3. 对 CD45 / CD45RA / CD45RB / CD45RC / CD45RO 应用规则。
4. 对 spike、IgFc、isotype、control、HTO 等应用技术特征规则。
5. 对 unresolved、ambiguous、多基因 validation feature 使用 LLM/curation。
6. 输出可用于后续绘图和验证的 validation_assay / validation_feature。
```

特殊规则：

```text
CD45RA -> CD45isoforms:PTPRC-RA
CD45RB -> CD45isoforms:PTPRC-RB
CD45RC -> CD45isoforms:PTPRC-RC
CD45RO -> CD45isoforms:PTPRC-RO
CD45   -> RNA:PTPRC
HTO-*  -> hashtag / HTO assay context
```

多基因例子：

```text
mAOC-CD16    初始 FCGR3A/FCGR3B，review 后 validation_feature = FCGR3A
mAOC-HLA-ABC 初始 HLA-A/HLA-B/HLA-C，review 后 validation_feature = HLA-A
```

复现说明：

```text
operation_notes/step1_feature_mapping_reproducible_workflow.md
```

---

## 7. Step 2: ADT Validation Feature Plots

脚本：

```text
scripts/plot_adt_validation_features.R
```

输入：

```text
tests/data/test_seurat_qc.rds
ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv
```

输出：

```text
ADTriage_output/02_feature_plots/adt_validation_plot_feature_top175_page001.png
...
ADTriage_output/02_feature_plots/adt_validation_plot_feature_top175_page008.png
ADTriage_output/02_feature_plots/adt_validation_plot_layout_index_top175.tsv
ADTriage_output/02_feature_plots/adt_validation_plot_status_top175.tsv
```

核心逻辑：

```text
1. 每个 ADT feature 绘制一对 panel：ADT panel + validation panel。
2. 普通 gene-level ADT 使用 RNA validation_feature。
3. CD45 isoform 使用 CD45isoforms validation_feature。
4. control / spike / hashtag 不做额外特殊处理；若绘图失败，则在原位置生成等尺寸占位图。
5. 所有单图先按固定 panel 尺寸输出，再拼接成大图。
6. 图标题包含序号，layout index TSV 用于在大图中定位 feature。
```

当前关键参数：

```text
reduction = standardpcaumap2d
ADT NormalizeData(..., normalization.method = "CLR", margin = 2)
CD45isoforms NormalizeData(..., normalization.method = "CLR")
HTO NormalizeData(..., assay = "HTO", normalization.method = "CLR")
panel_px = 800
只绘制 plot_feature，不再绘制 density feature plot
```

复现说明：

```text
operation_notes/step2_adt_validation_feature_plots_reproducible_workflow.md
```

---

## 8. Step 3: Group-Level Pseudobulk Summary

脚本：

```text
scripts/build_adt_group_summary.R
```

输入：

```text
tests/data/test_seurat_qc.rds
ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv
```

输出：

```text
ADTriage_output/03_group_summary/adt_group_summary.tsv
ADTriage_output/03_group_summary/rna_group_summary.tsv
ADTriage_output/03_group_summary/cd45isoform_group_summary.tsv
```

当前 group_by：

```text
Adult_cHSPCs_Ultima_majority_voting
Standardpcaclusters
```

核心逻辑：

```text
1. 对每个 group_by 的每个 group level 计算 mean_ADT_CLR。
2. 对 RNA validation feature 计算 mean_RNA_lognormalized。
3. 对 CD45isoforms validation feature 计算 mean_CD45isoform_signal。
4. 同时记录 n_cells、status 和 message。
```

复现说明：

```text
operation_notes/step3_group_summary_reproducible_workflow.md
```

---

## 9. Step 4: Signal Consistency

脚本：

```text
scripts/calculate_signal_consistency.R
```

输入：

```text
ADTriage_output/03_group_summary/rna_group_summary.tsv
ADTriage_output/03_group_summary/cd45isoform_group_summary.tsv
```

输出：

```text
ADTriage_output/04_signal_consistency/rna_adt_correlation_summary.tsv
ADTriage_output/04_signal_consistency/cd45isoform_adt_correlation_summary.tsv
ADTriage_output/04_signal_consistency/rna_adt_<group_by>_consistency_scatter_plots.pdf
ADTriage_output/04_signal_consistency/cd45isoform_adt_<group_by>_consistency_scatter_plots.pdf
```

统计量：

```text
pearson_correlation
spearman_correlation
linear_model_slope
linear_model_r_squared
n_groups
min_cells_per_group
```

散点图：

```text
x = validation feature group mean
y = ADT CLR group mean
point = group level
line = least-squares fit
```

复现说明：

```text
operation_notes/step4_signal_consistency_reproducible_workflow.md
```

---

## 10. Step 5: Marker Prior Review

脚本：

```text
scripts/build_marker_prior_review.py
scripts/apply_marker_prior_llm_review.py
```

输入：

```text
ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv
ADTriage_output/03_group_summary/adt_group_summary.tsv
ADTriage_output/03_group_summary/rna_group_summary.tsv
ADTriage_output/03_group_summary/cd45isoform_group_summary.tsv
priors/adult_cHSPC_marker_prior.strict.tsv
```

输出：

```text
ADTriage_output/05_prior_review/marker_prior_review.tsv
ADTriage_output/05_prior_review/marker_prior_review.llm.jsonl
ADTriage_output/05_prior_review/marker_prior_review.llm_responses.jsonl
ADTriage_output/05_prior_review/marker_prior_review.llm_resolved.tsv
```

当前只对生物学 metadata 运行：

```text
Adult_cHSPCs_Ultima_majority_voting
```

`Standardpcaclusters` 是数字 cluster，不在 Step 5 做 marker prior review。

strict prior 原则：

```text
1. 只有 lineage-defining 或 lineage-informative marker 才给 expected_positive_groups。
2. broad marker、activation marker、adhesion marker、checkpoint marker 不做 hard vote。
3. pDC/DC/basophil/mast/endothelial context marker 不强行归入现有 HSPC label。
4. mature monocyte marker 不强行支持 GMP-L。
```

先验文件：

```text
priors/adult_cHSPC_marker_prior.strict.tsv
```

后续更新先验时优先改这个 TSV，不需要改 Python 脚本。

核心输出字段：

```text
expected_positive_groups
observed_adt_positive_groups
observed_validation_high_groups
adt_validation_overlap_groups
expected_adt_overlap_groups
expected_validation_overlap_groups
prior_method
prior_confidence
needs_llm_review
needs_manual_review
```

复现说明：

```text
operation_notes/step5_prior_review_reproducible_workflow.md
```

---

## 11. Step 6: ADT Distribution Metrics

脚本：

```text
scripts/build_adt_distribution_metrics.R
```

输入：

```text
tests/data/test_seurat_qc.rds
ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv
ADTriage_output/05_prior_review/marker_prior_review.llm_resolved.tsv
```

输出：

```text
ADTriage_output/06_distribution_metrics/adt_distribution_metrics.tsv
ADTriage_output/06_distribution_metrics/adt_gmm_metrics.tsv
ADTriage_output/06_distribution_metrics/adt_density_plots.pdf
```

核心逻辑：

```text
1. 对 ADT assay 重新计算 CLR：NormalizeData(..., normalization.method = "CLR", margin = 2)。
2. 默认只用 CLR > 0 的值进行 GMM 拟合。
3. 默认密度图也只展示 CLR > 0，避免 0 值尖峰干扰看图。
4. threshold 应用于全体细胞，因此 positive_rate 包含 CLR = 0 的细胞。
5. TNR 使用 raw count 计算，包含 raw count = 0 的细胞。
6. Step 5 的 expected_positive_groups 用于 target / non-target 定义。
```

GMM/ThresholdR 风格设计：

```text
1. 用 mclust 对 k = 1,2,3 计算 BIC。
2. BIC delta 默认 >= 10 才升级 component 数。
3. 优先使用 mixtools::normalmixEM，失败时回退 mclust。
4. 默认 selected threshold = lowest component mean + 3 * sd。
5. 同时记录 component intersection threshold、overlap、separation 等指标。
```

`adt_signal_class` 是 Step 6 的解释性分布分类，不是最终 triage 标签：

```text
target_enriched
rare_target_enriched
rare_target_enriched_low_count
distribution_detected_but_off_target
no_signal
no_target_enrichment
no_target_detected_signal
target_raw_enriched_no_threshold
...
```

复现说明：

```text
operation_notes/step6_distribution_metrics_reproducible_workflow.md
```

---

## 12. Step 7: Final ADT Triage

脚本：

```text
scripts/build_final_adt_triage.py
```

输入：

```text
ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv
ADTriage_output/04_signal_consistency/rna_adt_correlation_summary.tsv
ADTriage_output/04_signal_consistency/cd45isoform_adt_correlation_summary.tsv
ADTriage_output/05_prior_review/marker_prior_review.llm_resolved.tsv
ADTriage_output/06_distribution_metrics/adt_distribution_metrics.tsv
```

输出：

```text
ADTriage_output/07_final_triage/final_adt_triage_table.tsv
ADTriage_output/07_final_triage/llm_review_queue.tsv
ADTriage_output/07_final_triage/final_triage.llm.jsonl
ADTriage_output/07_final_triage/final_triage.llm_responses.jsonl
ADTriage_output/07_final_triage/final_adt_triage_table.llm_resolved.tsv
```

本地规则只处理确定性场景：

```text
control/spike/HTO -> control
Step6 no_signal/no_target_no_signal -> no_signal
Step6 target_enriched 且无其它冲突 -> correct
mapping/validation 硬阻断 -> llm_review_needed
语义冲突、off-target、ambiguous、raw-only enrichment -> LLM 队列
```

LLM response 合并后生成终版：

```text
final_adt_triage_table.llm_resolved.tsv
```

本次终版标签分布：

```text
suspicious = 75
correct = 54
no_signal = 33
wrong = 7
control = 6
```

代表条目：

```text
mAOC-CD16  -> correct
mAOC-CD16a -> correct
mAOC-CD34  -> correct
mAOC-CD19  -> no_signal
mAOC-CD23  -> wrong
mAOC-CD49f -> wrong
HTO/spike  -> control
```

复现说明：

```text
operation_notes/step7_final_triage_reproducible_workflow.md
```

---

## 13. 一键式复现命令参考

以下命令展示当前项目中各步骤的主要调用方式。完整参数解释见 `operation_notes/`。

### 13.1 Step 1

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/build_adt_feature_mapping.py \
  --fb-ref tests/data/test_FB_ref.csv \
  --gene-map tests/data/mAOC_gene_symbol_map.csv \
  --outdir ADTriage_output/01_feature_mapping

/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/apply_llm_feature_mapping_review.py \
  --initial ADTriage_output/01_feature_mapping/adt_feature_mapping.initial.tsv \
  --output ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv \
  --threshold 0.8
```

### 13.2 Step 2

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript scripts/plot_adt_validation_features.R \
  --seurat-rds tests/data/test_seurat_qc.rds \
  --mapping ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv \
  --outdir ADTriage_output/02_feature_plots \
  --n-features 175 \
  --reduction standardpcaumap2d \
  --panel-px 800
```

### 13.3 Step 3

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript scripts/build_adt_group_summary.R \
  --seurat-rds tests/data/test_seurat_qc.rds \
  --mapping ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv \
  --outdir ADTriage_output/03_group_summary \
  --group-by Adult_cHSPCs_Ultima_majority_voting,Standardpcaclusters
```

### 13.4 Step 4

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript scripts/calculate_signal_consistency.R \
  --group-summary-dir ADTriage_output/03_group_summary \
  --outdir ADTriage_output/04_signal_consistency \
  --min-groups 3
```

### 13.5 Step 5

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/build_marker_prior_review.py \
  --mapping ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv \
  --group-summary-dir ADTriage_output/03_group_summary \
  --outdir ADTriage_output/05_prior_review \
  --prior-file priors/adult_cHSPC_marker_prior.strict.tsv \
  --group-by Adult_cHSPCs_Ultima_majority_voting \
  --top-quantile 0.75

/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/apply_marker_prior_llm_review.py \
  --review-table ADTriage_output/05_prior_review/marker_prior_review.tsv \
  --llm-responses ADTriage_output/05_prior_review/marker_prior_review.llm_responses.jsonl \
  --out ADTriage_output/05_prior_review/marker_prior_review.llm_resolved.tsv \
  --confidence-threshold 0.8
```

### 13.6 Step 6

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript scripts/build_adt_distribution_metrics.R \
  --seurat-rds tests/data/test_seurat_qc.rds \
  --mapping ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv \
  --prior-review ADTriage_output/05_prior_review/marker_prior_review.llm_resolved.tsv \
  --outdir ADTriage_output/06_distribution_metrics \
  --group-by Adult_cHSPCs_Ultima_majority_voting
```

### 13.7 Step 7

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/build_final_adt_triage.py \
  --outdir ADTriage_output/07_final_triage \
  --group-by Adult_cHSPCs_Ultima_majority_voting

/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/build_final_adt_triage.py \
  --outdir ADTriage_output/07_final_triage \
  --group-by Adult_cHSPCs_Ultima_majority_voting \
  --llm-responses ADTriage_output/07_final_triage/final_triage.llm_responses.jsonl
```

---

## 14. 输出解读

建议下游默认读取：

```text
ADTriage_output/07_final_triage/final_adt_triage_table.llm_resolved.tsv
```

常用字段：

```text
feature_id
feature_name
target_class
human_gene_symbol
validation_assay
validation_feature
mapping_confidence
expected_positive_groups
observed_adt_positive_groups
observed_validation_high_groups
expected_adt_overlap_groups
expected_validation_overlap_groups
spearman_correlation
positive_rate
positive_n_cells
target_positive_rate
non_target_positive_rate
tnr_raw
gmm_quality_status
adt_signal_class
final_label
final_confidence
decision_source
final_reason
recommended_action
```

标签解释：

```text
correct:
  当前 workflow 中较可信的 ADT feature。

suspicious:
  不建议作为强 target marker；可保留用于探索或后续单独复核。

wrong:
  不建议用于对应 target 的生物学解释，疑似 off-target、失败或错配。

no_signal:
  无可用 ADT 信号，不建议用于生物学解释。

control:
  技术对照，不进入 biological target correctness 排名。
```

---

## 15. 关键注意事项

```text
1. Step 6 的 adt_signal_class 不是最终标签，最终标签以 Step 7 的 final_label 为准。
2. Step 5 的 strict prior 是当前 adult cHSPC label 集合下的保守先验，不代表 marker 在所有生物学语境中都没有表达。
3. no-vote marker 不等于没有表达，只表示不作为当前 label 集合中的强先验证据。
4. CD45 isoform 必须用 CD45isoforms assay 评估，不能只用 RNA:PTPRC。
5. TNR 中的 target 由 Step 5 expected_positive_groups 定义。
6. Step 6 拟合和密度图默认只使用 CLR > 0，但 positive_rate 和 TNR 仍基于全体细胞。
7. LLM 只处理语义裁决和冲突解释；所有输入证据和输出 response 都保留为 JSONL/TSV。
```

---

## 16. 当前完成状态

当前项目已经完成并提交了 Step 1-7：

```text
Step 1  feature mapping
Step 2  validation feature plot
Step 3  group summary
Step 4  signal consistency
Step 5  marker prior review
Step 6  distribution metrics
Step 7  final triage
```

每一步均包含：

```text
1. 可复跑脚本
2. 输出目录
3. operation_notes 可复现说明
```
