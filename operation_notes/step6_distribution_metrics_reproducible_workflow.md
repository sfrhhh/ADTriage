# Step 6: ADT Distribution Metrics 可复现操作说明

## 1. 目标

Step 6 的目标是对每个 ADT feature 计算全局分布指标、候选阈值、GMM 拟合指标、target / non-target enrichment 指标，并输出每个 feature 的密度图用于人工复核。

本步骤不是单纯依赖 GMM 阈值判断抗体质量，而是整合三层信息：

```text
1. ADT CLR / raw count 的全局分布。
2. 基于 CLR > 0 的混合模型拟合和候选阈值。
3. Step 5 的 expected_positive_groups 与 observed group overlap。
```

这样可以区分：

```text
target_enriched:
  分布阈值和 target enrichment 均支持该 ADT 在预期细胞类型中表达。

rare_target_enriched / rare_target_enriched_low_count:
  阳性细胞很少，但阳性信号确实富集在预期 target 中。

distribution_detected_but_off_target:
  有可检测的 ADT-positive 分布，但不富集在 Step 5 预期 target 中。
  例如本次 mAOC-CD23 在密度图上有弱双峰，但 target=B 中 ADT 不富集。

no_signal:
  低阳性率、低全局 raw / CLR 信号，且没有 target enrichment。
```

## 2. 输入文件

```text
tests/data/test_seurat_qc.rds
ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv
ADTriage_output/05_prior_review/marker_prior_review.llm_resolved.tsv
```

Seurat 对象要求：

```text
assay:
  ADT

metadata:
  Adult_cHSPCs_Ultima_majority_voting
```

Step 5 prior review 表为 TNR 和 target enrichment 提供：

```text
expected_positive_groups
observed_adt_positive_groups
observed_validation_high_groups
adt_validation_overlap_groups
expected_adt_overlap_groups
expected_validation_overlap_groups
needs_manual_review
review_reason
```

## 3. 使用脚本

```text
scripts/build_adt_distribution_metrics.R
```

脚本主要功能：

```text
1. 读取 Seurat RDS。
2. 对 ADT assay 重新计算 CLR：
   NormalizeData(bm, assay = "ADT", normalization.method = "CLR", margin = 2)
3. 读取 Step 1 mapping 表，并只保留 Seurat ADT assay 中存在的 feature。
4. 读取 Step 5 prior review 表，按指定 group_by 取 expected target groups。
5. 对每个 ADT feature 计算全局 CLR / raw UMI 分布指标。
6. 默认仅用 CLR > 0 的细胞进行 GMM 拟合，避免大量 0 值主导拟合。
7. 默认仅绘制 CLR > 0 的密度图，保持拟合对象和图中展示对象一致。
8. 用 threshold 将全体细胞划分为 Positive / Negative，并计算 positive_rate。
9. 用 Step 5 的 expected_positive_groups 计算 target / non-target raw UMI 和 positive rate。
10. 输出 summary TSV、GMM TSV 和统一 PDF。
```

## 4. 关键算法

### 4.1 CLR 与 0 值

ADT 的 `data` 层在脚本中重新计算：

```r
NormalizeData(bm, assay = "ADT", normalization.method = "CLR", margin = 2)
```

默认设置：

```text
GMM fit input:
  CLR > 0

Density plot:
  CLR > 0

positive_rate / positive_n_cells:
  threshold 应用于全体细胞，包括 CLR = 0 的细胞。

target_positive_rate / non_target_positive_rate:
  threshold 应用于全体 target / non-target 细胞，包括 CLR = 0 的细胞。

TNR:
  使用 raw count 均值计算，包含 raw count = 0 的细胞。
```

这个设计的原因是：

```text
1. 拟合和看图时，0 值会形成很强的背景尖峰，干扰非零信号结构判断。
2. positive_rate / TNR 需要反映真实全体细胞中的阳性比例和 dropout，因此保留 0 值。
```

### 4.2 GMM 拟合

本步骤参考 ThresholdR 的思路，但不依赖 `AdaptGauss`。

模型选择：

```text
1. 对 k = 1, 2, 3 分别用 mclust 计算 BIC。
2. 如果 BIC(k=2) - BIC(k=1) >= min_bic_delta，则允许 k=2。
3. 如果 BIC(k=3) - BIC(k=2) >= min_bic_delta，则允许 k=3。
4. 默认 min_bic_delta = 10。
```

拟合：

```text
1. 优先使用 mixtools::normalmixEM。
2. 若 normalmixEM 失败，回退到 mclust::Mclust。
3. 默认每个 feature 最多抽样 3000 个 CLR > 0 细胞进行拟合。
4. 若 CLR > 0 细胞数 < 50 或 unique 值 < 20，则标记 fit_failed。
```

阈值：

```text
threshold_noise_3sd:
  lowest component mean + 3 * lowest component sd。

threshold_intersection_1_2:
  component 1 与 component 2 的加权正态曲线交点。

threshold_intersection_2_3:
  component 2 与 component 3 的加权正态曲线交点，仅 k=3 时存在。

selected_threshold_method:
  默认使用 noise_3sd。
```

低峰 sigma 校正：

```text
若最低 component 在其均值处的模型密度与 empirical KDE 密度差异超过 margin_den，
脚本会优化最低 component 的 sigma，使模型背景峰更贴近 KDE。

默认 margin_den = 0.1。
```

### 4.3 TNR 与 target enrichment

target cells 定义：

```text
meta[[group_by]] %in% expected_positive_groups
```

本次使用：

```text
group_by = Adult_cHSPCs_Ultima_majority_voting
```

TNR 计算：

```text
tnr_raw = mean(raw_count in target cells) / mean(raw_count in non-target cells)
```

脚本使用一个很小的 denominator floor，避免 non-target raw mean 为 0 时出现除零。

target threshold 指标：

```text
target_positive_rate:
  target cells 中 CLR > selected_threshold 的比例。

non_target_positive_rate:
  non-target cells 中 CLR > selected_threshold 的比例。

target_dropout_pct:
  target cells 中 raw count = 0 的比例。
```

## 5. 信号分类逻辑

`adt_signal_class` 是面向下游自动化筛选的解释性分类，不等价于最终生物学结论。

核心判定顺序：

```text
1. 若没有 expected_positive_groups：
   只报告 no_target_* 类型，不强行判断 target enrichment。

2. 若有 expected_positive_groups 且无可用 threshold：
   若 tnr_raw >= 1.5，则标记 target_raw_enriched_no_threshold。
   否则标记 no_target_enrichment。

3. 若有 threshold：
   优先判断 target_enriched：
     tnr_raw >= 1.5
     且满足 target_positive_rate 比 non_target_positive_rate 明显更高。

4. 若 target_enriched 且 positive_n_cells < 20：
   标记 rare_target_enriched_low_count。

5. 若存在 thresholded positive cells，
   但 Step 5 的 expected_adt_overlap_groups 为空，
   或 target_positive_rate <= non_target_positive_rate：
   标记 distribution_detected_but_off_target。

6. 若低阳性率、低 raw UMI、低 CLR dynamic range 且无 target enrichment：
   标记 no_signal。
```

本次关键例子：

```text
mAOC-CD16:
  positive_n_cells = 12
  expected_positive_groups = Monocyte;NK/T
  expected_adt_overlap_groups = Monocyte;NK/T
  tnr_raw = 10.253048
  adt_signal_class = rare_target_enriched_low_count

mAOC-CD23:
  expected_positive_groups = B
  expected_adt_overlap_groups 为空
  target_positive_rate = 0
  non_target_positive_rate = 0.150352
  adt_signal_class = distribution_detected_but_off_target

mAOC-CD34:
  expected_positive_groups = HSC_MPP;MEBEMP-L
  expected_adt_overlap_groups = HSC_MPP
  tnr_raw = 1.773551
  adt_signal_class = target_enriched
```

## 6. 复现命令

在项目根目录运行：

```bash
cd /data/project/yanhuixu/frsui/sc_project/ADTriage
```

### 6.1 R 脚本语法检查

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript -e "invisible(parse('scripts/build_adt_distribution_metrics.R')); cat('parse ok\n')"
```

### 6.2 计算 Step 6 指标并绘图

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript scripts/build_adt_distribution_metrics.R \
  --seurat-rds tests/data/test_seurat_qc.rds \
  --mapping ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv \
  --prior-review ADTriage_output/05_prior_review/marker_prior_review.llm_resolved.tsv \
  --outdir ADTriage_output/06_distribution_metrics \
  --group-by Adult_cHSPCs_Ultima_majority_voting
```

可选参数：

```text
--assay ADT
--max-k 3
--min-bic-delta 10
--fit-sample-size 3000
--min-fit-cells 50
--min-unique 20
--threshold-method noise_3sd
--plot-all-values
--no-positive-only
--margin-den 0.1
--plot-width 8
--plot-height 6
```

说明：

```text
--plot-all-values:
  绘图包含 CLR = 0。

--no-positive-only:
  GMM 拟合包含 CLR = 0。

默认不要开启这两个参数，除非要专门诊断 0 值背景峰。
```

## 7. 输出文件

输出目录：

```text
ADTriage_output/06_distribution_metrics/
```

本次正式输出：

```text
adt_distribution_metrics.tsv
adt_gmm_metrics.tsv
adt_density_plots.pdf
```

### 7.1 adt_distribution_metrics.tsv

每个 ADT feature 一行，包含全局分布指标、target enrichment 指标和最终解释性分类。

关键字段：

```text
feature_id
feature_name
target_class
assay
group_by
expected_positive_groups
prior_method
prior_confidence
observed_adt_positive_groups
observed_validation_high_groups
adt_validation_overlap_groups
expected_adt_overlap_groups
expected_validation_overlap_groups
step5_needs_manual_review
step5_review_reason
n_cells
n_nonzero_raw
zero_rate
nonzero_rate
mean_clr
median_clr
sd_clr
q05_clr
q25_clr
q75_clr
q95_clr
dynamic_range_clr
raw_mean_umi
raw_median_umi
raw_total_umi
selected_threshold_clr
selected_threshold_method
positive_rate
positive_n_cells
negative_n_cells
target_n_cells
non_target_n_cells
target_positive_rate
non_target_positive_rate
snr_raw
tnr_raw
target_dropout_pct
positive_cv_raw
distribution_status
gmm_quality_status
adt_signal_class
adt_signal_reason
```

### 7.2 adt_gmm_metrics.tsv

每个 ADT feature 一行，记录 GMM 拟合和阈值细节。

关键字段：

```text
fit_input
gmm_status
gmm_message
gmm_quality_status
selected_k
fit_method
fit_n_cells
fit_n_cells_full
fit_n_unique
bic_k1
bic_k2
bic_k3
bic_delta_2_vs_1
bic_delta_3_vs_2
sigma_corrected
selected_threshold_clr
threshold_noise_3sd
threshold_mid_3sd
threshold_intersection_1_2
threshold_intersection_2_3
component1_mean
component2_mean
component3_mean
component1_sd
component2_sd
component3_sd
component1_weight
component2_weight
component3_weight
separation_1_2
separation_2_3
overlap_1_2
overlap_2_3
positive_rate
positive_n_cells
target_positive_rate
non_target_positive_rate
tnr_raw
```

### 7.3 adt_density_plots.pdf

统一 PDF，每个 ADT feature 一页。

每页包含：

```text
1. CLR > 0 的 empirical density。
2. GMM component dashed curves。
3. selected threshold vertical line。
4. component intersection threshold vertical line。
5. feature id / feature name / fit status / selected k / quality / positive rate / TNR / signal class / expected target。
```

## 8. 本次运行统计

输出行数：

```text
adt_distribution_metrics.tsv rows = 170 + header
adt_gmm_metrics.tsv rows = 170 + header
adt_density_plots.pdf pages = 170
```

GMM 状态：

```text
ok = 141
fit_failed = 19
unimodal_like = 10
```

GMM quality：

```text
rare_positive = 100
trimodal_candidate = 38
fit_failed = 19
unimodal_like = 10
threshold_uncertain = 2
bimodal_candidate = 1
```

`adt_signal_class` 统计：

```text
no_target_detected_signal = 48
rare_target_enriched = 25
distribution_detected_but_off_target = 24
target_enriched = 23
no_target_enrichment = 14
rare_or_ambiguous = 13
no_target_no_signal = 9
target_raw_enriched_no_threshold = 7
no_target_threshold_unavailable = 2
rare_target_enriched_low_count = 2
ambiguous_trimodal = 1
no_target_broad_signal = 1
target_enriched_broad = 1
```

## 9. 人工检查建议

优先检查以下类别：

```text
distribution_detected_but_off_target:
  有 ADT 分布信号，但不在预期 target 中富集。需要判断是抗体非特异、标注问题、还是真实生物学偏移。

rare_target_enriched_low_count:
  target enrichment 强，但阳性细胞数 < 20。需要结合 feature plot 和原始实验设计判断是否保留。

target_raw_enriched_no_threshold:
  raw UMI 在 target 中富集，但 GMM 没有得到可用 threshold。需要检查低频稀有群体或拟合失败原因。

ambiguous_trimodal:
  k=3 且 target enrichment 不清楚。需要人工看密度图和 feature plot。

threshold_uncertain:
  component overlap 高，阈值不稳定。
```

本次已知需要关注：

```text
mAOC-CD16:
  稀有 target-enriched，当前分类符合“表达细胞很少但方向可信”的判断。

mAOC-CD23:
  密度图有弱双峰，但 ADT 不在 B 细胞中富集，当前分类为 off-target。

mAOC-CD34:
  数据来自 CD34 流式分选，CD34 多数细胞阳性是符合实验设计的；不要按普通低背景抗体理解。
```

