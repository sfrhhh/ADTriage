# Step 7: Final ADT Triage 可复现操作说明

## 1. 目标

Step 7 的目标是整合 Step 1-6 的结构化证据，为每个 ADT feature 生成最终 triage 标签。

本步骤不在本地维护人工审查流程。策略是：

```text
1. 本地规则只处理确定性场景。
2. 需要自然语言判断、证据冲突解释或生物学语义裁决的条目，统一写入 LLM JSONL 队列。
3. LLM 返回结构化 response JSONL 后，由同一个脚本合并回终版 TSV。
```

最终使用的主表是：

```text
ADTriage_output/07_final_triage/final_adt_triage_table.llm_resolved.tsv
```

## 2. 输入文件

```text
ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv
ADTriage_output/04_signal_consistency/rna_adt_correlation_summary.tsv
ADTriage_output/04_signal_consistency/cd45isoform_adt_correlation_summary.tsv
ADTriage_output/05_prior_review/marker_prior_review.llm_resolved.tsv
ADTriage_output/06_distribution_metrics/adt_distribution_metrics.tsv
```

本次使用的生物学分组列：

```text
Adult_cHSPCs_Ultima_majority_voting
```

## 3. 使用脚本

```text
scripts/build_final_adt_triage.py
```

脚本功能：

```text
1. 读取 Step 1 mapping 结果。
2. 读取 Step 4 RNA / CD45isoform group-level correlation 结果。
3. 读取 Step 5 marker prior review 结果。
4. 读取 Step 6 ADT distribution metrics 结果。
5. 按 feature_id 和 group_by 合并证据。
6. 应用确定性规则，生成 final_adt_triage_table.tsv。
7. 对无法确定的条目生成 llm_review_queue.tsv 和 final_triage.llm.jsonl。
8. 如果提供 --llm-responses，则合并 LLM response，生成 final_adt_triage_table.llm_resolved.tsv。
```

## 4. 最终标签定义

```text
correct:
  mapping 明确，ADT 信号存在，并且 ADT 阳性群或 target enrichment 与预期 target 一致。

suspicious:
  有 ADT 信号或部分证据支持，但 mapping、RNA/validation、prior 或分布证据之间存在不一致，尚不足以直接判为 wrong。

wrong:
  mapping/validation 指向某个明确 target，但 ADT-positive groups 与预期 target 方向冲突，或疑似抗体 target / feature id / ADT 信号错配。

no_signal:
  ADT 无可用阳性信号，或对预期 target 没有可解释的 target enrichment。

control:
  spike、HTO、IgG、IgFc、isotype、negative control 或其他技术对照。

llm_review_needed:
  规则初表中的临时标签，表示该 feature 需要 LLM 语义裁决。
  在 llm_resolved 终版表中应尽量不保留。
```

## 5. 本地确定性规则

### 5.1 control

以下条目直接标记为：

```text
final_label = control
decision_source = rule:control
```

规则：

```text
target_class in control / spike_control / hashtag
或 validation_assay = none 且 feature_name 包含 spike / isotype / control / IgG / IgFc / HTO 等技术关键词。
```

### 5.2 mapping / validation 硬阻断

以下情况进入 LLM 队列，而不是本地人工审查：

```text
mapping_confidence < 0.8
Step 1 needs_llm_review = TRUE
Step 1 needs_manual_review = TRUE
validation_assay = RNA 或 CD45isoforms 但 Step 4 没有可用 validation summary
```

### 5.3 no_signal

若 Step 6 的 `adt_signal_class` 为：

```text
no_signal
no_target_no_signal
```

且没有 mapping / validation 硬阻断，则直接：

```text
final_label = no_signal
decision_source = rule:no_signal
```

注意：

```text
RNA-ADT correlation 低不会覆盖明确的 no_signal 结论。
没有 ADT 信号时，RNA correlation 本身不应作为负面或正面 ADT 证据。
```

### 5.4 correct

若 Step 6 的 `adt_signal_class` 为：

```text
target_enriched
rare_target_enriched
rare_target_enriched_low_count
target_enriched_broad
```

且 mapping、Step 5、Step 4 没有其它需要语义裁决的冲突，则：

```text
final_label = correct
decision_source = rule:target_enriched
```

### 5.5 LLM 队列

以下情况默认进入 LLM：

```text
distribution_detected_but_off_target
rare_or_ambiguous
ambiguous_trimodal
target_raw_enriched_no_threshold
no_target_detected_signal
no_target_enrichment
no_target_threshold_unavailable
RNA-ADT Spearman correlation < 0.2
Step 5 needs_manual_review = TRUE
```

这些情况需要结合 marker 生物学、prior、observed group overlap、RNA/ADT 一致性和分布质量进行语义判断。

## 6. LLM 输入格式

LLM 输入文件：

```text
ADTriage_output/07_final_triage/final_triage.llm.jsonl
```

每行是一个 JSON object，核心字段：

```text
feature_id
feature_name
allowed_final_labels
candidate_labels
evidence.mapping
evidence.prior_and_groups
evidence.validation_correlation
evidence.adt_distribution
evidence.rule_queue_reason
response_schema
```

其中 `evidence` 包含：

```text
mapping:
  target_class
  human_gene_symbol
  validation_assay
  validation_feature
  mapping_method
  mapping_confidence

prior_and_groups:
  expected_positive_groups
  observed_adt_positive_groups
  observed_validation_high_groups
  adt_validation_overlap_groups
  expected_adt_overlap_groups
  expected_validation_overlap_groups

validation_correlation:
  status
  spearman_correlation
  pearson_correlation
  linear_model_r_squared
  n_groups

adt_distribution:
  adt_signal_class
  gmm_quality_status
  positive_rate
  positive_n_cells
  target_positive_rate
  non_target_positive_rate
  tnr_raw
  target_dropout_pct
```

## 7. LLM response 格式

LLM response 文件：

```text
ADTriage_output/07_final_triage/final_triage.llm_responses.jsonl
```

每行至少包含：

```json
{
  "feature_id": "mAOC0000",
  "final_label": "correct",
  "final_confidence": 0.9,
  "final_reason": "short rationale",
  "recommended_action": "short downstream action"
}
```

约束：

```text
final_label 必须属于：
  correct
  suspicious
  wrong
  no_signal
  control

final_confidence 必须 >= 0.8 才会被脚本接受。
低于阈值或非法 label 的 response 不会覆盖 final_label = llm_review_needed。
```

本次 LLM 裁决采用保守策略：

```text
1. 明确 target-enriched 才给 correct。
2. 有预期/RNA validation 支持，但 ADT 明显跑到非预期 group 时给 wrong。
3. 无有效先验或 broad/no-vote marker 有信号时给 suspicious。
4. 阈值、raw UMI 和 target enrichment 都不支持可用 ADT 信号时给 no_signal。
5. threshold failed 但 raw TNR 与 observed ADT groups 支持 expected target 时，可给 correct，但 recommended_action 中注明 threshold-unstable。
```

## 8. 复现命令

在项目根目录运行：

```bash
cd /data/project/yanhuixu/frsui/sc_project/ADTriage
```

### 8.1 Python 语法检查

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B -m py_compile \
  scripts/build_final_adt_triage.py
```

### 8.2 生成规则初表和 LLM 队列

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/build_final_adt_triage.py \
  --outdir ADTriage_output/07_final_triage \
  --group-by Adult_cHSPCs_Ultima_majority_voting
```

该命令生成：

```text
final_adt_triage_table.tsv
llm_review_queue.tsv
final_triage.llm.jsonl
```

### 8.3 生成或提供 LLM response

将：

```text
ADTriage_output/07_final_triage/final_triage.llm.jsonl
```

传递给 LLM，得到：

```text
ADTriage_output/07_final_triage/final_triage.llm_responses.jsonl
```

本次 response 已生成并保存在上述路径。

### 8.4 合并 LLM response

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/build_final_adt_triage.py \
  --outdir ADTriage_output/07_final_triage \
  --group-by Adult_cHSPCs_Ultima_majority_voting \
  --llm-responses ADTriage_output/07_final_triage/final_triage.llm_responses.jsonl
```

该命令会重新生成规则初表、LLM 队列，并额外输出：

```text
final_adt_triage_table.llm_resolved.tsv
```

## 9. 输出文件

输出目录：

```text
ADTriage_output/07_final_triage/
```

文件：

```text
final_adt_triage_table.tsv
  规则初表。未进入本地确定性规则的条目保留 final_label = llm_review_needed。

llm_review_queue.tsv
  LLM 队列的表格版，便于搜索和人工快速浏览。

final_triage.llm.jsonl
  LLM 输入 JSONL。

final_triage.llm_responses.jsonl
  LLM 输出 JSONL。

final_adt_triage_table.llm_resolved.tsv
  合并 LLM response 后的终版表。
```

## 10. 终版表字段说明

`final_adt_triage_table.llm_resolved.tsv` 关键字段：

```text
feature_id
feature_name
target_class
human_gene_symbol
validation_assay
validation_feature
mapping_method
mapping_confidence
mapping_needs_llm_review
mapping_needs_manual_review
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
correlation_status
spearman_correlation
pearson_correlation
linear_model_r_squared
n_groups
positive_rate
positive_n_cells
target_positive_rate
non_target_positive_rate
tnr_raw
target_dropout_pct
gmm_quality_status
adt_signal_class
adt_signal_reason
final_label
final_confidence
decision_source
final_reason
recommended_action
llm_review_needed
llm_review_reason
candidate_labels
```

其中：

```text
decision_source = rule:control
  control 确定性规则。

decision_source = rule:no_signal
  Step 6 no-signal 确定性规则。

decision_source = rule:target_enriched
  Step 6 target-enriched 确定性规则。

decision_source = llm
  来自 LLM response 合并。
```

## 11. 本次运行统计

规则初表：

```text
final_adt_triage_table.tsv rows = 175 + header

final_label:
  llm_review_needed = 124
  correct = 36
  no_signal = 9
  control = 6
```

LLM 队列：

```text
llm_review_queue.tsv rows = 124 + header
final_triage.llm.jsonl rows = 124
```

LLM response：

```text
final_triage.llm_responses.jsonl rows = 124

response final_label:
  suspicious = 75
  no_signal = 24
  correct = 18
  wrong = 7
```

终版表：

```text
final_adt_triage_table.llm_resolved.tsv rows = 175 + header

final_label:
  suspicious = 75
  correct = 54
  no_signal = 33
  wrong = 7
  control = 6

decision_source:
  llm = 124
  rule:target_enriched = 36
  rule:no_signal = 9
  rule:control = 6
```

## 12. 代表条目

```text
mAOC-CD16:
  final_label = correct
  decision_source = rule:target_enriched
  解释：Step 6 为 rare_target_enriched_low_count，方向与 Monocyte/NK/T 先验一致。

mAOC-CD16a:
  final_label = correct
  decision_source = llm
  解释：GMM threshold failed，但 raw ADT signal 在 expected target cells 中富集，且 observed ADT groups 与先验 overlap。

mAOC-CD34:
  final_label = correct
  decision_source = rule:target_enriched
  解释：CD34 在 HSPC 分选数据中多数细胞阳性符合实验设计。

mAOC-CD19:
  final_label = no_signal
  decision_source = llm
  解释：expected target groups 存在，但 ADT threshold/raw TNR 不支持 target enrichment。

mAOC-CD23:
  final_label = wrong
  decision_source = llm
  解释：RNA/validation 支持 B 细胞，但 ADT-positive groups 与 expected B 不 overlap。

mAOC-CD49f:
  final_label = wrong
  decision_source = llm
  解释：validation 支持 HSC_MPP，但 ADT-positive groups 不 overlap expected groups。

HTO / spike:
  final_label = control
  decision_source = rule:control
```

## 13. 下游使用建议

推荐默认使用：

```text
ADTriage_output/07_final_triage/final_adt_triage_table.llm_resolved.tsv
```

一般解释：

```text
correct:
  可作为当前 workflow 中较可信的 ADT feature。

suspicious:
  不建议作为强 target marker；可保留用于探索或后续单独复核。

wrong:
  不建议用于对应 target 的生物学解释，疑似 off-target、失败或错配。

no_signal:
  不建议用于生物学解释。

control:
  作为技术对照处理，不进入 biological target correctness 排名。
```

