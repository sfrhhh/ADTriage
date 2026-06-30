# Step 5: Marker Prior Review 可复现操作说明

## 1. 目标

Step 5 的目标是对每个非 control / spike / hashtag ADT feature 建立先验表达审查表，用于判断该 marker 在当前生物学分组中理论上应当支持哪些群体。

本步骤只对具有生物学含义的分组列运行。本次使用：

```text
Adult_cHSPCs_Ultima_majority_voting
```

`Standardpcaclusters` 等数字 cluster 不在本步骤中审查，因为数字 cluster 本身没有直接生物学语义。

## 2. 输入文件

```text
ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv
ADTriage_output/03_group_summary/adt_group_summary.tsv
ADTriage_output/03_group_summary/rna_group_summary.tsv
ADTriage_output/03_group_summary/cd45isoform_group_summary.tsv
priors/adult_cHSPC_marker_prior.strict.tsv
```

## 3. 使用脚本

```text
scripts/build_marker_prior_review.py
scripts/apply_marker_prior_llm_review.py
```

`build_marker_prior_review.py` 功能：

```text
1. 读取 Step 1 的 feature mapping。
2. 读取 Step 3 的 group-level ADT/RNA/CD45isoform summary。
3. 仅保留用户指定的 biological group_by。
4. 对非 control / spike / hashtag feature 计算 observed ADT-high groups。
5. 对 RNA/CD45isoform 验证特征计算 observed validation-high groups。
6. 用 strict local marker prior TSV 判断 expected_positive_groups。
7. 对无本地先验的 marker 生成 LLM review JSONL 队列。
```

`apply_marker_prior_llm_review.py` 功能：

```text
1. 读取 marker_prior_review.tsv。
2. 读取 LLM response JSONL。
3. 将 LLM 返回的 expected_positive_groups / confidence / reason 回填到 no-local-prior 条目。
4. 将回填后的 prior_method 改为 llm_marker_prior 或 llm_marker_prior:no_vote。
5. 重新计算 expected_adt_overlap_groups / expected_validation_overlap_groups。
6. 重新计算 needs_llm_review 与 needs_manual_review。
```

## 4. Strict Prior 原则

本步骤采用严格先验，不把“表达在某类细胞上”直接等价于“可作为该类群的先验正证据”。

规则：

```text
1. lineage-defining 或 lineage-informative marker 才能给 label 提供先验正证据。
2. broad marker 不提供先验正证据。
3. activation / checkpoint / adhesion marker 不提供先验正证据。
4. pDC / DC / basophil / mast / endothelial-context marker 不强行归入现有 label。
5. mature monocyte marker 不强行支持 GMP-L。
6. HSC_MPP 只保留相对明确的 primitive/progenitor marker。
```

先验表独立存放在：

```text
priors/adult_cHSPC_marker_prior.strict.tsv
```

字段含义：

```text
marker_symbol:
  先验 marker symbol，脚本会忽略大小写和非字母数字字符进行匹配。

expected_positive_groups:
  该 marker 在当前 label 集合中可支持的 expected positive groups，多个 group 用 ; 分隔。

prior_confidence:
  对本地先验条目的置信度。

prior_action:
  vote 或 no_vote。
  vote 表示 expected_positive_groups 可作为先验正证据。
  no_vote 表示该 marker 不适合作为当前 label 集合的先验正证据。

prior_evidence:
  简短依据说明。
```

后续如需更新先验，只需要修改该 TSV，不需要修改 Python 脚本。

## 5. no-vote 与 LLM review

本步骤区分三类情况：

```text
local_marker_prior
local_marker_prior:no_vote
llm_needed:no_local_prior
```

含义：

```text
local_marker_prior:
  有严格本地先验，且 expected_positive_groups 非空。

local_marker_prior:no_vote:
  有严格本地先验证据，但该 marker 不应在当前 label 集合中提供先验正证据。
  例如 broad marker、activation marker、pDC/DC-like marker。
  这类不进入 LLM 队列，默认也不自动进入人工审查。
  下游可继续使用 observed_adt_positive_groups / observed_validation_high_groups 作为结构化表达结果。

llm_needed:no_local_prior:
  没有本地先验命中，需要后续 LLM 或人工补充。
```

## 6. observed high group 定义

每个 feature 在每个 group_by 下独立计算。

默认使用：

```text
top_quantile = 0.75
```

逻辑：

```text
observed_adt_positive_groups:
  mean_CLR >= 当前 feature 在所有 group level 中的 75% 分位数，且 mean_CLR > 0。

observed_validation_high_groups:
  RNA 或 CD45isoform 验证信号同理使用 75% 分位数。
```

## 7. 复现命令

在项目根目录运行：

```bash
cd /data/project/yanhuixu/frsui/sc_project/ADTriage
```

### 7.1 Python 语法检查

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B -m py_compile \
  scripts/build_marker_prior_review.py \
  scripts/apply_marker_prior_llm_review.py
```

### 7.2 生成 prior review 表

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/build_marker_prior_review.py \
  --mapping ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv \
  --group-summary-dir ADTriage_output/03_group_summary \
  --outdir ADTriage_output/05_prior_review \
  --prior-file priors/adult_cHSPC_marker_prior.strict.tsv \
  --group-by Adult_cHSPCs_Ultima_majority_voting \
  --top-quantile 0.75
```

### 7.3 生成 LLM response JSONL

将 `marker_prior_review.llm.jsonl` 传递给 LLM。每条 response 必须至少包含：

```json
{
  "feature_id": "mAOC0000",
  "group_by": "Adult_cHSPCs_Ultima_majority_voting",
  "expected_positive_groups": ["Monocyte"],
  "confidence": 0.85,
  "reason": "short rationale"
}
```

本次输出文件为：

```text
ADTriage_output/05_prior_review/marker_prior_review.llm_responses.jsonl
```

约定：

```text
expected_positive_groups 非空:
  该 LLM response 合并后写作 llm_marker_prior。

expected_positive_groups 为空:
  该 LLM response 合并后写作 llm_marker_prior:no_vote。
  含义是该 marker 不适合作为当前 label 集合的先验正证据，但 observed groups 仍保留给下游。

confidence < 0.8:
  合并后仍标记 needs_manual_review = TRUE。
```

### 7.4 合并 LLM response

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/apply_marker_prior_llm_review.py \
  --review-table ADTriage_output/05_prior_review/marker_prior_review.tsv \
  --llm-responses ADTriage_output/05_prior_review/marker_prior_review.llm_responses.jsonl \
  --out ADTriage_output/05_prior_review/marker_prior_review.llm_resolved.tsv \
  --confidence-threshold 0.8
```

## 8. 输出文件

```text
ADTriage_output/05_prior_review/marker_prior_review.tsv
ADTriage_output/05_prior_review/marker_prior_review.llm.jsonl
ADTriage_output/05_prior_review/marker_prior_review.llm_responses.jsonl
ADTriage_output/05_prior_review/marker_prior_review.llm_resolved.tsv
```

### 8.1 marker_prior_review.tsv

关键字段：

```text
feature_id
feature_name
target_class
human_gene_symbol
validation_assay
validation_feature
group_by
group_levels
expected_positive_groups
prior_method
prior_confidence
prior_evidence
observed_adt_positive_groups
observed_validation_high_groups
adt_validation_overlap_groups
expected_adt_overlap_groups
expected_validation_overlap_groups
needs_llm_review
needs_manual_review
review_reason
```

### 8.2 marker_prior_review.llm.jsonl

只包含 `llm_needed:no_local_prior` 的条目。

每行是一个 JSON prompt record，包含：

```text
feature_id
feature_name
human_gene_symbol
validation_feature
group_by
group_levels
observed_adt_positive_groups
observed_validation_high_groups
required_output_schema
```

### 8.3 marker_prior_review.llm_responses.jsonl

包含 LLM 对 `marker_prior_review.llm.jsonl` 中每条记录的结构化 response。

关键字段：

```text
feature_id
group_by
expected_positive_groups
confidence
reason
reviewer
response_type
```

### 8.4 marker_prior_review.llm_resolved.tsv

这是 Step 5 的 LLM 合并后主输出表。

相对于 `marker_prior_review.tsv`：

```text
llm_needed:no_local_prior -> llm_marker_prior 或 llm_marker_prior:no_vote
needs_llm_review -> FALSE
expected_positive_groups / prior_confidence / prior_evidence 被 LLM response 回填
needs_manual_review 根据 LLM confidence 和 observed overlap 重新计算
```

## 9. 本次运行统计

```text
marker_prior_review.tsv rows = 169 + header
marker_prior_review.llm.jsonl rows = 53
marker_prior_review.llm_responses.jsonl rows = 53
marker_prior_review.llm_resolved.tsv rows = 169 + header
Standardpcaclusters residual rows = 0
```

原始 prior_method 统计：

```text
local_marker_prior = 89
local_marker_prior:no_vote = 27
llm_needed:no_local_prior = 53
```

合并后 prior_method 统计：

```text
local_marker_prior = 89
local_marker_prior:no_vote = 27
llm_marker_prior = 21
llm_marker_prior:no_vote = 32
```

合并后 needs_llm_review：

```text
FALSE = 169
TRUE = 0
```

合并后 needs_manual_review：

```text
FALSE = 133
TRUE = 36
```

LLM 回填后仍需人工检查的条目：

```text
mAOC-CD96 / CD96:
  LLM expected = NK/T
  reason = ADT-high groups did not overlap LLM expected groups.

mAOC-CD218a / IL18R1:
  LLM expected = Monocyte;NK/T
  reason = validation-high groups did not overlap LLM expected groups.

mAOC-IL1RAP / IL1RAP:
  LLM expected = GMP-L;HSC_MPP;MEBEMP-L
  reason = LLM confidence below 0.8.
```

## 10. 关键 strict prior 变更

以下 marker 在 strict prior 中不提供当前 label 的先验正证据，并标记为 `local_marker_prior:no_vote`：

```text
CD80
CD1D
CD49d / ITGA4
CD49b / ITGA2
CD155 / PVR
HLA-ABC / HLA-A
CD58
CD303 / CLEC4C
CD31 / PECAM1
GPC2
CD107a / LAMP1
TIM3 / HAVCR2
CD202b / TEK
CD46
FCER1A
CD50 / ICAM3
CD316 / IGSF8
CD130 / IL6ST
CD98 / SLC3A2
CD84
CD47
CD29 / ITGB1
CD326 / EPCAM
CD146 / MCAM
CD304 / NRP1
CD18 / ITGB2
CD39 / ENTPD1
```

部分 marker 被收窄或重分配：

```text
ALCAM -> HSC_MPP
ENG -> MEBEMP-L;ERYP
ICAM2 -> HSC_MPP
ITGA6 -> HSC_MPP
KIT -> HSC_MPP;MEBEMP-L;GMP-L
CD38 -> CLP;B;MEBEMP-L;GMP-L
TFRC -> ERYP
CD93 -> HSC_MPP;Monocyte
BST1 -> Monocyte
FCAR -> Monocyte
FCGR1A -> Monocyte
FCGR2A -> Monocyte
ITGAX -> Monocyte
LILRA5 -> Monocyte
LILRB3 -> Monocyte
SIRPA -> Monocyte
CCR7 -> NK/T
CXCR5 -> B;NK/T
SELL -> CLP;B;NK/T
```

## 11. 需要手动检查的部分

每次重新运行 Step 5 后，建议检查：

```text
1. local_marker_prior:no_vote 条目是否确实不应提供当前 label 的先验正证据。
2. needs_manual_review = TRUE 的本地先验冲突项是否需要调整 expected_positive_groups。
3. marker_prior_review.llm_resolved.tsv 中 llm_marker_prior 条目的 expected_positive_groups 是否符合生物学判断。
4. marker_prior_review.llm_resolved.tsv 中 llm_marker_prior:no_vote 条目是否确实只应保留 observed groups。
5. 如果后续加入 pDC/DC/basophil/mast/endothelial-like label，应重新评估 CLEC4C、IL3RA、NRP1、FCER1A 等 marker。
6. 如果修改 strict prior TSV，需要保证本说明中的 no-vote 语义被保留。
```
