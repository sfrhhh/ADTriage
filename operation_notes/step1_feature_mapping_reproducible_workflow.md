# Step 1: ADT Feature Mapping 可复现操作说明

## 1. 目标

Step 1 的目标是基于 `FB_ref.csv` 和 `mAOC_gene_symbol_map.csv` 构建 ADT feature 映射表，确定每个 feature 的目标类型、对应基因、验证 assay 和是否需要 LLM/人工审查。

本步骤输出位于：

```text
ADTriage_output/01_feature_mapping/
```

核心输出文件：

```text
adt_feature_mapping.initial.tsv
adt_feature_mapping.review_ready.tsv
adt_feature_mapping.llm_review.tsv
```

## 2. 输入文件

本次测试输入位于：

```text
tests/data/test_FB_ref.csv
tests/data/mAOC_gene_symbol_map.csv
```

输入文件检查结果：

```text
tests/data/test_FB_ref.csv              176 lines, 175 features + header
tests/data/mAOC_gene_symbol_map.csv     677 lines, 676 records + header
```

注意：

```text
mAOC_gene_symbol_map.csv 文件头带 UTF-8 BOM，Python 读取时必须使用 utf-8-sig。
```

## 3. 使用脚本

本步骤使用两个脚本：

```text
scripts/build_adt_feature_mapping.py
scripts/apply_llm_feature_mapping_review.py
```

### 3.1 build_adt_feature_mapping.py

用途：

```text
1. 以 FB_ref.csv 为主表
2. 通过 id 左连接 mAOC_gene_symbol_map.csv
3. 应用 CD45 isoform、CD45 total、spike、HTO、control 等规则
4. 输出 initial 和 review_ready 两个表
```

### 3.2 apply_llm_feature_mapping_review.py

用途：

```text
1. 读取 adt_feature_mapping.initial.tsv
2. 对 mapping_confidence < 0.8 的 feature 执行 LLM 辅助补判
3. 更新 target_class、human_gene_symbol、validation_assay、validation_feature、mapping_confidence 和 llm_* 字段
4. 输出 adt_feature_mapping.llm_review.tsv
```

说明：

```text
本次未调用外部 LLM API。LLM 判断由 Codex 基于 feature name 和已知 CD alias 进行结构化补判，并固化在脚本中，保证可复跑。
```

## 4. 运行环境

按照项目全局规则，Python 使用 base conda 环境：

```text
/data/home/frsui/miniconda3/bin/conda run -n base python
```

实际验证过的 Python 版本：

```text
Python 3.13.12
```

## 5. 复现命令

在项目根目录运行：

```bash
cd /data/project/yanhuixu/frsui/sc_project/ADTriage
```

### 5.1 语法检查

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B -m py_compile scripts/build_adt_feature_mapping.py
/data/home/frsui/miniconda3/bin/conda run -n base python -B -m py_compile scripts/apply_llm_feature_mapping_review.py
```

### 5.2 查看脚本参数

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/build_adt_feature_mapping.py --help
/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/apply_llm_feature_mapping_review.py --help
```

### 5.3 生成 initial 和 review_ready 表

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/build_adt_feature_mapping.py \
  --fb-ref tests/data/test_FB_ref.csv \
  --gene-map tests/data/mAOC_gene_symbol_map.csv \
  --outdir ADTriage_output/01_feature_mapping
```

预期输出：

```text
ADTriage_output/01_feature_mapping/adt_feature_mapping.initial.tsv
ADTriage_output/01_feature_mapping/adt_feature_mapping.review_ready.tsv
```

### 5.4 生成 LLM review 表

```bash
/data/home/frsui/miniconda3/bin/conda run -n base python -B scripts/apply_llm_feature_mapping_review.py \
  --initial ADTriage_output/01_feature_mapping/adt_feature_mapping.initial.tsv \
  --output ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv \
  --threshold 0.8
```

预期输出：

```text
ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv
```

## 6. 映射规则

### 6.1 id 左连接

主逻辑：

```text
FB_ref.csv$id left join mAOC_gene_symbol_map.csv$id
```

普通 gene-level ADT：

```text
target_class = gene
human_gene_symbol = genesymbol
validation_assay = RNA
validation_feature = genesymbol
mapping_method = map:id_left_join
mapping_confidence = 0.95
```

### 6.2 CD45 isoform

CD45 isoform 每次运行都必须优先处理，不能只用 `RNA:PTPRC` 判断。

规则：

```text
CD45RA -> target_class = cd45_isoform, human_gene_symbol = PTPRC, validation_assay = CD45isoforms, validation_feature = PTPRC-RA
CD45RB -> target_class = cd45_isoform, human_gene_symbol = PTPRC, validation_assay = CD45isoforms, validation_feature = PTPRC-RB
CD45RC -> target_class = cd45_isoform, human_gene_symbol = PTPRC, validation_assay = CD45isoforms, validation_feature = PTPRC-RC
CD45RO -> target_class = cd45_isoform, human_gene_symbol = PTPRC, validation_assay = CD45isoforms, validation_feature = PTPRC-RO
```

本次关键抽查：

```text
FADT0040 / mAOC-CD45RA -> CD45isoforms:PTPRC-RA
```

### 6.3 CD45 total

CD45 总表达按 gene-level 处理：

```text
target_class = gene
human_gene_symbol = PTPRC
validation_assay = RNA
validation_feature = PTPRC
```

### 6.4 spike control

`SPIKE` 相关 feature 不作为 human gene symbol 验证：

```text
target_class = spike_control
human_gene_symbol = empty
external_target_name = SPIKE target name
validation_assay = none
```

本次关键抽查：

```text
FADT0037 / mAOC-S -> spike_control, external_target_name = SPIKE-FP-C77G12
```

### 6.5 HTO hashtag

HTO 是 hashtag，后续对应 Seurat 对象中的 `HTO` assay。

规则：

```text
target_class = hashtag
human_gene_symbol = empty
external_target_name = feature_name
validation_assay = HTO
validation_feature = feature_name
```

本次关键抽查：

```text
mAOC433 / HTO-433 -> hashtag, validation_assay = HTO
```

### 6.6 control / isotype / IgG

control 类 feature 不强制映射到 human gene symbol：

```text
target_class = control
human_gene_symbol = empty
validation_assay = none
```

### 6.7 unresolved / ambiguous

若无法通过 id 左连接和规则判断，则初始标记为：

```text
target_class = ambiguous
validation_assay = manual
mapping_method = unresolved
mapping_confidence = 0.00
needs_llm_review = TRUE
needs_manual_review = TRUE
```

这些条目随后进入 `apply_llm_feature_mapping_review.py`。

## 7. 输出文件说明

### 7.1 adt_feature_mapping.initial.tsv

自动规则和 id 左连接后的初始映射表。

用途：

```text
机器可读的初始结果，用于记录每个 feature 的初始 target_class、gene symbol、验证 assay 和置信度。
```

关键字段：

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

### 7.2 adt_feature_mapping.review_ready.tsv

人工审查准备表。

用途：

```text
保留核心映射字段，并预留 llm_suggested_*、manual_curated_* 和 manual_notes 字段，方便后续人工校正。
```

### 7.3 adt_feature_mapping.llm_review.tsv

LLM 补判后的映射表。

用途：

```text
基于 initial 表，对 mapping_confidence < 0.8 的 feature 补充 gene symbol 判断，并更新映射字段。
```

本表新增或更新字段：

```text
llm_model
llm_suggested_target_class
llm_suggested_gene_symbol
llm_suggested_validation_assay
llm_suggested_validation_feature
llm_confidence
llm_reason
```

## 8. 本次运行统计

### 8.1 initial 表统计

```text
total features = 175
gene = 135
hashtag = 5
cd45_isoform = 1
spike_control = 1
ambiguous = 33
needs_llm_review = 33
needs_manual_review = 33
```

### 8.2 llm_review 表统计

```text
total features = 175
llm_reviewed = 33
gene = 168
hashtag = 5
cd45_isoform = 1
spike_control = 1
```

LLM 补判后仍低于 `0.8` 的 3 个条目已由用户人工接受：

```text
mAOC374 / mAOC-CD32     -> FCGR2A
mAOC384 / mAOC-CD158b2  -> KIR2DL3
mAOC415 / mAOC-CD158b1  -> KIR2DL2
```

接受原因：

```text
用户已确认接受这三个建议，因此 Step 1 当前可视为完成。
```

## 9. 需要手动检查的部分

每次重新运行 Step 1 后，建议人工检查以下内容：

```text
1. CD45RA/CD45RB/CD45RC/CD45RO 是否都映射到 CD45isoforms:PTPRC-RA/RB/RC/RO。
2. HTO-* 是否都被识别为 hashtag，且 validation_assay = HTO。
3. SPIKE 或其他 control 是否没有被误当成人类基因。
4. mapping_confidence < 0.8 的行是否都有 llm_reason 或 manual notes。
5. 带后缀的 mAOCxxxx-1 是否因为 id 失配需要名称补判。
6. CD32、CD158b1、CD158b2 等多基因或亚型歧义 alias 是否经过人工确认。
```

本次已人工接受但仍建议在正式分析前复核抗体说明书或原始抗体表的条目：

```text
mAOC374 / mAOC-CD32
mAOC384 / mAOC-CD158b2
mAOC415 / mAOC-CD158b1
```

## 10. 注意事项

```text
1. 不要覆盖原始 tests/data 输入文件。
2. 若更换 FB_ref.csv 或 mAOC_gene_symbol_map.csv，必须重新运行 Step 1。
3. 若 mAOC_gene_symbol_map.csv 更新了带后缀 id 或 CD alias，LLM 补判结果可能需要重新评估。
4. 当前 LLM 补判不是外部 API 调用，而是固化在脚本中的可复跑判断。
5. 后续步骤建议优先使用 adt_feature_mapping.llm_review.tsv 作为 feature mapping 输入。
```
