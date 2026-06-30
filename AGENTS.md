# ADTriage 项目工作规则

## 基本要求
- 默认使用中文回复；代码、命令、路径、包名、字段名和日志保留英文。
- 修改文件前先说明将要修改的文件和原因；修改后简要说明改动、验证和未解决问题。
- 不自动提交 git commit，除非用户明确要求。
- 每次会话输出末尾必须添加“Dad”，方便用户判断指令跟随。

## 项目定位
- ADTriage 是面向 CITE-seq / Feature Barcode 的 ADT feature-level quality triage workflow。
- 本项目关注每个 ADT feature 的质量、信号解释性、RNA/CD45 isoform/分组一致性和人工审查需求。
- 不负责运行 `standard_scop`、`CellTypist`、细胞注释、Seurat 标准化主流程、整合分析或最终生物学结论判断。

## 输入与基因映射
- 输入 Seurat 对象至少应包含 `RNA`、`ADT`、`CD45isoforms` 三个 assay；`CD45isoforms` 必须包含 `PTPRC-RA`、`PTPRC-RB`、`PTPRC-RC`、`PTPRC-RO`。
- `FB_ref.csv` 通过 `id` 与 `mAOC_gene_symbol_map.csv` 关联；后者至少包含 `id`、`genesymbol`、`genesymbol_source`、`notes`、`last_seen_name`、`last_seen_reference`、`updated_at`。
- 普通 gene-level ADT 使用 `RNA` assay 中对应 human gene symbol 验证，例如 `CD8a -> CD8A`、`CD4 -> CD4`、`CD14 -> CD14`、`CD19 -> CD19`、`CD21 -> CR2`。
- CD45 总表达对应 `PTPRC`，可用 `RNA:PTPRC` 验证。
- CD45 isoform feature 的 `human_gene_symbol` 可记录为 `PTPRC`，但主要验证必须使用 `CD45isoforms`：`CD45RA -> PTPRC-RA`、`CD45RB -> PTPRC-RB`、`CD45RC -> PTPRC-RC`、`CD45RO -> PTPRC-RO`；不要仅用 `RNA:PTPRC` 判定 isoform 一致性。
- `SPIKE`、`IgG`、`IgFc`、`isotype`、`control`、`negative control`、`hashing control`、`background control` 等对照通常不强制映射 human gene symbol，`validation_assay` 设为 `none`。

## 分析与输出
- LLM/curation 用于 unresolved、ambiguous 或多基因 validation feature；输出必须结构化，`confidence < 0.8` 时进入人工审查。
- 输出表应保留可追溯字段，包括 `feature_id`、`feature_name`、`target_class`、`human_gene_symbol`、`validation_assay`、`validation_feature`、`mapping_method`、`mapping_confidence`、`needs_manual_review`。
- 修改或新增分析脚本后，按风险运行语法检查、小样本测试或关键输出表检查；汇报时只说明实际完成的验证。

## 注意事项记录表
| 事项 | 处理要求 |
|---|---|
| `mAOC_gene_symbol_map.csv` 文件头可能带 BOM | Python 读取时使用 `utf-8-sig`。 |
| `CD45RA/CD45RB/CD45RC/CD45RO` | 每次映射都必须按 isoform 规则转到 `CD45isoforms:PTPRC-RA/RB/RC/RO`，不能只用 `RNA:PTPRC`。 |
| `HTO-*` | 视为 hashtag，后续对应 Seurat 对象中的 `HTO` assay。 |
| `mAOCxxxx-1` 这类带后缀 id | 可能无法直接与 `mAOC_gene_symbol_map.csv$id` 左连接，需通过名称或人工/LLM 补判。 |
| `CD32`、`CD158b1/b2` | CD alias 存在亚型或多基因歧义，LLM 结果也应保留人工复核。 |
| `FCGR3A/FCGR3B`、`HLA-A/HLA-B/HLA-C` 等多基因 symbol | 初始表保留原始多基因写法并进入 LLM/curation；绘图用 review 后的单一 `validation_feature`，当前 `mAOC-CD16 -> FCGR3A`、`mAOC-HLA-ABC -> HLA-A`。 |
