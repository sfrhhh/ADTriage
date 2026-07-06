---
name: adtriage
description: Run or adapt the ADTriage workflow for CITE-seq and Feature Barcode Seurat objects, including ADT antibody quality, ADT feature triage, ADT/RNA consistency, CD45 isoform validation, marker prior review, Seurat object checks, new ADT panel evaluation, and final ADT feature-level quality labels. Use when the user asks to run ADTriage, use $adtriage, evaluate CITE-seq antibody quality, triage Feature Barcode ADT features, review ADT/RNA consistency, process a new ADT panel, or repeat this repo's ADT quality workflow on an existing Seurat object. This skill does not run upstream Cell Ranger, does not run the main Seurat preprocessing/integration workflow, and does not perform cell annotation; it operates on an already prepared Seurat object.
---

# ADTriage Skill

Use this skill to run the repo-scoped ADTriage workflow on a new CITE-seq / Feature Barcode Seurat object. The goal is to reuse the existing project scripts and outputs, not to rewrite the ADTriage algorithms.

## Core Purpose

ADTriage produces feature-level triage labels for ADT antibodies by combining feature-reference mapping, RNA/CD45isoform validation, UMAP feature plots, group summaries, ADT-vs-validation consistency, strict marker prior review, ADT distribution metrics, TNR, target enrichment, and final rule-gated/LLM-resolved labels.

The final result must be:

```text
<outdir>/07_final_triage/final_adt_triage_table.llm_resolved.tsv
```

## Required Inputs

Require these files from the user or from a project config:

```yaml
seurat_rds: path/to/prepared_seurat_object.rds
fb_ref: path/to/FB_ref.csv
gene_symbol_map: path/to/mAOC_gene_symbol_map.csv
outdir: path/to/ADTriage_output
```

The Seurat object must already contain assays, metadata, and reductions. Do not run upstream Cell Ranger, standard Seurat preprocessing, integration, or cell annotation.

Required assays:

```text
RNA
ADT
CD45isoforms
```

`CD45isoforms` must include:

```text
PTPRC-RA
PTPRC-RB
PTPRC-RC
PTPRC-RO
```

Required mapping columns:

```text
FB_ref.csv: id, name, read, pattern, sequence, feature_type
mAOC_gene_symbol_map.csv: id, genesymbol, genesymbol_source, notes, last_seen_name, last_seen_reference, updated_at
```

Read `mAOC_gene_symbol_map.csv` with `utf-8-sig` behavior because the header may contain a BOM.

## Optional Inputs And Parameters

Use these when provided by the user:

```yaml
prior_file: priors/adult_cHSPC_marker_prior.strict.tsv
step5_llm_responses: <outdir>/05_prior_review/marker_prior_review.llm_responses.jsonl
step7_llm_responses: <outdir>/07_final_triage/final_triage.llm_responses.jsonl
n_features: all mapping rows, or a smaller test number
panel_px: 800
panel_dpi: 200
pt_size: 0.5
min_groups: 3
top_quantile: 0.75
```

If optional LLM responses are not yet available, generate the queue JSONL first and stop for the user to provide or approve LLM response generation. Do not invent final labels in chat.

## Metadata And Reduction Configuration

Use these exact terms. Do not use `plot_group_by` as a generic umbrella term.

```yaml
metadata:
  biological_group_by: Adult_cHSPCs_Ultima_majority_voting
  summary_group_by:
    - Standardpcacluster
    - Adult_cHSPCs_Ultima_majority_voting
  feature_plot_group_by:
    - Standardpcacluster
    - Adult_cHSPCs_Ultima_majority_voting

reduction: StandardpcaUMAP2D
```

Definitions:

- `biological_group_by`: exactly one metadata column used for biological target/non-target interpretation in Step 5, Step 6, and Step 7.
- `summary_group_by`: one or more metadata columns used for Step 3 group summaries and Step 4 ADT/RNA or ADT/CD45isoform correlation plots.
- `feature_plot_group_by`: optional future setting for Step 2 metadata grouping. The current `scripts/plot_adt_validation_features.R` does not support a metadata group-by argument. Do not pass unsupported parameters.
- `reduction`: a Seurat reduction from `Reductions(obj)`, used by Step 2 feature plots. This is not a metadata column.

If a requested metadata column is missing, stop and report available metadata columns. Do not guess a similar column name.

If a requested reduction is missing, stop and report available reductions. Do not guess a replacement. In particular, if the user requests `StandardpcaUMAP2D` and it is not present in `Reductions(obj)`, stop.

## Environment

Use the repo's existing environments and tools:

```bash
PY="/data/home/frsui/miniconda3/bin/conda run -n base python -B"
RS="/data/home/frsui/miniconda3/bin/conda run -n scop Rscript"
MAGICK="/data/project/yanhuixu/frsui/software/bin/magick"
```

Do not install packages and do not modify conda environments unless the user explicitly asks.

## Mandatory Input Inspection

Before running any full workflow on a new RDS, inspect object structure:

```bash
/data/home/frsui/miniconda3/bin/conda run -n scop Rscript -e '
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("path/to/object.rds")
cat("Assays:\n"); print(Assays(obj))
cat("Reductions:\n"); print(Reductions(obj))
cat("Metadata columns:\n"); print(colnames(obj@meta.data))
'
```

Then verify:

- `RNA`, `ADT`, and `CD45isoforms` are in `Assays(obj)`.
- `PTPRC-RA`, `PTPRC-RB`, `PTPRC-RC`, `PTPRC-RO` are in `CD45isoforms`.
- configured `reduction` is in `Reductions(obj)`.
- configured `biological_group_by` is in `colnames(obj@meta.data)`.
- all configured `summary_group_by` columns are in `colnames(obj@meta.data)`.
- `feature_plot_group_by` columns are only checked if a future plotting script supports them.

If any check fails, stop and report the actual available assays/reductions/metadata. Do not continue with guessed names.

## Workflow Order

Run steps in order. Prefer a clean, user-selected output directory for each new dataset, for example `ADTriage_output_<dataset_id>`.

### Step 1: ADT Feature Mapping

Scripts:

```text
scripts/build_adt_feature_mapping.py
scripts/apply_llm_feature_mapping_review.py
```

Commands:

```bash
$PY scripts/build_adt_feature_mapping.py \
  --fb-ref "$FB_REF" \
  --gene-map "$GENE_SYMBOL_MAP" \
  --outdir "$OUTDIR/01_feature_mapping"

$PY scripts/apply_llm_feature_mapping_review.py \
  --initial "$OUTDIR/01_feature_mapping/adt_feature_mapping.initial.tsv" \
  --output "$OUTDIR/01_feature_mapping/adt_feature_mapping.llm_review.tsv" \
  --threshold 0.8
```

Check:

```text
01_feature_mapping/adt_feature_mapping.initial.tsv
01_feature_mapping/adt_feature_mapping.review_ready.tsv
01_feature_mapping/adt_feature_mapping.llm_review.tsv
```

Confirm CD45 isoforms map to `CD45isoforms:PTPRC-RA/RB/RC/RO`, not only `RNA:PTPRC`. Confirm multi-gene features such as `FCGR3A/FCGR3B` or `HLA-A/HLA-B/HLA-C` are resolved to a single validation feature before plotting.

### Step 2: ADT Validation Feature Plots

Script:

```text
scripts/plot_adt_validation_features.R
```

Command:

```bash
$RS scripts/plot_adt_validation_features.R \
  --seurat-rds "$SEURAT_RDS" \
  --mapping "$OUTDIR/01_feature_mapping/adt_feature_mapping.llm_review.tsv" \
  --outdir "$OUTDIR/02_feature_plots" \
  --n-features <number_of_mapping_rows> \
  --reduction "$REDUCTION" \
  --panel-px 800 \
  --panel-dpi 200 \
  --pt-size 0.5 \
  --max-page-height-px 9600 \
  --magick "$MAGICK"
```

Do not pass metadata group-by parameters to this script; it currently does not support `feature_plot_group_by`.

Check:

```text
02_feature_plots/adt_validation_plot_feature_*_page*.png
02_feature_plots/adt_validation_plot_layout_index_*.tsv
02_feature_plots/adt_validation_plot_status_*.tsv
```

Use the layout index TSV to locate features in the large PNG pages.

### Step 3: Group-Level Summary

Script:

```text
scripts/build_adt_group_summary.R
```

Command:

```bash
SUMMARY_GROUP_BY_CSV="$(IFS=,; echo "${SUMMARY_GROUP_BY[*]}")"

$RS scripts/build_adt_group_summary.R \
  --seurat-rds "$SEURAT_RDS" \
  --mapping "$OUTDIR/01_feature_mapping/adt_feature_mapping.llm_review.tsv" \
  --outdir "$OUTDIR/03_group_summary" \
  --group-by "$SUMMARY_GROUP_BY_CSV"
```

Check:

```text
03_group_summary/adt_group_summary.tsv
03_group_summary/rna_group_summary.tsv
03_group_summary/cd45isoform_group_summary.tsv
```

Verify rows exist for every requested `summary_group_by` column and that group cell counts match the Seurat metadata.

### Step 4: Signal Consistency

Script:

```text
scripts/calculate_signal_consistency.R
```

Command:

```bash
$RS scripts/calculate_signal_consistency.R \
  --group-summary-dir "$OUTDIR/03_group_summary" \
  --outdir "$OUTDIR/04_signal_consistency" \
  --min-groups 3
```

Check:

```text
04_signal_consistency/rna_adt_correlation_summary.tsv
04_signal_consistency/cd45isoform_adt_correlation_summary.tsv
04_signal_consistency/rna_adt_<group_by>_consistency_scatter_plots.pdf
04_signal_consistency/cd45isoform_adt_<group_by>_consistency_scatter_plots.pdf
```

Step 4 uses all `summary_group_by` values created by Step 3.

### Step 5: Marker Prior Review

Scripts:

```text
scripts/build_marker_prior_review.py
scripts/apply_marker_prior_llm_review.py
```

Command:

```bash
$PY scripts/build_marker_prior_review.py \
  --mapping "$OUTDIR/01_feature_mapping/adt_feature_mapping.llm_review.tsv" \
  --group-summary-dir "$OUTDIR/03_group_summary" \
  --outdir "$OUTDIR/05_prior_review" \
  --prior-file priors/adult_cHSPC_marker_prior.strict.tsv \
  --group-by "$BIOLOGICAL_GROUP_BY" \
  --top-quantile 0.75
```

If `marker_prior_review.llm.jsonl` contains rows, send them to LLM and write responses to:

```text
05_prior_review/marker_prior_review.llm_responses.jsonl
```

Then merge:

```bash
$PY scripts/apply_marker_prior_llm_review.py \
  --review-table "$OUTDIR/05_prior_review/marker_prior_review.tsv" \
  --llm-responses "$OUTDIR/05_prior_review/marker_prior_review.llm_responses.jsonl" \
  --out "$OUTDIR/05_prior_review/marker_prior_review.llm_resolved.tsv" \
  --confidence-threshold 0.8
```

Check:

```text
05_prior_review/marker_prior_review.tsv
05_prior_review/marker_prior_review.llm.jsonl
05_prior_review/marker_prior_review.llm_responses.jsonl
05_prior_review/marker_prior_review.llm_resolved.tsv
```

Use only `biological_group_by` here. Do not run Step 5 on numeric clusters unless the user explicitly defines them as biological target labels.

Step 5 LLM response schema:

```json
{
  "feature_id": "mAOC0000",
  "group_by": "Adult_cHSPCs_Ultima_majority_voting",
  "expected_positive_groups": ["Monocyte"],
  "confidence": 0.85,
  "reason": "short rationale"
}
```

### Step 6: ADT Distribution Metrics

Script:

```text
scripts/build_adt_distribution_metrics.R
```

Command:

```bash
$RS scripts/build_adt_distribution_metrics.R \
  --seurat-rds "$SEURAT_RDS" \
  --mapping "$OUTDIR/01_feature_mapping/adt_feature_mapping.llm_review.tsv" \
  --prior-review "$OUTDIR/05_prior_review/marker_prior_review.llm_resolved.tsv" \
  --outdir "$OUTDIR/06_distribution_metrics" \
  --group-by "$BIOLOGICAL_GROUP_BY"
```

Check:

```text
06_distribution_metrics/adt_distribution_metrics.tsv
06_distribution_metrics/adt_gmm_metrics.tsv
06_distribution_metrics/adt_density_plots.pdf
```

Step 6 fits and plots CLR > 0 values by default, but positive rate and TNR are computed over all cells. TNR target cells are defined by Step 5 `expected_positive_groups` under `biological_group_by`.

### Step 7: Final ADT Triage

Script:

```text
scripts/build_final_adt_triage.py
```

Generate rule table and LLM queue:

```bash
$PY scripts/build_final_adt_triage.py \
  --mapping "$OUTDIR/01_feature_mapping/adt_feature_mapping.llm_review.tsv" \
  --prior-review "$OUTDIR/05_prior_review/marker_prior_review.llm_resolved.tsv" \
  --distribution "$OUTDIR/06_distribution_metrics/adt_distribution_metrics.tsv" \
  --rna-correlation "$OUTDIR/04_signal_consistency/rna_adt_correlation_summary.tsv" \
  --cd45-correlation "$OUTDIR/04_signal_consistency/cd45isoform_adt_correlation_summary.tsv" \
  --outdir "$OUTDIR/07_final_triage" \
  --group-by "$BIOLOGICAL_GROUP_BY"
```

If `07_final_triage/final_triage.llm.jsonl` contains rows, send them to LLM and write responses to:

```text
07_final_triage/final_triage.llm_responses.jsonl
```

Merge final LLM responses:

```bash
$PY scripts/build_final_adt_triage.py \
  --mapping "$OUTDIR/01_feature_mapping/adt_feature_mapping.llm_review.tsv" \
  --prior-review "$OUTDIR/05_prior_review/marker_prior_review.llm_resolved.tsv" \
  --distribution "$OUTDIR/06_distribution_metrics/adt_distribution_metrics.tsv" \
  --rna-correlation "$OUTDIR/04_signal_consistency/rna_adt_correlation_summary.tsv" \
  --cd45-correlation "$OUTDIR/04_signal_consistency/cd45isoform_adt_correlation_summary.tsv" \
  --outdir "$OUTDIR/07_final_triage" \
  --group-by "$BIOLOGICAL_GROUP_BY" \
  --llm-responses "$OUTDIR/07_final_triage/final_triage.llm_responses.jsonl"
```

Check:

```text
07_final_triage/final_adt_triage_table.tsv
07_final_triage/llm_review_queue.tsv
07_final_triage/final_triage.llm.jsonl
07_final_triage/final_triage.llm_responses.jsonl
07_final_triage/final_adt_triage_table.llm_resolved.tsv
```

The resolved table is the final result.

Step 7 LLM response schema:

```json
{
  "feature_id": "mAOC0000",
  "final_label": "correct",
  "final_confidence": 0.9,
  "final_reason": "short rationale",
  "recommended_action": "short downstream action"
}
```

Allowed `final_label` values: `correct`, `suspicious`, `wrong`, `no_signal`, `control`.

## LLM Review Rules

Do not let Codex decide final labels only in conversation.

All LLM decisions must be written to `*.llm_responses.jsonl`, then merged by the appropriate script.

For Step 5:

```text
marker_prior_review.llm.jsonl -> marker_prior_review.llm_responses.jsonl -> marker_prior_review.llm_resolved.tsv
```

For Step 7:

```text
final_triage.llm.jsonl -> final_triage.llm_responses.jsonl -> final_adt_triage_table.llm_resolved.tsv
```

Never edit `final_adt_triage_table.tsv` or `final_adt_triage_table.llm_resolved.tsv` by hand. Regenerate them by rerunning the merge step.

Only `final_adt_triage_table.llm_resolved.tsv` is the final triage output after Step 7 LLM review.

## Final Result Checklist

Before reporting completion, verify:

- `final_adt_triage_table.llm_resolved.tsv` exists.
- It has one row per mapped feature.
- `final_label` contains only `correct`, `suspicious`, `wrong`, `no_signal`, or `control`.
- no unresolved `llm_review_needed` remains unless the user explicitly accepts unresolved rows.
- control/HTO/spike features are labeled `control`.
- CD45 isoform features used `CD45isoforms` validation.
- representative rows such as CD16, CD34, CD23, controls, and any user-specified markers have sensible labels.
- all LLM response files are preserved in the output directory.

## Failure Troubleshooting Order

1. Check input paths and file existence.
2. Inspect `Assays(obj)`, `Reductions(obj)`, and `colnames(obj@meta.data)`.
3. Confirm `biological_group_by`, `summary_group_by`, `reduction`, assays, and CD45isoform features exactly exist.
4. Check Step 1 unresolved/multi-gene mappings and Step 2 plot status TSV.
5. Check Step 3 status/message columns and Step 4 `status = insufficient_groups`.
6. Check Step 5 prior rows/LLM response schema, Step 6 fit status/target counts, and Step 7 LLM response schema/confidence.

## Forbidden Actions

- Do not run Cell Ranger.
- Do not run standard Seurat preprocessing, integration, or cell annotation.
- Do not install packages or modify conda environments unless the user explicitly asks.
- Do not invent metadata column names or substitute similar-looking columns.
- Do not substitute a different reduction when the configured reduction is missing.
- Do not pass unsupported `feature_plot_group_by` arguments to `scripts/plot_adt_validation_features.R`.
- Do not manually edit final TSV result tables.
- Do not produce LLM triage decisions only in chat.
- Do not discard or overwrite LLM JSONL queues/responses.
- Do not automatically git commit.
- Do not treat `adt_signal_class` from Step 6 as the final label; Step 7 `final_label` is authoritative.

## Current Customization Gaps

These are known script-level limits to mention when relevant:

- `feature_plot_group_by` is not currently supported by `scripts/plot_adt_validation_features.R`; it is optional/future enhancement only.
- Step 5, Step 6, and Step 7 accept only one biological group column through `--group-by`.
- Step 4 uses all group summaries generated by Step 3; it does not take a separate metadata selector.
- Step 1 mapping review uses the current script's explicit internal LLM-assisted suggestions for low-confidence mappings, not a general external API call.
