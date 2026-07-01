#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path


VALID_LABELS = {
    "correct",
    "suspicious",
    "wrong",
    "no_signal",
    "control",
    "llm_review_needed",
}

CONTROL_TARGET_CLASSES = {
    "control",
    "spike_control",
    "hashtag",
}

TARGET_ENRICHED_CLASSES = {
    "target_enriched",
    "rare_target_enriched",
    "rare_target_enriched_low_count",
    "target_enriched_broad",
}

NO_SIGNAL_CLASSES = {
    "no_signal",
    "no_target_no_signal",
}

LLM_SIGNAL_CLASSES = {
    "distribution_detected_but_off_target",
    "ambiguous",
    "ambiguous_trimodal",
    "rare_or_ambiguous",
    "no_target_detected_signal",
    "no_target_broad_signal",
    "no_target_enrichment",
    "no_target_threshold_unavailable",
    "target_raw_enriched_no_threshold",
    "target_signal_unresolved",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build final ADT triage table and LLM review queue."
    )
    parser.add_argument(
        "--mapping",
        default="ADTriage_output/01_feature_mapping/adt_feature_mapping.llm_review.tsv",
        help="Step 1 resolved feature mapping TSV.",
    )
    parser.add_argument(
        "--prior-review",
        default="ADTriage_output/05_prior_review/marker_prior_review.llm_resolved.tsv",
        help="Step 5 resolved marker prior review TSV.",
    )
    parser.add_argument(
        "--distribution",
        default="ADTriage_output/06_distribution_metrics/adt_distribution_metrics.tsv",
        help="Step 6 ADT distribution metrics TSV.",
    )
    parser.add_argument(
        "--rna-correlation",
        default="ADTriage_output/04_signal_consistency/rna_adt_correlation_summary.tsv",
        help="Step 4 RNA-ADT correlation summary TSV.",
    )
    parser.add_argument(
        "--cd45-correlation",
        default="ADTriage_output/04_signal_consistency/cd45isoform_adt_correlation_summary.tsv",
        help="Step 4 CD45isoform-ADT correlation summary TSV.",
    )
    parser.add_argument(
        "--llm-responses",
        default=None,
        help="Optional JSONL responses for LLM-needed rows.",
    )
    parser.add_argument(
        "--outdir",
        default="ADTriage_output/07_final_triage",
        help="Output directory.",
    )
    parser.add_argument(
        "--group-by",
        default="Adult_cHSPCs_Ultima_majority_voting",
        help="Biological metadata column used for Step 5-7 decisions.",
    )
    parser.add_argument(
        "--mapping-confidence-threshold",
        type=float,
        default=0.8,
        help="Rows below this mapping confidence are sent to LLM review.",
    )
    parser.add_argument(
        "--min-rna-spearman",
        type=float,
        default=0.2,
        help="RNA-ADT Spearman correlation below this value is sent to LLM review.",
    )
    parser.add_argument(
        "--llm-confidence-threshold",
        type=float,
        default=0.8,
        help="Accepted LLM response confidence threshold.",
    )
    return parser.parse_args()


def read_tsv(path: str) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: normalize_value(row.get(field, "")) for field in fieldnames})


def normalize_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    return str(value)


def to_float(value: object) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.upper() == "NA":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def to_bool(value: object) -> bool:
    return str(value).strip().upper() in {"TRUE", "T", "1", "YES", "Y"}


def split_groups(value: str | None) -> list[str]:
    if value is None:
        return []
    return [part.strip() for part in str(value).split(";") if part.strip()]


def by_feature(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["feature_id"]: row for row in rows}


def corr_by_feature(rows: list[dict[str, str]], group_by: str) -> dict[str, dict[str, str]]:
    return {row["feature_id"]: row for row in rows if row.get("group_by") == group_by}


def is_control_row(mapping: dict[str, str]) -> bool:
    target_class = mapping.get("target_class", "").strip().lower()
    feature_name = mapping.get("feature_name", "").strip().lower()
    validation_assay = mapping.get("validation_assay", "").strip().lower()
    if target_class in CONTROL_TARGET_CLASSES:
        return True
    control_tokens = ("spike", "isotype", "control", "igg", "igfc", "hto")
    return validation_assay == "none" and any(token in feature_name for token in control_tokens)


def get_correlation_row(
    mapping: dict[str, str],
    rna_corr: dict[str, dict[str, str]],
    cd45_corr: dict[str, dict[str, str]],
) -> dict[str, str]:
    assay = mapping.get("validation_assay", "")
    if assay == "RNA":
        return rna_corr.get(mapping["feature_id"], {})
    if assay == "CD45isoforms":
        return cd45_corr.get(mapping["feature_id"], {})
    return {}


def validation_issue(mapping: dict[str, str], corr_row: dict[str, str]) -> str:
    assay = mapping.get("validation_assay", "")
    if assay not in {"RNA", "CD45isoforms"}:
        return ""
    if not corr_row:
        return f"{assay} validation summary row is missing."
    if corr_row.get("status") != "ok":
        return f"{assay} validation summary status is {corr_row.get('status', 'missing')}."
    return ""


def rna_correlation_issue(mapping: dict[str, str], corr_row: dict[str, str], min_spearman: float) -> str:
    if mapping.get("validation_assay") != "RNA":
        return ""
    if corr_row.get("status") != "ok":
        return ""
    spearman = to_float(corr_row.get("spearman_correlation"))
    if spearman is not None and spearman < min_spearman:
        return f"RNA-ADT Spearman correlation is below threshold: {spearman:.3f} < {min_spearman:.3f}."
    return ""


def make_decision(
    mapping: dict[str, str],
    prior: dict[str, str],
    dist: dict[str, str],
    corr_row: dict[str, str],
    args: argparse.Namespace,
) -> dict[str, object]:
    reasons: list[str] = []
    blocking_reasons: list[str] = []
    candidate_labels: list[str] = []
    final_confidence = ""
    recommended_action = ""
    source = "rule"

    mapping_confidence = to_float(mapping.get("mapping_confidence"))
    step1_needs_llm = to_bool(mapping.get("needs_llm_review"))
    step1_needs_manual = to_bool(mapping.get("needs_manual_review"))
    step5_needs_manual = to_bool(prior.get("needs_manual_review") or dist.get("step5_needs_manual_review"))
    signal_class = dist.get("adt_signal_class", "")
    gmm_quality_status = dist.get("gmm_quality_status", "")
    corr_issue = rna_correlation_issue(mapping, corr_row, args.min_rna_spearman)
    val_issue = validation_issue(mapping, corr_row)

    if is_control_row(mapping):
        return {
            "final_label": "control",
            "final_confidence": "1.0",
            "decision_source": "rule:control",
            "final_reason": "Feature is a control, spike, hashtag, isotype, IgG/IgFc, or validation_assay=none technical feature.",
            "recommended_action": "Keep as technical control; exclude from biological target correctness calls.",
            "llm_review_needed": False,
            "llm_review_reason": "",
            "candidate_labels": "",
        }

    if mapping_confidence is None or mapping_confidence < args.mapping_confidence_threshold:
        reason = (
            f"mapping_confidence is below threshold: {mapping.get('mapping_confidence', '')} < {args.mapping_confidence_threshold}."
        )
        reasons.append(reason)
        blocking_reasons.append(reason)
        candidate_labels.extend(["correct", "suspicious", "wrong"])
    if step1_needs_llm:
        reason = "Step 1 mapping still has needs_llm_review=TRUE."
        reasons.append(reason)
        blocking_reasons.append(reason)
        candidate_labels.extend(["correct", "suspicious", "wrong"])
    if step1_needs_manual:
        reason = "Step 1 mapping has needs_manual_review=TRUE."
        reasons.append(reason)
        blocking_reasons.append(reason)
        candidate_labels.extend(["correct", "suspicious", "wrong"])
    if val_issue:
        reasons.append(val_issue)
        blocking_reasons.append(val_issue)
        candidate_labels.extend(["correct", "suspicious", "wrong"])
    if step5_needs_manual:
        reasons.append("Step 5 prior review marked this row as needing manual review.")
        candidate_labels.extend(["correct", "suspicious", "wrong"])
    if corr_issue:
        reasons.append(corr_issue)
        candidate_labels.extend(["correct", "suspicious"])

    if signal_class in NO_SIGNAL_CLASSES:
        if blocking_reasons:
            candidate_labels.append("no_signal")
            return llm_decision(reasons, candidate_labels)
        return {
            "final_label": "no_signal",
            "final_confidence": "0.95",
            "decision_source": "rule:no_signal",
            "final_reason": f"Step 6 classified ADT signal as {signal_class}.",
            "recommended_action": "Treat as no detectable ADT signal unless later biological context suggests an ultra-rare target.",
            "llm_review_needed": False,
            "llm_review_reason": "",
            "candidate_labels": "",
        }

    if signal_class in TARGET_ENRICHED_CLASSES and not reasons:
        return {
            "final_label": "correct",
            "final_confidence": "0.9",
            "decision_source": "rule:target_enriched",
            "final_reason": f"Mapping is resolved and Step 6 classified ADT signal as {signal_class}.",
            "recommended_action": "Use as target-consistent ADT feature.",
            "llm_review_needed": False,
            "llm_review_reason": "",
            "candidate_labels": "",
        }

    if signal_class == "distribution_detected_but_off_target":
        reasons.append("Step 6 detected an ADT-positive distribution, but it is off-target relative to expected groups.")
        candidate_labels.extend(["suspicious", "wrong"])
    elif signal_class in LLM_SIGNAL_CLASSES:
        reasons.append(f"Step 6 signal class is {signal_class}, which requires semantic review.")
        candidate_labels.extend(["correct", "suspicious", "wrong", "no_signal"])
    elif signal_class in TARGET_ENRICHED_CLASSES and reasons:
        reasons.append(f"Step 6 signal is target-enriched ({signal_class}), but other evidence needs semantic adjudication.")
        candidate_labels.extend(["correct", "suspicious"])
    elif gmm_quality_status in {"fit_failed", "threshold_uncertain", "unimodal_like"}:
        reasons.append(f"GMM quality status is {gmm_quality_status}.")
        candidate_labels.extend(["correct", "suspicious", "no_signal"])
    else:
        reasons.append("No deterministic final rule matched this feature.")
        candidate_labels.extend(["correct", "suspicious", "wrong", "no_signal"])

    return llm_decision(reasons, candidate_labels)


def llm_decision(reasons: list[str], candidate_labels: list[str]) -> dict[str, object]:
    labels = [label for label in dict.fromkeys(candidate_labels) if label in VALID_LABELS and label != "llm_review_needed"]
    return {
        "final_label": "llm_review_needed",
        "final_confidence": "",
        "decision_source": "rule:llm_queue",
        "final_reason": " | ".join(reasons),
        "recommended_action": "Send structured evidence to LLM for final triage label.",
        "llm_review_needed": True,
        "llm_review_reason": " | ".join(reasons),
        "candidate_labels": ";".join(labels),
    }


def build_output_row(
    mapping: dict[str, str],
    prior: dict[str, str],
    dist: dict[str, str],
    corr_row: dict[str, str],
    decision: dict[str, object],
) -> dict[str, object]:
    return {
        "feature_id": mapping.get("feature_id", ""),
        "feature_name": mapping.get("feature_name", ""),
        "target_class": mapping.get("target_class", ""),
        "human_gene_symbol": mapping.get("human_gene_symbol", ""),
        "validation_assay": mapping.get("validation_assay", ""),
        "validation_feature": mapping.get("validation_feature", ""),
        "mapping_method": mapping.get("mapping_method", ""),
        "mapping_confidence": mapping.get("mapping_confidence", ""),
        "mapping_needs_llm_review": mapping.get("needs_llm_review", ""),
        "mapping_needs_manual_review": mapping.get("needs_manual_review", ""),
        "group_by": dist.get("group_by") or prior.get("group_by", ""),
        "expected_positive_groups": dist.get("expected_positive_groups") or prior.get("expected_positive_groups", ""),
        "prior_method": dist.get("prior_method") or prior.get("prior_method", ""),
        "prior_confidence": dist.get("prior_confidence") or prior.get("prior_confidence", ""),
        "observed_adt_positive_groups": dist.get("observed_adt_positive_groups") or prior.get("observed_adt_positive_groups", ""),
        "observed_validation_high_groups": dist.get("observed_validation_high_groups") or prior.get("observed_validation_high_groups", ""),
        "adt_validation_overlap_groups": dist.get("adt_validation_overlap_groups") or prior.get("adt_validation_overlap_groups", ""),
        "expected_adt_overlap_groups": dist.get("expected_adt_overlap_groups") or prior.get("expected_adt_overlap_groups", ""),
        "expected_validation_overlap_groups": dist.get("expected_validation_overlap_groups") or prior.get("expected_validation_overlap_groups", ""),
        "step5_needs_manual_review": dist.get("step5_needs_manual_review") or prior.get("needs_manual_review", ""),
        "step5_review_reason": dist.get("step5_review_reason") or prior.get("review_reason", ""),
        "correlation_status": corr_row.get("status", ""),
        "spearman_correlation": corr_row.get("spearman_correlation", ""),
        "pearson_correlation": corr_row.get("pearson_correlation", ""),
        "linear_model_r_squared": corr_row.get("linear_model_r_squared", ""),
        "n_groups": corr_row.get("n_groups", ""),
        "positive_rate": dist.get("positive_rate", ""),
        "positive_n_cells": dist.get("positive_n_cells", ""),
        "target_positive_rate": dist.get("target_positive_rate", ""),
        "non_target_positive_rate": dist.get("non_target_positive_rate", ""),
        "tnr_raw": dist.get("tnr_raw", ""),
        "target_dropout_pct": dist.get("target_dropout_pct", ""),
        "gmm_quality_status": dist.get("gmm_quality_status", ""),
        "adt_signal_class": dist.get("adt_signal_class", ""),
        "adt_signal_reason": dist.get("adt_signal_reason", ""),
        "final_label": decision["final_label"],
        "final_confidence": decision["final_confidence"],
        "decision_source": decision["decision_source"],
        "final_reason": decision["final_reason"],
        "recommended_action": decision["recommended_action"],
        "llm_review_needed": decision["llm_review_needed"],
        "llm_review_reason": decision["llm_review_reason"],
        "candidate_labels": decision["candidate_labels"],
    }


def make_llm_payload(row: dict[str, object]) -> dict[str, object]:
    return {
        "feature_id": row["feature_id"],
        "feature_name": row["feature_name"],
        "task": "Assign final ADT triage label from structured evidence.",
        "allowed_final_labels": ["correct", "suspicious", "wrong", "no_signal", "control"],
        "candidate_labels": split_groups(str(row.get("candidate_labels", "")).replace(";", ";")),
        "evidence": {
            "mapping": {
                "target_class": row["target_class"],
                "human_gene_symbol": row["human_gene_symbol"],
                "validation_assay": row["validation_assay"],
                "validation_feature": row["validation_feature"],
                "mapping_method": row["mapping_method"],
                "mapping_confidence": row["mapping_confidence"],
            },
            "prior_and_groups": {
                "group_by": row["group_by"],
                "expected_positive_groups": row["expected_positive_groups"],
                "prior_method": row["prior_method"],
                "prior_confidence": row["prior_confidence"],
                "observed_adt_positive_groups": row["observed_adt_positive_groups"],
                "observed_validation_high_groups": row["observed_validation_high_groups"],
                "adt_validation_overlap_groups": row["adt_validation_overlap_groups"],
                "expected_adt_overlap_groups": row["expected_adt_overlap_groups"],
                "expected_validation_overlap_groups": row["expected_validation_overlap_groups"],
            },
            "validation_correlation": {
                "status": row["correlation_status"],
                "spearman_correlation": row["spearman_correlation"],
                "pearson_correlation": row["pearson_correlation"],
                "linear_model_r_squared": row["linear_model_r_squared"],
                "n_groups": row["n_groups"],
            },
            "adt_distribution": {
                "adt_signal_class": row["adt_signal_class"],
                "adt_signal_reason": row["adt_signal_reason"],
                "gmm_quality_status": row["gmm_quality_status"],
                "positive_rate": row["positive_rate"],
                "positive_n_cells": row["positive_n_cells"],
                "target_positive_rate": row["target_positive_rate"],
                "non_target_positive_rate": row["non_target_positive_rate"],
                "tnr_raw": row["tnr_raw"],
                "target_dropout_pct": row["target_dropout_pct"],
            },
            "rule_queue_reason": row["llm_review_reason"],
        },
        "response_schema": {
            "feature_id": "same feature_id",
            "final_label": "one of correct, suspicious, wrong, no_signal, control",
            "final_confidence": "number from 0 to 1",
            "final_reason": "short rationale",
            "recommended_action": "short action for downstream use",
        },
    }


def write_jsonl(path: Path, payloads: list[dict[str, object]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for payload in payloads:
            handle.write(json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n")


def read_llm_responses(path: str) -> dict[str, dict[str, object]]:
    responses: dict[str, dict[str, object]] = {}
    with open(path, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            data = json.loads(line)
            feature_id = str(data.get("feature_id", "")).strip()
            if not feature_id:
                raise ValueError(f"LLM response line {line_no} is missing feature_id.")
            responses[feature_id] = data
    return responses


def apply_llm_responses(
    rows: list[dict[str, object]],
    responses: dict[str, dict[str, object]],
    confidence_threshold: float,
) -> list[dict[str, object]]:
    resolved = []
    for row in rows:
        out = dict(row)
        feature_id = str(row["feature_id"])
        response = responses.get(feature_id)
        if not response or row.get("final_label") != "llm_review_needed":
            resolved.append(out)
            continue
        label = str(response.get("final_label", "")).strip()
        confidence = to_float(response.get("final_confidence"))
        if label not in VALID_LABELS - {"llm_review_needed"}:
            out["final_reason"] = f"{out['final_reason']} | LLM response returned invalid final_label={label}."
            resolved.append(out)
            continue
        if confidence is None or confidence < confidence_threshold:
            out["final_reason"] = (
                f"{out['final_reason']} | LLM confidence below threshold: "
                f"{response.get('final_confidence', '')} < {confidence_threshold}."
            )
            resolved.append(out)
            continue
        out["final_label"] = label
        out["final_confidence"] = str(confidence)
        out["decision_source"] = "llm"
        out["final_reason"] = str(response.get("final_reason", "")).strip()
        out["recommended_action"] = str(response.get("recommended_action", "")).strip()
        out["llm_review_needed"] = False
        out["llm_review_reason"] = ""
        resolved.append(out)
    return resolved


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mapping_rows = read_tsv(args.mapping)
    prior_rows = [
        row for row in read_tsv(args.prior_review)
        if row.get("group_by") == args.group_by
    ]
    dist_rows = [
        row for row in read_tsv(args.distribution)
        if row.get("group_by") == args.group_by
    ]
    rna_corr = corr_by_feature(read_tsv(args.rna_correlation), args.group_by)
    cd45_corr = corr_by_feature(read_tsv(args.cd45_correlation), args.group_by)

    prior_by_id = by_feature(prior_rows)
    dist_by_id = by_feature(dist_rows)

    output_rows: list[dict[str, object]] = []
    for mapping in mapping_rows:
        feature_id = mapping["feature_id"]
        prior = prior_by_id.get(feature_id, {})
        dist = dist_by_id.get(feature_id, {})
        corr_row = get_correlation_row(mapping, rna_corr, cd45_corr)
        decision = make_decision(mapping, prior, dist, corr_row, args)
        output_rows.append(build_output_row(mapping, prior, dist, corr_row, decision))

    fieldnames = [
        "feature_id",
        "feature_name",
        "target_class",
        "human_gene_symbol",
        "validation_assay",
        "validation_feature",
        "mapping_method",
        "mapping_confidence",
        "mapping_needs_llm_review",
        "mapping_needs_manual_review",
        "group_by",
        "expected_positive_groups",
        "prior_method",
        "prior_confidence",
        "observed_adt_positive_groups",
        "observed_validation_high_groups",
        "adt_validation_overlap_groups",
        "expected_adt_overlap_groups",
        "expected_validation_overlap_groups",
        "step5_needs_manual_review",
        "step5_review_reason",
        "correlation_status",
        "spearman_correlation",
        "pearson_correlation",
        "linear_model_r_squared",
        "n_groups",
        "positive_rate",
        "positive_n_cells",
        "target_positive_rate",
        "non_target_positive_rate",
        "tnr_raw",
        "target_dropout_pct",
        "gmm_quality_status",
        "adt_signal_class",
        "adt_signal_reason",
        "final_label",
        "final_confidence",
        "decision_source",
        "final_reason",
        "recommended_action",
        "llm_review_needed",
        "llm_review_reason",
        "candidate_labels",
    ]

    final_path = outdir / "final_adt_triage_table.tsv"
    queue_rows = [row for row in output_rows if row["llm_review_needed"] is True]
    queue_tsv_path = outdir / "llm_review_queue.tsv"
    queue_jsonl_path = outdir / "final_triage.llm.jsonl"

    write_tsv(final_path, output_rows, fieldnames)
    write_tsv(queue_tsv_path, queue_rows, fieldnames)
    write_jsonl(queue_jsonl_path, [make_llm_payload(row) for row in queue_rows])

    print(f"Wrote: {final_path} rows={len(output_rows)}")
    print(f"Wrote: {queue_tsv_path} rows={len(queue_rows)}")
    print(f"Wrote: {queue_jsonl_path} rows={len(queue_rows)}")

    if args.llm_responses:
        responses = read_llm_responses(args.llm_responses)
        resolved_rows = apply_llm_responses(output_rows, responses, args.llm_confidence_threshold)
        resolved_path = outdir / "final_adt_triage_table.llm_resolved.tsv"
        write_tsv(resolved_path, resolved_rows, fieldnames)
        print(f"Wrote: {resolved_path} rows={len(resolved_rows)}")


if __name__ == "__main__":
    main()
