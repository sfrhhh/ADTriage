#!/usr/bin/env python3
"""Build marker-prior review tables for ADTriage Step 5."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path


EXCLUDED_TARGET_CLASSES = {"control", "spike_control", "hashtag"}
KNOWN_ADULT_GROUP_LEVELS = {
    "B",
    "CLP",
    "ERYP",
    "GMP-L",
    "HSC_MPP",
    "MEBEMP-L",
    "MKP",
    "Monocyte",
    "NK/T",
}
DEFAULT_PRIOR_FILE = Path("priors/adult_cHSPC_marker_prior.strict.tsv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Step 5 marker prior review tables.")
    parser.add_argument("--mapping", required=True, type=Path, help="Step 1 llm_review mapping TSV.")
    parser.add_argument(
        "--group-summary-dir",
        required=True,
        type=Path,
        help="Directory containing Step 3 group summary TSVs.",
    )
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    parser.add_argument(
        "--prior-file",
        default=DEFAULT_PRIOR_FILE,
        type=Path,
        help=f"Marker prior TSV. Default: {DEFAULT_PRIOR_FILE}",
    )
    parser.add_argument(
        "--group-by",
        default="",
        help="Optional comma-separated group_by columns to include in prior review.",
    )
    parser.add_argument(
        "--top-quantile",
        type=float,
        default=0.75,
        help="Quantile used to define observed high groups. Default: 0.75.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Input not found: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Input has no header: {path}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_jsonl(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")


def parse_float(value: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def quantile(values: list[float], q: float) -> float:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return math.nan
    if len(clean) == 1:
        return clean[0]
    position = (len(clean) - 1) * q
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return clean[int(position)]
    weight = position - lower
    return clean[lower] * (1 - weight) + clean[upper] * weight


def split_groups(value: str) -> set[str]:
    if not value:
        return set()
    return {part.strip() for part in value.split(";") if part.strip()}


def parse_group_by_filter(value: str) -> set[str]:
    if not value:
        return set()
    return {part.strip() for part in value.split(",") if part.strip()}


def join_groups(groups: set[str] | list[str]) -> str:
    return ";".join(sorted(groups))


def compact_symbol(value: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "", value.upper())


def group_levels(rows: list[dict[str, str]], group_by: str) -> list[str]:
    seen: list[str] = []
    for row in rows:
        if row["group_by"] != group_by:
            continue
        level = row["group_level"]
        if level not in seen:
            seen.append(level)
    return seen


def load_marker_prior(path: Path) -> dict[str, dict[str, object]]:
    required = {
        "marker_symbol",
        "expected_positive_groups",
        "prior_confidence",
        "prior_action",
        "prior_evidence",
    }
    rows = read_tsv(path)
    if not rows:
        raise ValueError(f"Prior file is empty: {path}")
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"Prior file missing columns: {', '.join(sorted(missing))}")

    prior: dict[str, dict[str, object]] = {}
    for row in rows:
        symbol = row["marker_symbol"]
        if not symbol:
            raise ValueError(f"Prior file contains empty marker_symbol: {path}")
        action = row["prior_action"] or "vote"
        if action not in {"vote", "no_vote"}:
            raise ValueError(f"Invalid prior_action for {symbol}: {action}")
        groups = split_groups(row["expected_positive_groups"])
        if action == "no_vote" and groups:
            raise ValueError(f"no_vote marker must not have expected groups: {symbol}")
        confidence = parse_float(row["prior_confidence"])
        if not math.isfinite(confidence):
            raise ValueError(f"Invalid prior_confidence for {symbol}: {row['prior_confidence']}")
        prior[compact_symbol(symbol)] = {
            "groups": groups,
            "confidence": confidence,
            "action": action,
            "evidence": row["prior_evidence"],
        }
    return prior


def observed_high_groups(
    rows: list[dict[str, str]],
    group_by: str,
    feature_id: str,
    value_col: str,
    top_quantile: float,
) -> set[str]:
    subset = [
        row
        for row in rows
        if row["group_by"] == group_by
        and row["feature_id"] == feature_id
        and row.get("status") == "ok"
    ]
    values = [parse_float(row.get(value_col, "")) for row in subset]
    finite_values = [value for value in values if math.isfinite(value)]
    if not finite_values:
        return set()
    cutoff = quantile(finite_values, top_quantile)
    return {
        row["group_level"]
        for row in subset
        if math.isfinite(parse_float(row.get(value_col, "")))
        and parse_float(row.get(value_col, "")) >= cutoff
        and parse_float(row.get(value_col, "")) > 0
    }


def table_by_feature(rows: list[dict[str, str]]) -> dict[tuple[str, str], list[dict[str, str]]]:
    result: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        result.setdefault((row["group_by"], row["feature_id"]), []).append(row)
    return result


def expected_groups_for(
    row: dict[str, str],
    levels: list[str],
    prior: dict[str, dict[str, object]],
) -> tuple[set[str], str, float, str, bool]:
    level_set = set(levels)
    if not level_set:
        return set(), "none", 0.0, "No group levels available.", True
    if not level_set.issubset(KNOWN_ADULT_GROUP_LEVELS):
        return (
            set(),
            "manual:uninterpretable_group_levels",
            0.0,
            "Group levels are not biological labels in the local marker prior.",
            True,
        )

    candidates = [
        row.get("validation_feature", ""),
        row.get("human_gene_symbol", ""),
        row.get("feature_name", ""),
    ]
    for candidate in candidates:
        key = compact_symbol(candidate)
        if key in prior:
            entry = prior[key]
            groups = set(entry["groups"])
            evidence = str(entry["evidence"])
            confidence = float(entry["confidence"])
            if entry["action"] == "no_vote":
                return set(), "local_marker_prior:no_vote", confidence, evidence, False
            return groups & level_set, "local_marker_prior", confidence, evidence, False

    return (
        set(),
        "llm_needed:no_local_prior",
        0.0,
        "No local marker prior found for this feature and biological group labels.",
        True,
    )


def build_review_rows(
    mapping: list[dict[str, str]],
    adt_rows: list[dict[str, str]],
    rna_rows: list[dict[str, str]],
    cd45_rows: list[dict[str, str]],
    top_quantile: float,
    group_by_filter: set[str],
    prior: dict[str, dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    group_bys = sorted({row["group_by"] for row in adt_rows})
    if group_by_filter:
        group_bys = [group_by for group_by in group_bys if group_by in group_by_filter]
    if not group_bys:
        raise ValueError("No group_by columns selected for prior review.")
    rna_by_feature = table_by_feature(rna_rows)
    cd45_by_feature = table_by_feature(cd45_rows)

    rows: list[dict[str, object]] = []
    llm_queue: list[dict[str, object]] = []
    features = [row for row in mapping if row.get("target_class", "") not in EXCLUDED_TARGET_CLASSES]

    for feature in features:
        for group_by in group_bys:
            levels = group_levels(adt_rows, group_by)
            expected, method, confidence, evidence, needs_llm = expected_groups_for(
                feature, levels, prior
            )
            adt_high = observed_high_groups(
                adt_rows, group_by, feature["feature_id"], "mean_CLR", top_quantile
            )

            validation_assay = feature.get("validation_assay", "")
            validation_high: set[str]
            if validation_assay == "RNA":
                validation_high = observed_high_groups(
                    rna_by_feature.get((group_by, feature["feature_id"]), []),
                    group_by,
                    feature["feature_id"],
                    "mean_RNA_lognormalized",
                    top_quantile,
                )
            elif validation_assay == "CD45isoforms":
                validation_high = observed_high_groups(
                    cd45_by_feature.get((group_by, feature["feature_id"]), []),
                    group_by,
                    feature["feature_id"],
                    "mean_CD45isoform_signal",
                    top_quantile,
                )
            else:
                validation_high = set()

            expected_adt_overlap = expected & adt_high
            expected_validation_overlap = expected & validation_high
            review_reason = ""
            needs_manual = confidence < 0.8
            if needs_llm:
                review_reason = evidence
            elif expected and not expected_adt_overlap:
                review_reason = "No observed ADT-high group overlaps expected prior groups."
                needs_manual = True
            elif expected and validation_high and not expected_validation_overlap:
                review_reason = "No observed validation-high group overlaps expected prior groups."
                needs_manual = True

            output_row: dict[str, object] = {
                "feature_id": feature["feature_id"],
                "feature_name": feature["feature_name"],
                "target_class": feature["target_class"],
                "human_gene_symbol": feature.get("human_gene_symbol", ""),
                "validation_assay": validation_assay,
                "validation_feature": feature.get("validation_feature", ""),
                "group_by": group_by,
                "group_levels": join_groups(set(levels)),
                "expected_positive_groups": join_groups(expected),
                "prior_method": method,
                "prior_confidence": f"{confidence:.2f}",
                "prior_evidence": evidence,
                "observed_adt_positive_groups": join_groups(adt_high),
                "observed_validation_high_groups": join_groups(validation_high),
                "adt_validation_overlap_groups": join_groups(adt_high & validation_high),
                "expected_adt_overlap_groups": join_groups(expected_adt_overlap),
                "expected_validation_overlap_groups": join_groups(expected_validation_overlap),
                "n_expected_groups": len(expected),
                "n_observed_adt_groups": len(adt_high),
                "n_observed_validation_groups": len(validation_high),
                "needs_llm_review": "TRUE" if needs_llm else "FALSE",
                "needs_manual_review": "TRUE" if needs_manual else "FALSE",
                "review_reason": review_reason,
            }
            rows.append(output_row)

            if needs_llm:
                llm_queue.append(
                    {
                        "feature_id": feature["feature_id"],
                        "feature_name": feature["feature_name"],
                        "target_class": feature["target_class"],
                        "human_gene_symbol": feature.get("human_gene_symbol", ""),
                        "validation_assay": validation_assay,
                        "validation_feature": feature.get("validation_feature", ""),
                        "group_by": group_by,
                        "group_levels": levels,
                        "observed_adt_positive_groups": sorted(adt_high),
                        "observed_validation_high_groups": sorted(validation_high),
                        "reason": evidence,
                        "required_output_schema": {
                            "expected_positive_groups": "list[str]",
                            "confidence": "float between 0 and 1",
                            "reason": "short rationale",
                        },
                    }
                )

    return rows, llm_queue


def main() -> int:
    args = parse_args()
    try:
        mapping = read_tsv(args.mapping)
        adt_rows = read_tsv(args.group_summary_dir / "adt_group_summary.tsv")
        rna_rows = read_tsv(args.group_summary_dir / "rna_group_summary.tsv")
        cd45_rows = read_tsv(args.group_summary_dir / "cd45isoform_group_summary.tsv")
        marker_prior = load_marker_prior(args.prior_file)

        review_rows, llm_queue = build_review_rows(
            mapping,
            adt_rows,
            rna_rows,
            cd45_rows,
            args.top_quantile,
            parse_group_by_filter(args.group_by),
            marker_prior,
        )

        columns = [
            "feature_id",
            "feature_name",
            "target_class",
            "human_gene_symbol",
            "validation_assay",
            "validation_feature",
            "group_by",
            "group_levels",
            "expected_positive_groups",
            "prior_method",
            "prior_confidence",
            "prior_evidence",
            "observed_adt_positive_groups",
            "observed_validation_high_groups",
            "adt_validation_overlap_groups",
            "expected_adt_overlap_groups",
            "expected_validation_overlap_groups",
            "n_expected_groups",
            "n_observed_adt_groups",
            "n_observed_validation_groups",
            "needs_llm_review",
            "needs_manual_review",
            "review_reason",
        ]
        review_path = args.outdir / "marker_prior_review.tsv"
        llm_path = args.outdir / "marker_prior_review.llm.jsonl"
        write_tsv(review_path, review_rows, columns)
        write_jsonl(llm_path, llm_queue)

        print(f"Wrote {len(review_rows)} rows: {review_path}")
        print(f"Wrote {len(llm_queue)} LLM review prompts: {llm_path}")
    except Exception as exc:
        print(f"Error: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
