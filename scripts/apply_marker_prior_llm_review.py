#!/usr/bin/env python3
"""Merge Step 5 LLM marker-prior responses back into the review table."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Apply Step 5 LLM prior review responses.")
    parser.add_argument("--review-table", required=True, type=Path, help="Input marker_prior_review TSV.")
    parser.add_argument("--llm-responses", required=True, type=Path, help="LLM response JSONL.")
    parser.add_argument("--out", required=True, type=Path, help="Merged output TSV.")
    parser.add_argument(
        "--confidence-threshold",
        type=float,
        default=0.8,
        help="LLM confidence below this threshold remains manual review. Default: 0.8.",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Allow missing responses for rows where needs_llm_review is TRUE.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(f"Input not found: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Input has no header: {path}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_jsonl(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        raise FileNotFoundError(f"Input not found: {path}")
    rows: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSON at {path}:{line_number}: {exc}") from exc
    return rows


def split_groups(value: str) -> set[str]:
    if not value:
        return set()
    return {part.strip() for part in value.split(";") if part.strip()}


def join_groups(groups: set[str] | list[str]) -> str:
    return ";".join(sorted(groups))


def parse_confidence(value: object, feature_id: str) -> float:
    try:
        confidence = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid confidence for {feature_id}: {value}") from exc
    if not math.isfinite(confidence) or confidence < 0 or confidence > 1:
        raise ValueError(f"Confidence must be between 0 and 1 for {feature_id}: {value}")
    return confidence


def response_key(row: dict[str, object]) -> tuple[str, str]:
    feature_id = str(row.get("feature_id", "")).strip()
    group_by = str(row.get("group_by", "")).strip()
    if not feature_id or not group_by:
        raise ValueError(f"Response missing feature_id or group_by: {row}")
    return feature_id, group_by


def build_response_index(rows: list[dict[str, object]]) -> dict[tuple[str, str], dict[str, object]]:
    result: dict[tuple[str, str], dict[str, object]] = {}
    for row in rows:
        key = response_key(row)
        if key in result:
            raise ValueError(f"Duplicate LLM response for {key[0]} / {key[1]}")
        groups = row.get("expected_positive_groups")
        if not isinstance(groups, list) or not all(isinstance(group, str) for group in groups):
            raise ValueError(f"expected_positive_groups must be list[str] for {key[0]}")
        parse_confidence(row.get("confidence"), key[0])
        if not str(row.get("reason", "")).strip():
            raise ValueError(f"Response missing reason for {key[0]}")
        result[key] = row
    return result


def apply_response(
    review_row: dict[str, str],
    response: dict[str, object],
    confidence_threshold: float,
) -> dict[str, str]:
    row = dict(review_row)
    group_levels = split_groups(row["group_levels"])
    expected = set(str(group).strip() for group in response["expected_positive_groups"] if str(group).strip())
    unknown = expected - group_levels
    if unknown:
        raise ValueError(
            f"LLM response for {row['feature_id']} contains groups not present in group_levels: "
            f"{join_groups(unknown)}"
        )

    confidence = parse_confidence(response["confidence"], row["feature_id"])
    reason = str(response["reason"]).strip()
    adt_high = split_groups(row["observed_adt_positive_groups"])
    validation_high = split_groups(row["observed_validation_high_groups"])
    expected_adt_overlap = expected & adt_high
    expected_validation_overlap = expected & validation_high

    row["expected_positive_groups"] = join_groups(expected)
    row["prior_method"] = "llm_marker_prior" if expected else "llm_marker_prior:no_vote"
    row["prior_confidence"] = f"{confidence:.2f}"
    row["prior_evidence"] = reason
    row["expected_adt_overlap_groups"] = join_groups(expected_adt_overlap)
    row["expected_validation_overlap_groups"] = join_groups(expected_validation_overlap)
    row["n_expected_groups"] = str(len(expected))
    row["needs_llm_review"] = "FALSE"

    needs_manual = confidence < confidence_threshold
    review_reason = ""
    if needs_manual:
        review_reason = "LLM confidence below threshold."
    elif expected and not expected_adt_overlap:
        needs_manual = True
        review_reason = "No observed ADT-high group overlaps LLM expected groups."
    elif expected and validation_high and not expected_validation_overlap:
        needs_manual = True
        review_reason = "No observed validation-high group overlaps LLM expected groups."

    row["needs_manual_review"] = "TRUE" if needs_manual else "FALSE"
    row["review_reason"] = review_reason
    return row


def main() -> int:
    args = parse_args()
    try:
        fieldnames, review_rows = read_tsv(args.review_table)
        responses = build_response_index(read_jsonl(args.llm_responses))
        merged_rows: list[dict[str, str]] = []
        missing: list[tuple[str, str]] = []

        for row in review_rows:
            key = (row["feature_id"], row["group_by"])
            if row.get("needs_llm_review") == "TRUE":
                response = responses.get(key)
                if response is None:
                    missing.append(key)
                    merged_rows.append(row)
                else:
                    merged_rows.append(apply_response(row, response, args.confidence_threshold))
            else:
                merged_rows.append(row)

        if missing and not args.allow_missing:
            missing_text = ", ".join(f"{feature_id}/{group_by}" for feature_id, group_by in missing)
            raise ValueError(f"Missing LLM responses: {missing_text}")

        write_tsv(args.out, fieldnames, merged_rows)
        print(f"Wrote {len(merged_rows)} rows: {args.out}")
        print(f"Applied {len(responses)} LLM responses")
        if missing:
            print(f"Missing responses left unchanged: {len(missing)}")
    except Exception as exc:
        print(f"Error: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
