#!/usr/bin/env python3
"""Build ADT feature mapping tables for ADTriage Step 1."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path


FB_REQUIRED_COLUMNS = ["id", "name", "read", "pattern", "sequence", "feature_type"]
MAP_REQUIRED_COLUMNS = [
    "id",
    "genesymbol",
    "genesymbol_source",
    "notes",
    "last_seen_name",
    "last_seen_reference",
    "updated_at",
]

INITIAL_COLUMNS = [
    "feature_id",
    "feature_name",
    "feature_type",
    "genesymbol_from_map",
    "genesymbol_source",
    "last_seen_name",
    "last_seen_reference",
    "target_class",
    "human_gene_symbol",
    "external_target_name",
    "validation_assay",
    "validation_feature",
    "mapping_method",
    "mapping_confidence",
    "mapping_evidence",
    "needs_llm_review",
    "needs_manual_review",
]

REVIEW_COLUMNS = [
    "feature_id",
    "feature_name",
    "target_class",
    "human_gene_symbol",
    "external_target_name",
    "validation_assay",
    "validation_feature",
    "mapping_method",
    "mapping_confidence",
    "llm_suggested_target_class",
    "llm_suggested_gene_symbol",
    "llm_confidence",
    "llm_reason",
    "needs_manual_review",
    "manual_curated_target_class",
    "manual_curated_gene_symbol",
    "manual_notes",
]

CD45_ISOFORM_MAP = {
    "CD45RA": "PTPRC-RA",
    "CD45RB": "PTPRC-RB",
    "CD45RC": "PTPRC-RC",
    "CD45RO": "PTPRC-RO",
}

CONTROL_RE = re.compile(
    r"(^|[^A-Z0-9])(IGG|IGG1|IGFC|ISOTYPE|CONTROL|NEGATIVE CONTROL|BACKGROUND CONTROL)([^A-Z0-9]|$)",
    re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build initial and review-ready ADT feature mapping tables."
    )
    parser.add_argument("--fb-ref", required=True, type=Path, help="Input FB_ref.csv file.")
    parser.add_argument(
        "--gene-map",
        required=True,
        type=Path,
        help="Input mAOC_gene_symbol_map.csv file.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for Step 1 feature mapping tables.",
    )
    parser.add_argument(
        "--llm-review-jsonl",
        type=Path,
        default=None,
        help="Optional JSONL file with LLM suggestions for unresolved/ambiguous features.",
    )
    return parser.parse_args()


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def bool_str(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def read_csv(path: Path, required_columns: list[str], label: str) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{label} has no header: {path}")
        missing = [col for col in required_columns if col not in reader.fieldnames]
        if missing:
            raise ValueError(f"{label} missing required columns: {', '.join(missing)}")
        return [{key: clean(value) for key, value in row.items()} for row in reader]


def build_gene_map(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    by_id: dict[str, dict[str, str]] = {}
    duplicate_ids: set[str] = set()
    for row in rows:
        feature_id = clean(row.get("id"))
        if not feature_id:
            continue
        if feature_id in by_id:
            duplicate_ids.add(feature_id)
            continue
        by_id[feature_id] = row

    if duplicate_ids:
        print(
            f"Warning: ignored duplicate gene-map ids: {', '.join(sorted(duplicate_ids))}",
            file=sys.stderr,
        )
    return by_id


def load_llm_suggestions(path: Path | None) -> dict[str, dict[str, object]]:
    if path is None:
        return {}
    if not path.exists():
        raise FileNotFoundError(f"LLM review JSONL not found: {path}")

    suggestions: dict[str, dict[str, object]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            try:
                record = json.loads(stripped)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSON at {path}:{line_number}: {exc}") from exc
            feature_id = clean(record.get("feature_id"))
            if not feature_id:
                raise ValueError(f"Missing feature_id at {path}:{line_number}")
            suggestions[feature_id] = record
    return suggestions


def joined_text(*values: str) -> str:
    return " | ".join(clean(value) for value in values if clean(value))


def detect_cd45_isoform(text: str) -> tuple[str, str] | None:
    normalized = re.sub(r"[^A-Z0-9]+", "", text.upper())
    for isoform, validation_feature in CD45_ISOFORM_MAP.items():
        if isoform in normalized:
            return isoform, validation_feature
    return None


def is_cd45_total(text: str) -> bool:
    normalized = re.sub(r"[^A-Z0-9]+", "", text.upper())
    if any(isoform in normalized for isoform in CD45_ISOFORM_MAP):
        return False
    return "CD45" in normalized or "PTPRC" in normalized


def is_spike(feature_name: str, genesymbol: str, last_seen_name: str) -> bool:
    text = joined_text(feature_name, genesymbol, last_seen_name).upper()
    compact = re.sub(r"[^A-Z0-9]+", "", text)
    return "SPIKE" in compact or bool(re.search(r"(^|[^A-Z0-9])S[-_]", text))


def is_hashtag(feature_name: str, last_seen_name: str) -> bool:
    text = joined_text(feature_name, last_seen_name).upper()
    return bool(re.search(r"(^|[^A-Z0-9])HTO[-_A-Z0-9]*", text)) or "HASHTAG" in text


def is_control(feature_name: str, genesymbol: str, last_seen_name: str, notes: str) -> bool:
    text = joined_text(feature_name, genesymbol, last_seen_name, notes)
    return bool(CONTROL_RE.search(text))


def is_multi_gene_symbol(genesymbol: str) -> bool:
    return "/" in genesymbol


def classify_feature(fb_row: dict[str, str], map_row: dict[str, str] | None) -> dict[str, str]:
    feature_id = clean(fb_row.get("id"))
    feature_name = clean(fb_row.get("name"))
    feature_type = clean(fb_row.get("feature_type"))

    map_row = map_row or {}
    genesymbol = clean(map_row.get("genesymbol"))
    genesymbol_source = clean(map_row.get("genesymbol_source"))
    notes = clean(map_row.get("notes"))
    last_seen_name = clean(map_row.get("last_seen_name"))
    last_seen_reference = clean(map_row.get("last_seen_reference"))
    text_for_rules = joined_text(feature_name, genesymbol, last_seen_name, notes)

    base = {
        "feature_id": feature_id,
        "feature_name": feature_name,
        "feature_type": feature_type,
        "genesymbol_from_map": genesymbol,
        "genesymbol_source": genesymbol_source,
        "last_seen_name": last_seen_name,
        "last_seen_reference": last_seen_reference,
        "target_class": "",
        "human_gene_symbol": "",
        "external_target_name": "",
        "validation_assay": "",
        "validation_feature": "",
        "mapping_method": "",
        "mapping_confidence": "",
        "mapping_evidence": "",
        "needs_llm_review": "FALSE",
        "needs_manual_review": "FALSE",
    }

    isoform = detect_cd45_isoform(text_for_rules)
    if isoform is not None:
        isoform_name, validation_feature = isoform
        base.update(
            {
                "target_class": "cd45_isoform",
                "human_gene_symbol": "PTPRC",
                "validation_assay": "CD45isoforms",
                "validation_feature": validation_feature,
                "mapping_method": "rule:cd45_isoform",
                "mapping_confidence": "1.00",
                "mapping_evidence": f"{isoform_name} detected in feature/map fields",
            }
        )
        return base

    if is_cd45_total(text_for_rules):
        base.update(
            {
                "target_class": "gene",
                "human_gene_symbol": "PTPRC",
                "validation_assay": "RNA",
                "validation_feature": "PTPRC",
                "mapping_method": "rule:cd45_total",
                "mapping_confidence": "1.00",
                "mapping_evidence": "CD45/PTPRC total expression detected",
            }
        )
        return base

    if is_spike(feature_name, genesymbol, last_seen_name):
        external_name = genesymbol or last_seen_name or feature_name
        base.update(
            {
                "target_class": "spike_control",
                "external_target_name": external_name,
                "validation_assay": "none",
                "mapping_method": "rule:spike_control",
                "mapping_confidence": "1.00",
                "mapping_evidence": "Spike target detected",
            }
        )
        return base

    if is_hashtag(feature_name, last_seen_name):
        base.update(
            {
                "target_class": "hashtag",
                "external_target_name": feature_name,
                "validation_assay": "HTO",
                "validation_feature": feature_name,
                "mapping_method": "rule:hashtag",
                "mapping_confidence": "0.95",
                "mapping_evidence": "HTO hashtag feature detected",
            }
        )
        return base

    if is_control(feature_name, genesymbol, last_seen_name, notes):
        external_name = last_seen_name or feature_name or genesymbol
        base.update(
            {
                "target_class": "control",
                "external_target_name": external_name,
                "validation_assay": "none",
                "mapping_method": "rule:control",
                "mapping_confidence": "1.00",
                "mapping_evidence": "Control/isotype/IgG feature detected",
            }
        )
        return base

    if genesymbol:
        needs_llm_review = is_multi_gene_symbol(genesymbol)
        base.update(
            {
                "target_class": "gene",
                "human_gene_symbol": genesymbol,
                "validation_assay": "RNA",
                "validation_feature": genesymbol,
                "mapping_method": "map:id_left_join",
                "mapping_confidence": "0.70" if needs_llm_review else "0.95",
                "mapping_evidence": (
                    "multi-gene symbol found by id left join; requires LLM/manual selection of one RNA validation feature"
                    if needs_llm_review
                    else "genesymbol found by FB_ref.id to mAOC_gene_symbol_map.id left join"
                ),
                "needs_llm_review": bool_str(needs_llm_review),
                "needs_manual_review": bool_str(needs_llm_review),
            }
        )
        return base

    base.update(
        {
            "target_class": "ambiguous",
            "validation_assay": "manual",
            "mapping_method": "unresolved",
            "mapping_confidence": "0.00",
            "mapping_evidence": "No map hit or rule-based target assignment",
            "needs_llm_review": "TRUE",
            "needs_manual_review": "TRUE",
        }
    )
    return base


def apply_llm_suggestion(
    initial_row: dict[str, str], suggestion: dict[str, object] | None
) -> dict[str, str]:
    llm_confidence = ""
    if suggestion is not None and suggestion.get("confidence") is not None:
        llm_confidence = str(suggestion.get("confidence"))

    needs_manual_review = initial_row["needs_manual_review"] == "TRUE"
    try:
        if llm_confidence and float(llm_confidence) < 0.8:
            needs_manual_review = True
    except ValueError:
        needs_manual_review = True

    if initial_row["mapping_confidence"]:
        try:
            if float(initial_row["mapping_confidence"]) < 0.8:
                needs_manual_review = True
        except ValueError:
            needs_manual_review = True

    return {
        "feature_id": initial_row["feature_id"],
        "feature_name": initial_row["feature_name"],
        "target_class": initial_row["target_class"],
        "human_gene_symbol": initial_row["human_gene_symbol"],
        "external_target_name": initial_row["external_target_name"],
        "validation_assay": initial_row["validation_assay"],
        "validation_feature": initial_row["validation_feature"],
        "mapping_method": initial_row["mapping_method"],
        "mapping_confidence": initial_row["mapping_confidence"],
        "llm_suggested_target_class": clean(suggestion.get("target_class")) if suggestion else "",
        "llm_suggested_gene_symbol": clean(suggestion.get("human_gene_symbol")) if suggestion else "",
        "llm_confidence": llm_confidence,
        "llm_reason": clean(suggestion.get("reason")) if suggestion else "",
        "needs_manual_review": bool_str(needs_manual_review),
        "manual_curated_target_class": "",
        "manual_curated_gene_symbol": "",
        "manual_notes": "",
    }


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
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


def main() -> int:
    args = parse_args()
    try:
        fb_rows = read_csv(args.fb_ref, FB_REQUIRED_COLUMNS, "FB_ref.csv")
        map_rows = read_csv(args.gene_map, MAP_REQUIRED_COLUMNS, "mAOC_gene_symbol_map.csv")
        gene_map = build_gene_map(map_rows)
        llm_suggestions = load_llm_suggestions(args.llm_review_jsonl)

        args.outdir.mkdir(parents=True, exist_ok=True)
        initial_rows = [
            classify_feature(fb_row, gene_map.get(clean(fb_row.get("id")))) for fb_row in fb_rows
        ]
        review_rows = [
            apply_llm_suggestion(row, llm_suggestions.get(row["feature_id"])) for row in initial_rows
        ]

        initial_path = args.outdir / "adt_feature_mapping.initial.tsv"
        review_path = args.outdir / "adt_feature_mapping.review_ready.tsv"
        write_tsv(initial_path, initial_rows, INITIAL_COLUMNS)
        write_tsv(review_path, review_rows, REVIEW_COLUMNS)

        print(f"Wrote {len(initial_rows)} rows: {initial_path}")
        print(f"Wrote {len(review_rows)} rows: {review_path}")
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
