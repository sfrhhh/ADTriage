#!/usr/bin/env python3
"""Apply LLM-assisted review to low-confidence ADT feature mappings."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


THRESHOLD_DEFAULT = 0.8

OUTPUT_COLUMNS = [
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
    "llm_model",
    "llm_suggested_target_class",
    "llm_suggested_gene_symbol",
    "llm_suggested_validation_assay",
    "llm_suggested_validation_feature",
    "llm_confidence",
    "llm_reason",
]

# This table records this run's LLM-assisted judgments for low-confidence rows.
# It is intentionally explicit so the result is reproducible and reviewable.
LLM_SUGGESTIONS = {
    "mAOC0301": ("gene", "FCGR3A", "RNA", "FCGR3A", 0.86, "CD16 maps to FCGR3A/FCGR3B in the reference map; select FCGR3A as the single RNA validation feature for this workflow."),
    "mAOC0453": ("gene", "HLA-A", "RNA", "HLA-A", 0.86, "HLA-ABC maps to HLA-A/HLA-B/HLA-C in the reference map; select HLA-A as the single RNA validation feature for this workflow."),
    "mAOC0518-1": ("gene", "ACVR2B", "RNA", "ACVR2B", 0.96, "Feature name directly names the human gene ACVR2B."),
    "mAOC0519-1": ("gene", "ALK", "RNA", "ALK", 0.90, "CD246 is commonly used for ALK."),
    "mAOC0522-1": ("gene", "FZD10", "RNA", "FZD10", 0.92, "CD350 corresponds to FZD10."),
    "mAOC0523-1": ("gene", "GLP1R", "RNA", "GLP1R", 0.96, "Feature name directly names GLP1R."),
    "mAOC363": ("gene", "CD47", "RNA", "CD47", 0.96, "CD47 is already the human gene symbol."),
    "mAOC365": ("gene", "ITGB1", "RNA", "ITGB1", 0.93, "CD29 corresponds to integrin beta 1, ITGB1."),
    "mAOC366": ("gene", "CD244", "RNA", "CD244", 0.95, "CD244 is already the human gene symbol."),
    "mAOC373": ("gene", "CD81", "RNA", "CD81", 0.96, "CD81 is already the human gene symbol."),
    "mAOC374": ("gene", "FCGR2A", "RNA", "FCGR2A", 0.72, "CD32 can indicate FCGR2A, FCGR2B, or FCGR2C; project map often uses FCGR2A for unsuffixed CD32, but this remains ambiguous."),
    "mAOC375": ("gene", "SELPLG", "RNA", "SELPLG", 0.93, "CD162 corresponds to PSGL-1, encoded by SELPLG."),
    "mAOC376": ("gene", "TNFRSF18", "RNA", "TNFRSF18", 0.93, "CD357 corresponds to GITR, encoded by TNFRSF18."),
    "mAOC377": ("gene", "EPCAM", "RNA", "EPCAM", 0.95, "CD326 corresponds to EPCAM."),
    "mAOC378": ("gene", "TGFBI", "RNA", "TGFBI", 0.96, "Feature name directly names TGFBI."),
    "mAOC379": ("gene", "CD68", "RNA", "CD68", 0.96, "CD68 is already the human gene symbol."),
    "mAOC381": ("gene", "ICAM2", "RNA", "ICAM2", 0.92, "CD102 corresponds to ICAM2."),
    "mAOC382": ("gene", "TIGIT", "RNA", "TIGIT", 0.96, "Feature name directly names TIGIT."),
    "mAOC383": ("gene", "MCAM", "RNA", "MCAM", 0.94, "CD146 corresponds to MCAM."),
    "mAOC384": ("gene", "KIR2DL3", "RNA", "KIR2DL3", 0.74, "CD158b2 is commonly associated with KIR2DL3, but CD158 family naming is subtype-sensitive."),
    "mAOC386": ("gene", "LILRB1", "RNA", "LILRB1", 0.93, "CD85j corresponds to LILRB1."),
    "mAOC387": ("gene", "IFNGR1", "RNA", "IFNGR1", 0.93, "CD119 corresponds to IFNGR1."),
    "mAOC388": ("gene", "LAIR1", "RNA", "LAIR1", 0.94, "CD305 corresponds to LAIR1."),
    "mAOC391": ("gene", "NRP1", "RNA", "NRP1", 0.94, "CD304 corresponds to NRP1."),
    "mAOC392": ("gene", "SLAMF6", "RNA", "SLAMF6", 0.93, "CD352 corresponds to SLAMF6."),
    "mAOC411": ("gene", "IL2RB", "RNA", "IL2RB", 0.93, "CD122 corresponds to IL2RB."),
    "mAOC413": ("gene", "CSF2RA", "RNA", "CSF2RA", 0.93, "CD116 corresponds to CSF2RA."),
    "mAOC415": ("gene", "KIR2DL2", "RNA", "KIR2DL2", 0.74, "CD158b1 is commonly associated with KIR2DL2, but CD158 family naming is subtype-sensitive."),
    "mAOC416": ("gene", "ICOS", "RNA", "ICOS", 0.94, "CD278 corresponds to ICOS."),
    "mAOC417": ("gene", "ITGB2", "RNA", "ITGB2", 0.93, "CD18 corresponds to ITGB2."),
    "mAOC419": ("gene", "ITGAL", "RNA", "ITGAL", 0.93, "CD11a corresponds to ITGAL."),
    "mAOC421": ("gene", "CD37", "RNA", "CD37", 0.96, "CD37 is already the human gene symbol."),
    "mAOC422": ("gene", "DPP4", "RNA", "DPP4", 0.94, "CD26 corresponds to DPP4."),
    "mAOC423": ("gene", "ENTPD1", "RNA", "ENTPD1", 0.94, "CD39 corresponds to ENTPD1."),
    "mAOC426": ("gene", "SIRPA", "RNA", "SIRPA", 0.93, "CD172a corresponds to SIRPA."),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create an LLM-reviewed mapping table from low-confidence initial mappings."
    )
    parser.add_argument(
        "--initial",
        required=True,
        type=Path,
        help="Input adt_feature_mapping.initial.tsv file.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output adt_feature_mapping.llm_review.tsv file.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=THRESHOLD_DEFAULT,
        help="Rows with mapping_confidence below this value receive LLM review.",
    )
    return parser.parse_args()


def parse_confidence(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def bool_str(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Input table not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Input table has no header: {path}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def apply_review(row: dict[str, str], threshold: float) -> dict[str, str]:
    reviewed = dict(row)
    for column in OUTPUT_COLUMNS:
        reviewed.setdefault(column, "")

    feature_id = row["feature_id"]
    suggestion = LLM_SUGGESTIONS.get(feature_id)

    confidence = parse_confidence(row.get("mapping_confidence", ""))
    if confidence >= threshold and suggestion is None:
        reviewed["llm_model"] = ""
        reviewed["llm_suggested_target_class"] = ""
        reviewed["llm_suggested_gene_symbol"] = ""
        reviewed["llm_suggested_validation_assay"] = ""
        reviewed["llm_suggested_validation_feature"] = ""
        reviewed["llm_confidence"] = ""
        reviewed["llm_reason"] = ""
        return reviewed

    if suggestion is None:
        reviewed["llm_model"] = "codex_manual_review"
        reviewed["llm_confidence"] = "0.00"
        reviewed["llm_reason"] = "No LLM suggestion available for this low-confidence feature."
        reviewed["needs_llm_review"] = "TRUE"
        reviewed["needs_manual_review"] = "TRUE"
        return reviewed

    target_class, gene_symbol, validation_assay, validation_feature, llm_confidence, reason = suggestion
    llm_confidence_str = f"{llm_confidence:.2f}"
    needs_manual_review = llm_confidence < threshold

    reviewed.update(
        {
            "target_class": target_class,
            "human_gene_symbol": gene_symbol,
            "external_target_name": "",
            "validation_assay": validation_assay,
            "validation_feature": validation_feature,
            "mapping_method": "llm:name_to_gene",
            "mapping_confidence": llm_confidence_str,
            "mapping_evidence": reason,
            "needs_llm_review": "FALSE",
            "needs_manual_review": bool_str(needs_manual_review),
            "llm_model": "codex_manual_review",
            "llm_suggested_target_class": target_class,
            "llm_suggested_gene_symbol": gene_symbol,
            "llm_suggested_validation_assay": validation_assay,
            "llm_suggested_validation_feature": validation_feature,
            "llm_confidence": llm_confidence_str,
            "llm_reason": reason,
        }
    )
    return reviewed


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in OUTPUT_COLUMNS})


def main() -> int:
    args = parse_args()
    try:
        rows = read_tsv(args.initial)
        reviewed_rows = [apply_review(row, args.threshold) for row in rows]
        write_tsv(args.output, reviewed_rows)

        reviewed_count = sum(1 for row in reviewed_rows if row["llm_model"])
        manual_count = sum(1 for row in reviewed_rows if row["needs_manual_review"] == "TRUE")
        print(f"Wrote {len(reviewed_rows)} rows: {args.output}")
        print(f"LLM-reviewed rows: {reviewed_count}")
        print(f"Rows still needing manual review: {manual_count}")
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
