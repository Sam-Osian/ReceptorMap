#!/usr/bin/env python
"""Join PDSP Ki data with GtoPdb targets/ligands/interactions and emit a CSV plus diagnostics."""
from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd


DATA_DIR = Path("data")
REPORT_DIR = Path("reports")
REPORT_DIR.mkdir(exist_ok=True)

logger = logging.getLogger("join_data")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)


def read_gtp_csv(path: Path) -> pd.DataFrame:
    """Read GtoPdb CSVs that include a metadata row as the first line."""
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        first_line = handle.readline().strip()
        if not first_line.startswith('"# GtoPdb Version'):
            logger.warning("%s: missing or unexpected metadata row: %s", path.name, first_line)
    df = pd.read_csv(path, skiprows=1, low_memory=False)
    df.columns = [col.strip() for col in df.columns]
    return df


def read_pdsp_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    df.columns = [col.strip() for col in df.columns]
    return df


def require_columns(df: pd.DataFrame, required: Iterable[str], name: str) -> None:
    missing = [col for col in required if col not in df.columns]
    if missing:
        logger.error("%s missing required columns: %s", name, ", ".join(missing))
        raise ValueError(f"{name} missing required columns: {', '.join(missing)}")


def normalise_text(value: object) -> Optional[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    return str(value).strip().lower()


def normalise_symbol(value: object) -> Optional[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    return str(value).strip().upper()


def build_lookup_map(values: pd.Series, ligand_ids: pd.Series) -> Dict[str, List[str]]:
    mapping: Dict[str, List[str]] = {}
    for value, ligand_id in zip(values, ligand_ids):
        if value is None or (isinstance(value, float) and pd.isna(value)):
            continue
        norm = normalise_text(value)
        if not norm:
            continue
        mapping.setdefault(norm, []).append(str(ligand_id))
    return mapping


def build_name_map(ligands: pd.DataFrame) -> Dict[str, List[str]]:
    name_map: Dict[str, List[str]] = {}
    for _, row in ligands.iterrows():
        ligand_id = row.get("Ligand ID")
        if pd.isna(ligand_id):
            continue
        for field in ("Name", "INN"):
            value = row.get(field)
            norm = normalise_text(value)
            if norm:
                name_map.setdefault(norm, []).append(str(ligand_id))
        synonyms = row.get("Synonyms")
        if isinstance(synonyms, str):
            for part in synonyms.split("|"):
                norm = normalise_text(part)
                if norm:
                    name_map.setdefault(norm, []).append(str(ligand_id))
    return name_map


def summarise_matches(match_ids: pd.Series) -> Dict[str, int]:
    def has_match(value: object) -> bool:
        return isinstance(value, list) and len(value) > 0

    def is_ambiguous(value: object) -> bool:
        return isinstance(value, list) and len(value) > 1

    return {
        "matched": int(match_ids.apply(has_match).sum()),
        "ambiguous": int(match_ids.apply(is_ambiguous).sum()),
    }


def join_unique(series: pd.Series) -> str:
    values = [str(v).strip() for v in series.dropna().unique().tolist() if str(v).strip()]
    if not values:
        return ""
    return "|".join(sorted(set(values)))


def choose_ligand_match(row: pd.Series) -> pd.Series:
    for field, method in (
        ("match_smiles_ids", "SMILES"),
        ("match_cas_ids", "CAS"),
        ("match_name_ids", "NAME"),
    ):
        value = row.get(field)
        if isinstance(value, list) and value:
            if len(value) == 1:
                return pd.Series(
                    {
                        "ligand_match_method": method,
                        "ligand_match_status": "matched",
                        "gtp_ligand_id": value[0],
                        "ligand_match_ids": "|".join(value),
                    }
                )
            return pd.Series(
                {
                    "ligand_match_method": method,
                    "ligand_match_status": "ambiguous",
                    "gtp_ligand_id": None,
                    "ligand_match_ids": "|".join(value),
                }
            )
    return pd.Series(
        {
            "ligand_match_method": None,
            "ligand_match_status": "unmatched",
            "gtp_ligand_id": None,
            "ligand_match_ids": "",
        }
    )


def main() -> None:
    ki_path = DATA_DIR / "KiDatabase.csv"
    interactions_path = DATA_DIR / "interactions.csv"
    targets_path = DATA_DIR / "targets_and_families.csv"
    ligands_path = DATA_DIR / "ligands.csv"
    ligand_map_path = DATA_DIR / "ligand_id_mapping.csv"

    for path in [ki_path, interactions_path, targets_path, ligands_path, ligand_map_path]:
        if not path.exists():
            logger.error("Missing data file: %s", path)
            raise FileNotFoundError(path)

    ki = read_pdsp_csv(ki_path)
    interactions = read_gtp_csv(interactions_path)
    targets = read_gtp_csv(targets_path)
    ligands = read_gtp_csv(ligands_path)
    ligand_map = read_gtp_csv(ligand_map_path)

    require_columns(
        ki,
        ["Unigene", "Ligand Name", "SMILES", "CAS", "ki Val"],
        "PDSP KiDatabase",
    )
    require_columns(
        interactions,
        ["Target UniProt ID", "Ligand ID", "Action", "Type", "PubMed ID"],
        "GtoPdb interactions",
    )
    require_columns(
        targets,
        ["Target id", "HGNC symbol", "Human SwissProt", "Target name"],
        "GtoPdb targets",
    )
    require_columns(
        ligands,
        ["Ligand ID", "Name", "SMILES", "InChIKey", "Synonyms"],
        "GtoPdb ligands",
    )
    require_columns(
        ligand_map,
        ["Ligand id", "CAS", "DrugBank ID", "Drug Central ID"],
        "GtoPdb ligand_id_mapping",
    )

    # Normalise key column names
    ligand_map = ligand_map.rename(columns={"Ligand id": "Ligand ID"})

    # Target join (PDSP -> GtoPdb targets)
    target_map = targets[["HGNC symbol", "Human SwissProt", "Target id", "Target name"]].copy()
    target_map = target_map.dropna(subset=["HGNC symbol", "Human SwissProt"])
    target_map["HGNC symbol_norm"] = target_map["HGNC symbol"].map(normalise_symbol)

    ki_targets = ki.copy()
    ki_targets["Unigene_norm"] = ki_targets["Unigene"].map(normalise_symbol)
    ki_targets = ki_targets.merge(
        target_map,
        left_on="Unigene_norm",
        right_on="HGNC symbol_norm",
        how="left",
        suffixes=("", "_gtp"),
    )
    ki_targets = ki_targets.rename(
        columns={
            "Human SwissProt": "target_uniprot",
            "Target id": "gtp_target_id",
            "Target name": "gtp_target_name",
        }
    )
    ki_targets["target_uniprot"] = ki_targets["target_uniprot"].astype("string")

    # Ligand join (PDSP -> GtoPdb ligands)
    ligands_combined = ligands.merge(ligand_map, on="Ligand ID", how="left", suffixes=("", "_map"))

    gtp_smiles_map = build_lookup_map(ligands_combined["SMILES"], ligands_combined["Ligand ID"])
    gtp_cas_map = build_lookup_map(ligands_combined["CAS"], ligands_combined["Ligand ID"])
    gtp_name_map = build_name_map(ligands_combined)

    ki_targets["SMILES_norm"] = ki_targets["SMILES"].map(normalise_text)
    ki_targets["CAS_norm"] = ki_targets["CAS"].map(normalise_text)
    ki_targets["Ligand_name_norm"] = ki_targets["Ligand Name"].map(normalise_text)

    ki_targets["match_smiles_ids"] = ki_targets["SMILES_norm"].map(gtp_smiles_map)
    ki_targets["match_cas_ids"] = ki_targets["CAS_norm"].map(gtp_cas_map)
    ki_targets["match_name_ids"] = ki_targets["Ligand_name_norm"].map(gtp_name_map)

    match_info = ki_targets.apply(choose_ligand_match, axis=1)
    ki_targets = pd.concat([ki_targets, match_info], axis=1)

    # Attach GtoPdb ligand metadata
    ligands_meta = ligands[["Ligand ID", "Name", "SMILES", "InChIKey", "Type", "Approved", "Withdrawn"]].copy()
    ligands_meta = ligands_meta.rename(
        columns={
            "Ligand ID": "gtp_ligand_id",
            "Name": "gtp_ligand_name",
            "SMILES": "gtp_ligand_smiles",
            "InChIKey": "gtp_ligand_inchikey",
            "Type": "gtp_ligand_type",
            "Approved": "gtp_ligand_approved",
            "Withdrawn": "gtp_ligand_withdrawn",
        }
    )

    ki_targets["gtp_ligand_id"] = ki_targets["gtp_ligand_id"].astype("string")
    ligands_meta["gtp_ligand_id"] = ligands_meta["gtp_ligand_id"].astype("string")
    ki_targets = ki_targets.merge(ligands_meta, on="gtp_ligand_id", how="left")

    # Aggregate functional interactions
    interactions_use = interactions.dropna(subset=["Target UniProt ID", "Ligand ID"]).copy()
    interactions_use["Ligand ID"] = interactions_use["Ligand ID"].astype("string")
    interactions_use["Target UniProt ID"] = interactions_use["Target UniProt ID"].astype("string")

    interactions_agg = (
        interactions_use.groupby(["Ligand ID", "Target UniProt ID"], as_index=False)
        .agg(
            gtp_interaction_count=("Ligand ID", "size"),
            gtp_action_set=("Action", join_unique),
            gtp_type_set=("Type", join_unique),
            gtp_selectivity_set=("Selectivity", join_unique),
            gtp_primary_target_set=("Primary Target", join_unique),
            gtp_affinity_units_set=("Affinity Units", join_unique),
            gtp_affinity_median_set=("Affinity Median", join_unique),
            gtp_pubmed_ids=("PubMed ID", join_unique),
        )
    )

    joined = ki_targets.merge(
        interactions_agg,
        how="left",
        left_on=["gtp_ligand_id", "target_uniprot"],
        right_on=["Ligand ID", "Target UniProt ID"],
    )
    joined = joined.drop(columns=["Ligand ID", "Target UniProt ID"], errors="ignore")

    # Diagnostics
    total_rows = len(joined)
    target_matched = joined["target_uniprot"].notna().sum()
    target_unmatched = total_rows - target_matched

    ligand_match_counts = joined["ligand_match_status"].value_counts(dropna=False).to_dict()
    smiles_summary = summarise_matches(joined["match_smiles_ids"])
    cas_summary = summarise_matches(joined["match_cas_ids"])
    name_summary = summarise_matches(joined["match_name_ids"])

    functional_joined = joined["gtp_action_set"].astype(str).str.len().gt(0).sum()
    eligible_for_functional = joined[["gtp_ligand_id", "target_uniprot"]].notna().all(axis=1).sum()

    ambiguous_ligands = joined[joined["ligand_match_status"] == "ambiguous"]
    unmatched_ligands = joined[joined["ligand_match_status"] == "unmatched"]

    ambiguous_sample = (
        ambiguous_ligands[["Ligand Name", "SMILES", "CAS", "ligand_match_method", "ligand_match_ids"]]
        .head(20)
        .to_dict(orient="records")
    )
    unmatched_sample = (
        unmatched_ligands[["Ligand Name", "SMILES", "CAS"]].head(20).to_dict(orient="records")
    )

    if not ambiguous_ligands.empty:
        logger.warning("Ambiguous ligand matches detected: %s", len(ambiguous_ligands))

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_txt = REPORT_DIR / f"join_metrics_{timestamp}.txt"
    output_csv = REPORT_DIR / f"joined_pdsp_gtopdb_{timestamp}.csv"

    with report_txt.open("w", encoding="utf-8") as handle:
        handle.write("Join Metrics Report\n")
        handle.write(f"Generated: {datetime.now().isoformat()}\n\n")

        handle.write("Targets\n")
        handle.write(f"Total PDSP records: {total_rows}\n")
        handle.write(f"Matched to UniProt: {target_matched}\n")
        handle.write(f"Unmatched: {target_unmatched}\n\n")

        handle.write("Ligands\n")
        handle.write(f"Matched: {ligand_match_counts.get('matched', 0)}\n")
        handle.write(f"Ambiguous: {ligand_match_counts.get('ambiguous', 0)}\n")
        handle.write(f"Unmatched: {ligand_match_counts.get('unmatched', 0)}\n")
        handle.write("\nMatch method coverage (includes ambiguous matches):\n")
        handle.write(
            f"SMILES matched: {smiles_summary['matched']} (ambiguous: {smiles_summary['ambiguous']})\n"
        )
        handle.write(
            f"CAS matched: {cas_summary['matched']} (ambiguous: {cas_summary['ambiguous']})\n"
        )
        handle.write(
            f"Name matched: {name_summary['matched']} (ambiguous: {name_summary['ambiguous']})\n\n"
        )

        handle.write("Functional annotations\n")
        handle.write(f"Eligible for functional join: {eligible_for_functional}\n")
        handle.write(f"Rows with functional data: {functional_joined}\n\n")

        handle.write("Sample ambiguous ligand matches (first 20)\n")
        if not ambiguous_sample:
            handle.write("(none)\n\n")
        else:
            for row in ambiguous_sample:
                handle.write(f"- {row}\n")
            handle.write("\n")

        handle.write("Sample unmatched ligand records (first 20)\n")
        if not unmatched_sample:
            handle.write("(none)\n\n")
        else:
            for row in unmatched_sample:
                handle.write(f"- {row}\n")
            handle.write("\n")

    joined.to_csv(output_csv, index=False)

    logger.info("Wrote metrics report: %s", report_txt)
    logger.info("Wrote joined CSV: %s", output_csv)

    print("Join complete")
    print(f"Targets matched to UniProt: {target_matched}/{total_rows}")
    print(f"Ligands matched: {ligand_match_counts.get('matched', 0)}/{total_rows}")
    print(f"Rows with functional data: {functional_joined}/{total_rows}")
    print(f"Report (txt): {report_txt}")
    print(f"Joined CSV: {output_csv}")


if __name__ == "__main__":
    main()
