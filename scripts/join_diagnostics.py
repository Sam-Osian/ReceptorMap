#!/usr/bin/env python
"""Generate join diagnostics for PDSP + GtoPdb datasets.

Outputs a text metrics report.
"""
from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd


DATA_DIR = Path("data")
REPORT_DIR = Path("reports")
REPORT_DIR.mkdir(exist_ok=True)

logger = logging.getLogger("join_diagnostics")
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
        mapping.setdefault(norm, []).append(ligand_id)
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
                name_map.setdefault(norm, []).append(ligand_id)
        synonyms = row.get("Synonyms")
        if isinstance(synonyms, str):
            for part in synonyms.split("|"):
                norm = normalise_text(part)
                if norm:
                    name_map.setdefault(norm, []).append(ligand_id)
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


def main() -> None:
    ki_path = DATA_DIR / "KiDatabase.csv"
    interactions_path = DATA_DIR / "interactions.csv"
    targets_path = DATA_DIR / "targets_and_families.csv"
    ligands_path = DATA_DIR / "ligands.csv"
    ligand_map_path = DATA_DIR / "ligand_id_mapping.csv"
    uniprot_map_path = DATA_DIR / "GtP_to_UniProt_mapping.csv"

    for path in [
        ki_path,
        interactions_path,
        targets_path,
        ligands_path,
        ligand_map_path,
        uniprot_map_path,
    ]:
        if not path.exists():
            logger.error("Missing data file: %s", path)
            raise FileNotFoundError(path)

    ki = read_pdsp_csv(ki_path)
    interactions = read_gtp_csv(interactions_path)
    targets = read_gtp_csv(targets_path)
    ligands = read_gtp_csv(ligands_path)
    ligand_map = read_gtp_csv(ligand_map_path)
    uniprot_map = read_gtp_csv(uniprot_map_path)

    require_columns(
        ki,
        ["Unigene", "Ligand Name", "SMILES", "CAS"],
        "PDSP KiDatabase",
    )
    require_columns(
        interactions,
        ["Target ID", "Target UniProt ID", "Ligand ID", "Ligand"],
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
    require_columns(
        uniprot_map,
        ["UniProtKB ID", "Species", "GtoPdb IUPHAR ID"],
        "GtoPdb UniProt mapping",
    )

    # Normalise key column names
    ligand_map = ligand_map.rename(columns={"Ligand id": "Ligand ID"})

    # Target join diagnostics (PDSP -> GtoPdb -> UniProt)
    symbol_map = targets[["HGNC symbol", "Human SwissProt", "Target id", "Target name"]].copy()
    symbol_map = symbol_map.dropna(subset=["HGNC symbol", "Human SwissProt"])
    symbol_map["HGNC symbol_norm"] = symbol_map["HGNC symbol"].map(normalise_symbol)

    ki_targets = ki.copy()
    ki_targets["Unigene_norm"] = ki_targets["Unigene"].map(normalise_symbol)
    ki_targets = ki_targets.merge(
        symbol_map,
        left_on="Unigene_norm",
        right_on="HGNC symbol_norm",
        how="left",
        suffixes=("", "_gtp"),
    )

    target_total = len(ki_targets)
    target_symbol_present = ki_targets["Unigene"].notna().sum()
    target_matched = ki_targets["Human SwissProt"].notna().sum()
    target_unmatched = ki_targets["Human SwissProt"].isna().sum()
    target_unmatched_symbols = (
        ki_targets.loc[ki_targets["Human SwissProt"].isna(), "Unigene_norm"]
        .dropna()
        .value_counts()
        .head(20)
    )

    # GtoPdb interactions target coverage
    interactions_target_missing = interactions["Target UniProt ID"].isna().sum()

    # Ligand join diagnostics (PDSP -> GtoPdb ligands)
    ligands_combined = ligands.merge(ligand_map, on="Ligand ID", how="left", suffixes=("", "_map"))

    ligands_combined["SMILES_norm"] = ligands_combined["SMILES"].map(normalise_text)
    ligands_combined["CAS_norm"] = ligands_combined["CAS"].map(normalise_text)

    gtp_smiles_map = build_lookup_map(ligands_combined["SMILES"], ligands_combined["Ligand ID"])
    gtp_cas_map = build_lookup_map(ligands_combined["CAS"], ligands_combined["Ligand ID"])
    gtp_name_map = build_name_map(ligands_combined)

    ki_ligands = ki.copy()
    ki_ligands["SMILES_norm"] = ki_ligands["SMILES"].map(normalise_text)
    ki_ligands["CAS_norm"] = ki_ligands["CAS"].map(normalise_text)
    ki_ligands["Ligand_name_norm"] = ki_ligands["Ligand Name"].map(normalise_text)

    ki_ligands["match_smiles_ids"] = ki_ligands["SMILES_norm"].map(gtp_smiles_map)
    ki_ligands["match_cas_ids"] = ki_ligands["CAS_norm"].map(gtp_cas_map)
    ki_ligands["match_name_ids"] = ki_ligands["Ligand_name_norm"].map(gtp_name_map)

    def choose_match(row: pd.Series) -> pd.Series:
        for field, method in (
            ("match_smiles_ids", "SMILES"),
            ("match_cas_ids", "CAS"),
            ("match_name_ids", "NAME"),
        ):
            value = row.get(field)
            if isinstance(value, list) and value:
                return pd.Series({"match_method": method, "match_ids": value})
        return pd.Series({"match_method": None, "match_ids": []})

    match_info = ki_ligands.apply(choose_match, axis=1)
    ki_ligands = pd.concat([ki_ligands, match_info], axis=1)

    ligand_total = len(ki_ligands)
    ligand_matched = (ki_ligands["match_method"].notna()).sum()
    ligand_unmatched = ligand_total - ligand_matched

    smiles_summary = summarise_matches(ki_ligands["match_smiles_ids"])
    cas_summary = summarise_matches(ki_ligands["match_cas_ids"])
    name_summary = summarise_matches(ki_ligands["match_name_ids"])

    ambiguous_ligands = ki_ligands[ki_ligands["match_ids"].apply(lambda x: isinstance(x, list) and len(x) > 1)]
    ambiguous_sample = (
        ambiguous_ligands[["Ligand Name", "SMILES", "CAS", "match_method", "match_ids"]]
        .head(20)
        .to_dict(orient="records")
    )

    unmatched_sample = (
        ki_ligands.loc[ki_ligands["match_method"].isna(), ["Ligand Name", "SMILES", "CAS"]]
        .head(20)
        .to_dict(orient="records")
    )

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_txt = REPORT_DIR / f"join_report_{timestamp}.txt"

    with report_txt.open("w", encoding="utf-8") as handle:
        handle.write("Join Diagnostics Report\n")
        handle.write(f"Generated: {datetime.now().isoformat()}\n\n")

        handle.write("Targets (PDSP -> GtoPdb -> UniProt)\n")
        handle.write(f"Total PDSP records: {target_total}\n")
        handle.write(f"PDSP records with Unigene: {target_symbol_present}\n")
        handle.write(f"Matched to Human SwissProt: {target_matched}\n")
        handle.write(f"Unmatched: {target_unmatched}\n")
        handle.write(f"GtoPdb interactions missing Target UniProt ID: {interactions_target_missing}\n\n")

        handle.write("Top unmatched target symbols (PDSP Unigene)\n")
        if target_unmatched_symbols.empty:
            handle.write("(none)\n\n")
        else:
            for symbol, count in target_unmatched_symbols.items():
                handle.write(f"- {symbol}: {count}\n")
            handle.write("\n")

        handle.write("Ligands (PDSP -> GtoPdb)\n")
        handle.write(f"Total PDSP records: {ligand_total}\n")
        handle.write(f"Matched by any method: {ligand_matched}\n")
        handle.write(f"Unmatched: {ligand_unmatched}\n")
        handle.write("\nMatch method coverage (includes ambiguous matches):\n")
        handle.write(f"SMILES matched: {smiles_summary['matched']} (ambiguous: {smiles_summary['ambiguous']})\n")
        handle.write(f"CAS matched: {cas_summary['matched']} (ambiguous: {cas_summary['ambiguous']})\n")
        handle.write(f"Name matched: {name_summary['matched']} (ambiguous: {name_summary['ambiguous']})\n\n")

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

    logger.info("Wrote report: %s", report_txt)

    print("Join diagnostics complete")
    print(f"Targets matched to UniProt: {target_matched}/{target_total}")
    print(f"Ligands matched to GtoPdb: {ligand_matched}/{ligand_total}")
    print(f"Report (txt): {report_txt}")


if __name__ == "__main__":
    main()
