#!/usr/bin/env python
"""Join PDSP Ki data with GtoPdb targets/ligands/interactions and emit datasets plus diagnostics."""
from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
import numpy as np
try:
    from rdkit import Chem as _Chem
    from rdkit.Chem import inchi as _inchi
    from rdkit.Chem.SaltRemover import SaltRemover as _SaltRemover
except ImportError:  # pragma: no cover - optional dependency
    _Chem = None
    _inchi = None
    _SaltRemover = None


DATA_DIR = Path("data")
REPORT_DIR = Path("reports")
METRICS_DIR = REPORT_DIR / "metrics"
OUTPUT_DIR = DATA_DIR / "joined"
REPORT_DIR.mkdir(exist_ok=True)
METRICS_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)

logger = logging.getLogger("join_data")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

Chem = _Chem
inchi = _inchi
SaltRemover = _SaltRemover

HAS_RDKIT = Chem is not None and inchi is not None and SaltRemover is not None
if HAS_RDKIT and SaltRemover is not None:
    SALT_REMOVER = SaltRemover()
else:
    SALT_REMOVER = None


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


def canonicalise_smiles(value: object) -> Tuple[Optional[str], Optional[str], bool]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None, None, False
    if not HAS_RDKIT or Chem is None or inchi is None:
        return None, None, False
    assert Chem is not None and inchi is not None
    smiles = str(value).strip()
    if not smiles:
        return None, None, False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, False
    stripped = False
    if SALT_REMOVER is not None:
        stripped_mol = SALT_REMOVER.StripMol(mol, dontRemoveEverything=True)
        if stripped_mol is not None:
            stripped = stripped_mol.GetNumAtoms() < mol.GetNumAtoms()
            mol = stripped_mol
    canon = Chem.MolToSmiles(mol, isomericSmiles=True)
    key = inchi.MolToInchiKey(mol)
    return canon, key, stripped


def build_lookup_map(values: pd.Series, ligand_ids: pd.Series) -> Dict[str, List[str]]:
    mapping: Dict[str, set[str]] = {}
    for value, ligand_id in zip(values, ligand_ids):
        if value is None or (isinstance(value, float) and pd.isna(value)):
            continue
        norm = normalise_text(value)
        if not norm:
            continue
        mapping.setdefault(norm, set()).add(str(ligand_id))
    return {key: sorted(values) for key, values in mapping.items()}


def build_name_map(ligands: pd.DataFrame) -> Dict[str, List[str]]:
    name_map: Dict[str, set[str]] = {}
    for _, row in ligands.iterrows():
        ligand_id = row.get("Ligand ID")
        if pd.isna(ligand_id):
            continue
        for field in ("Name", "INN"):
            value = row.get(field)
            norm = normalise_text(value)
            if norm:
                name_map.setdefault(norm, set()).add(str(ligand_id))
        synonyms = row.get("Synonyms")
        if isinstance(synonyms, str):
            for part in synonyms.split("|"):
                norm = normalise_text(part)
                if norm:
                    name_map.setdefault(norm, set()).add(str(ligand_id))
    return {key: sorted(values) for key, values in name_map.items()}


def summarise_matches(match_ids: pd.Series) -> Dict[str, int]:
    def normalise_ids(value: object) -> List[str]:
        if not isinstance(value, list):
            return []
        return sorted({str(item) for item in value if str(item)})

    def has_match(value: object) -> bool:
        return len(normalise_ids(value)) > 0

    def is_ambiguous(value: object) -> bool:
        return len(normalise_ids(value)) > 1

    return {
        "matched": int(match_ids.apply(has_match).sum()),
        "ambiguous": int(match_ids.apply(is_ambiguous).sum()),
    }


def join_unique(series: pd.Series) -> str:
    values = [str(v).strip() for v in series.dropna().unique().tolist() if str(v).strip()]
    if not values:
        return ""
    return "|".join(sorted(set(values)))


def dedupe_ids(values: object) -> List[str]:
    if not isinstance(values, list):
        return []
    return sorted({str(item) for item in values if str(item)})


def select_cols(df: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    available = [column for column in columns if column in df.columns]
    return df[available].copy()


def write_csv_overwrite(df: pd.DataFrame, path: Path) -> None:
    if path.exists():
        logger.warning("Overwriting existing file: %s", path)
    df.to_csv(path, index=False)


def choose_ligand_match(row: pd.Series) -> pd.Series:
    for field, method in (
        ("match_inchikey_ids", "INCHIKEY"),
        ("match_canon_smiles_ids", "CANON_SMILES"),
        ("match_smiles_ids", "SMILES"),
        ("match_cas_ids", "CAS"),
        ("match_name_ids", "NAME"),
    ):
        value = row.get(field)
        ids = dedupe_ids(value)
        if ids:
            if len(ids) == 1:
                return pd.Series(
                    {
                        "ligand_match_method": method,
                        "ligand_match_status": "matched",
                        "gtp_ligand_id": ids[0],
                        "ligand_match_ids": "|".join(ids),
                    }
                )
            return pd.Series(
                {
                    "ligand_match_method": method,
                    "ligand_match_status": "ambiguous",
                    "gtp_ligand_id": None,
                    "ligand_match_ids": "|".join(ids),
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

    if not HAS_RDKIT:
        logger.warning("RDKit not available. Canonical SMILES and InChIKey matching will be skipped.")
    else:
        assert Chem is not None and inchi is not None
        test_mol = Chem.MolFromSmiles("CCO")
        test_key = inchi.MolToInchiKey(test_mol) if test_mol is not None else None
        if not test_key:
            logger.warning(
                "RDKit InChIKey generation appears unavailable. InChIKey matching will rely on raw GtoPdb values."
            )

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

    if HAS_RDKIT:
        ligands_combined[["gtp_canon_smiles", "gtp_inchikey_canon", "gtp_salt_stripped"]] = ligands_combined[
            "SMILES"
        ].apply(lambda value: pd.Series(canonicalise_smiles(value)))
    else:
        ligands_combined["gtp_canon_smiles"] = None
        ligands_combined["gtp_inchikey_canon"] = None
        ligands_combined["gtp_salt_stripped"] = False

    ligands_combined["gtp_inchikey_raw"] = ligands_combined["InChIKey"].map(normalise_symbol)
    ligands_combined["gtp_inchikey_for_match"] = ligands_combined["gtp_inchikey_canon"].fillna(
        ligands_combined["gtp_inchikey_raw"]
    )

    gtp_inchikey_map = build_lookup_map(ligands_combined["gtp_inchikey_for_match"], ligands_combined["Ligand ID"])
    gtp_canon_smiles_map = build_lookup_map(ligands_combined["gtp_canon_smiles"], ligands_combined["Ligand ID"])
    gtp_smiles_map = build_lookup_map(ligands_combined["SMILES"], ligands_combined["Ligand ID"])
    gtp_cas_map = build_lookup_map(ligands_combined["CAS"], ligands_combined["Ligand ID"])
    gtp_name_map = build_name_map(ligands_combined)

    ki_targets["SMILES_norm"] = ki_targets["SMILES"].map(normalise_text)
    ki_targets["CAS_norm"] = ki_targets["CAS"].map(normalise_text)
    ki_targets["Ligand_name_norm"] = ki_targets["Ligand Name"].map(normalise_text)
    if HAS_RDKIT:
        ki_targets[["pdsp_canon_smiles", "pdsp_inchikey_canon", "pdsp_salt_stripped"]] = ki_targets[
            "SMILES"
        ].apply(lambda value: pd.Series(canonicalise_smiles(value)))
    else:
        ki_targets["pdsp_canon_smiles"] = None
        ki_targets["pdsp_inchikey_canon"] = None
        ki_targets["pdsp_salt_stripped"] = False

    ki_targets["pdsp_inchikey_norm"] = ki_targets["pdsp_inchikey_canon"].map(normalise_text)
    ki_targets["pdsp_canon_smiles_norm"] = ki_targets["pdsp_canon_smiles"].map(normalise_text)
    ki_targets["match_inchikey_ids"] = ki_targets["pdsp_inchikey_norm"].map(gtp_inchikey_map)
    ki_targets["match_canon_smiles_ids"] = ki_targets["pdsp_canon_smiles_norm"].map(gtp_canon_smiles_map)
    ki_targets["match_smiles_ids"] = ki_targets["SMILES_norm"].map(gtp_smiles_map)
    ki_targets["match_cas_ids"] = ki_targets["CAS_norm"].map(gtp_cas_map)
    ki_targets["match_name_ids"] = ki_targets["Ligand_name_norm"].map(gtp_name_map)

    match_info = ki_targets.apply(choose_ligand_match, axis=1)
    ki_targets = pd.concat([ki_targets, match_info], axis=1)

    # Attach GtoPdb ligand metadata
    ligands_meta = ligands_combined[
        [
            "Ligand ID",
            "Name",
            "SMILES",
            "InChIKey",
            "gtp_canon_smiles",
            "gtp_inchikey_raw",
            "gtp_inchikey_canon",
            "gtp_salt_stripped",
            "Type",
            "Approved",
            "Withdrawn",
        ]
    ].copy()
    ligands_meta = ligands_meta.rename(
        columns={
            "Ligand ID": "gtp_ligand_id",
            "Name": "gtp_ligand_name",
            "SMILES": "gtp_ligand_smiles",
            "InChIKey": "gtp_ligand_inchikey_raw",
            "gtp_canon_smiles": "gtp_ligand_canon_smiles",
            "gtp_inchikey_raw": "gtp_ligand_inchikey_norm",
            "gtp_inchikey_canon": "gtp_ligand_inchikey_canon",
            "gtp_salt_stripped": "gtp_ligand_salt_stripped",
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
    if "Ligand ID_x" in joined.columns:
        joined = joined.rename(columns={"Ligand ID_x": "pdsp_ligand_id"})
    if "Ligand ID_y" in joined.columns:
        joined = joined.drop(columns=["Ligand ID_y"], errors="ignore")

    joined["ki_val_numeric"] = pd.to_numeric(joined["ki Val"], errors="coerce")
    pki = pd.Series(np.nan, index=joined.index, dtype="float64")
    valid_ki = joined["ki_val_numeric"] > 0
    if valid_ki.any():
        pki.loc[valid_ki] = 9 - np.log10(joined.loc[valid_ki, "ki_val_numeric"].astype(float))
    joined["pKi_assuming_nM"] = pki

    binding_base = joined[joined["gtp_ligand_id"].notna() & joined["target_uniprot"].notna()].copy()
    binding_summary = (
        binding_base.groupby(["gtp_ligand_id", "target_uniprot"], as_index=False)
        .agg(
            gtp_ligand_name=("gtp_ligand_name", join_unique),
            gtp_target_name=("gtp_target_name", join_unique),
            pdsp_ligand_names=("Ligand Name", join_unique),
            pdsp_unigene_symbols=("Unigene", join_unique),
            ligand_match_method=("ligand_match_method", join_unique),
            ki_count=("ki_val_numeric", "count"),
            ki_median=("ki_val_numeric", "median"),
            ki_min=("ki_val_numeric", "min"),
            ki_max=("ki_val_numeric", "max"),
            pKi_median=("pKi_assuming_nM", "median"),
            pKi_min=("pKi_assuming_nM", "min"),
            pKi_max=("pKi_assuming_nM", "max"),
        )
    )

    gtp_smiles_present = ligands_combined["SMILES"].notna() & (
        ligands_combined["SMILES"].astype(str).str.strip() != ""
    )
    pdsp_smiles_present = ki_targets["SMILES"].notna() & (
        ki_targets["SMILES"].astype(str).str.strip() != ""
    )
    if HAS_RDKIT:
        gtp_smiles_parse_failures = (gtp_smiles_present & ligands_combined["gtp_canon_smiles"].isna()).sum()
        pdsp_smiles_parse_failures = (pdsp_smiles_present & ki_targets["pdsp_canon_smiles"].isna()).sum()
        gtp_inchikey_raw_count = ligands_combined["gtp_inchikey_raw"].notna().sum()
        gtp_inchikey_canon_count = ligands_combined["gtp_inchikey_canon"].notna().sum()
        pdsp_inchikey_canon_count = ki_targets["pdsp_inchikey_canon"].notna().sum()
    else:
        gtp_smiles_parse_failures = None
        pdsp_smiles_parse_failures = None
        gtp_inchikey_raw_count = None
        gtp_inchikey_canon_count = None
        pdsp_inchikey_canon_count = None

    targets_out = select_cols(
        targets,
        [
            "Type",
            "Family id",
            "Family name",
            "Target id",
            "Target name",
            "HGNC symbol",
            "Human SwissProt",
        ],
    ).rename(
        columns={
            "Target id": "gtp_target_id",
            "Target name": "gtp_target_name",
            "HGNC symbol": "hgnc_symbol",
            "Human SwissProt": "target_uniprot",
            "Family id": "gtp_family_id",
            "Family name": "gtp_family_name",
        }
    )
    if "target_uniprot" in targets_out.columns:
        targets_out = targets_out.dropna(subset=["target_uniprot"])

    ligands_out = select_cols(
        ligands_combined,
        [
            "Ligand ID",
            "Name",
            "SMILES",
            "InChIKey",
            "gtp_canon_smiles",
            "gtp_inchikey_raw",
            "gtp_inchikey_canon",
            "gtp_salt_stripped",
            "Type",
            "Approved",
            "Withdrawn",
            "INN",
            "Synonyms",
            "CAS",
            "DrugBank ID",
            "Drug Central ID",
        ],
    ).rename(
        columns={
            "Ligand ID": "gtp_ligand_id",
            "Name": "gtp_ligand_name",
            "SMILES": "gtp_ligand_smiles",
            "InChIKey": "gtp_ligand_inchikey_raw",
            "gtp_canon_smiles": "gtp_ligand_canon_smiles",
            "gtp_inchikey_raw": "gtp_ligand_inchikey_norm",
            "gtp_inchikey_canon": "gtp_ligand_inchikey_canon",
            "gtp_salt_stripped": "gtp_ligand_salt_stripped",
            "Type": "gtp_ligand_type",
            "Approved": "gtp_ligand_approved",
            "Withdrawn": "gtp_ligand_withdrawn",
        }
    )
    if "gtp_ligand_id" in ligands_out.columns:
        ligands_out = ligands_out.dropna(subset=["gtp_ligand_id"])

    functional_out = interactions_agg.rename(
        columns={"Ligand ID": "gtp_ligand_id", "Target UniProt ID": "target_uniprot"}
    )
    functional_out["gtp_ligand_id"] = functional_out["gtp_ligand_id"].astype("string")
    functional_out["target_uniprot"] = functional_out["target_uniprot"].astype("string")

    # Diagnostics
    total_rows = len(joined)
    target_matched = joined["target_uniprot"].notna().sum()
    target_unmatched = total_rows - target_matched

    ligand_match_counts = joined["ligand_match_status"].value_counts(dropna=False).to_dict()
    inchikey_summary = summarise_matches(joined["match_inchikey_ids"])
    canon_smiles_summary = summarise_matches(joined["match_canon_smiles_ids"])
    smiles_summary = summarise_matches(joined["match_smiles_ids"])
    cas_summary = summarise_matches(joined["match_cas_ids"])
    name_summary = summarise_matches(joined["match_name_ids"])

    functional_joined = joined["gtp_interaction_count"].fillna(0).astype(int).gt(0).sum()
    eligible_for_functional = joined[["gtp_ligand_id", "target_uniprot"]].notna().all(axis=1).sum()

    ambiguous_ligands = joined[joined["ligand_match_status"] == "ambiguous"]
    unmatched_ligands = joined[joined["ligand_match_status"] == "unmatched"]

    ambiguous_sample = (
        ambiguous_ligands[["Ligand Name", "SMILES", "CAS", "ligand_match_method", "ligand_match_ids"]]
        .drop_duplicates()
        .head(20)
        .to_dict(orient="records")
    )
    unmatched_sample = (
        unmatched_ligands[["Ligand Name", "SMILES", "CAS"]]
        .drop_duplicates()
        .head(20)
        .to_dict(orient="records")
    )

    if not ambiguous_ligands.empty:
        logger.warning("Ambiguous ligand matches detected: %s", len(ambiguous_ligands))

    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    report_txt = METRICS_DIR / f"join_metrics_{timestamp}.txt"
    output_csv = OUTPUT_DIR / "joined_pdsp_gtopdb.csv"
    targets_csv = OUTPUT_DIR / "targets_normalised.csv"
    ligands_csv = OUTPUT_DIR / "ligands_normalised.csv"
    binding_csv = OUTPUT_DIR / "binding_measurements.csv"
    functional_csv = OUTPUT_DIR / "functional_interactions.csv"

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
            f"InChIKey matched: {inchikey_summary['matched']} (ambiguous: {inchikey_summary['ambiguous']})\n"
        )
        handle.write(
            "Canonical SMILES matched: "
            f"{canon_smiles_summary['matched']} (ambiguous: {canon_smiles_summary['ambiguous']})\n"
        )
        handle.write(
            f"SMILES matched: {smiles_summary['matched']} (ambiguous: {smiles_summary['ambiguous']})\n"
        )
        handle.write(
            f"CAS matched: {cas_summary['matched']} (ambiguous: {cas_summary['ambiguous']})\n"
        )
        handle.write(
            f"Name matched: {name_summary['matched']} (ambiguous: {name_summary['ambiguous']})\n\n"
        )

        handle.write("Structure canonicalisation (RDKit)\n")
        handle.write(f"GtoPdb SMILES present: {int(gtp_smiles_present.sum())}\n")
        handle.write(f"PDSP SMILES present: {int(pdsp_smiles_present.sum())}\n")
        if HAS_RDKIT:
            handle.write(f"GtoPdb SMILES parse failures: {gtp_smiles_parse_failures}\n")
            handle.write(f"PDSP SMILES parse failures: {pdsp_smiles_parse_failures}\n")
            handle.write(f"GtoPdb InChIKey raw present: {gtp_inchikey_raw_count}\n")
            handle.write(f"GtoPdb InChIKey canonical: {gtp_inchikey_canon_count}\n")
            handle.write(f"PDSP InChIKey canonical: {pdsp_inchikey_canon_count}\n\n")
        else:
            handle.write("RDKit not available; canonicalisation skipped.\n\n")

        handle.write("Functional annotations\n")
        handle.write(f"Eligible for functional join: {eligible_for_functional}\n")
        handle.write(f"Rows with functional data: {functional_joined}\n\n")
        handle.write("Binding summaries\n")
        handle.write("pKi computed as 9 - log10(Ki) assuming Ki values are in nM.\n\n")
        handle.write(f"Binding summary rows: {len(binding_summary)}\n\n")

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

    write_csv_overwrite(joined, output_csv)
    write_csv_overwrite(targets_out.drop_duplicates(), targets_csv)
    write_csv_overwrite(ligands_out.drop_duplicates(), ligands_csv)
    write_csv_overwrite(binding_summary, binding_csv)
    write_csv_overwrite(functional_out.drop_duplicates(), functional_csv)

    logger.info("Wrote metrics report: %s", report_txt)
    logger.info("Wrote joined CSV: %s", output_csv)
    logger.info("Wrote targets CSV: %s", targets_csv)
    logger.info("Wrote ligands CSV: %s", ligands_csv)
    logger.info("Wrote binding CSV: %s", binding_csv)
    logger.info("Wrote functional CSV: %s", functional_csv)

    print("Join complete")
    print(f"Targets matched to UniProt: {target_matched}/{total_rows}")
    print(f"Ligands matched: {ligand_match_counts.get('matched', 0)}/{total_rows}")
    print(f"Rows with functional data: {functional_joined}/{total_rows}")
    print(f"Report (txt): {report_txt}")
    print(f"Joined CSV: {output_csv}")
    print(f"Targets CSV: {targets_csv}")
    print(f"Ligands CSV: {ligands_csv}")
    print(f"Binding CSV: {binding_csv}")
    print(f"Functional CSV: {functional_csv}")


if __name__ == "__main__":
    main()
