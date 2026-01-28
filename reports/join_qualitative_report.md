# PDSP + GtoPdb Join Notes

This note documents the datasets in this repository and the join logic used by the scripted pipeline. It is written to support the brief for a receptor‑informed medication switching tool, with explicit separation between binding evidence (PDSP) and functional direction (GtoPdb).

## Data sources in this repo

- **PDSP Ki Database** (`data/KiDatabase.csv`)
  - Long‑format binding measurements (ligand × target × Ki).
  - Target identifiers are provided as gene symbols in the `Unigene` column.
  - Ligand identifiers are variable; SMILES and CAS are often missing.

- **GtoPdb interactions** (`data/interactions.csv`)
  - Curated ligand–target interactions with functional actions (e.g. agonist/antagonist), affinity summaries, and PubMed IDs.
  - Includes `Target UniProt ID` and `Ligand ID` columns that act as stable join keys.
  - The first line is a metadata row and must be skipped when parsing.

- **GtoPdb targets & families** (`data/targets_and_families.csv`)
  - Target dictionary, including HGNC symbols and human UniProt accessions (`Human SwissProt`).
  - The first line is a metadata row and must be skipped when parsing.

- **GtoPdb ligands** (`data/ligands.csv`) and **Ligand ID mapping** (`data/ligand_id_mapping.csv`)
  - Ligand dictionary, names, structures (SMILES, InChIKey), and identifier cross‑references (CAS, DrugBank, etc.).
  - Each file begins with a metadata row and requires `skiprows=1`.

## Join strategy (implemented in `scripts/join_data.py`)

### Targets (PDSP → UniProt)
1. PDSP target symbols (`Unigene`) are normalised to uppercase.
2. They are matched to GtoPdb `HGNC symbol`.
3. The corresponding `Human SwissProt` entry is used as the UniProt‑normalised target ID.

This yields a robust path to UniProt for most PDSP targets and enables a clean join to GtoPdb interactions via `Target UniProt ID`.

### Ligands (PDSP → GtoPdb)
Because PDSP ligand identifiers are sparse, matches are attempted in descending priority:

1. **Exact SMILES match** (PDSP `SMILES` → GtoPdb `SMILES`)
2. **CAS match** (PDSP `CAS` → GtoPdb `CAS` via the ligand mapping file)
3. **Name/synonym match** (PDSP `Ligand Name` → GtoPdb `Name`, `INN`, and `Synonyms`)

If multiple GtoPdb ligands match a single PDSP record, the match is labelled **ambiguous** and no ligand ID is assigned. The pipeline de‑duplicates repeated IDs before deciding whether a match is truly ambiguous. Ambiguities are reported in the diagnostics output for manual review.

### Functional annotations (GtoPdb interactions)
After target and ligand normalisation, PDSP records are left‑joined to an aggregated GtoPdb interactions table keyed on:

- `gtp_ligand_id`
- `target_uniprot`

Interaction records are aggregated to preserve provenance while avoiding row explosions:

- Unique sets of `Action`, `Type`, `Selectivity`, and `PubMed ID` are collapsed into pipe‑delimited strings.
- A `gtp_interaction_count` is included to show how many interaction records were collapsed.

### Binding summaries (PDSP)
The join pipeline additionally produces an aggregated binding table, grouped by `(gtp_ligand_id, target_uniprot)`, with:

- `ki_count`, `ki_median`, `ki_min`, `ki_max`
- `pKi_median`, `pKi_min`, `pKi_max`

`pKi` is computed as `9 - log10(Ki)` under the assumption that Ki values are in nM. This assumption is documented in the diagnostics report for transparency.

## Outputs

- **Joined dataset**: a wide CSV produced by the script with PDSP binding rows enriched by:
  - UniProt‑normalised targets
  - GtoPdb ligand metadata
  - Aggregated functional annotations (where a match exists)

- **Normalised targets**: a CSV of GtoPdb targets with UniProt IDs.
- **Normalised ligands**: a CSV of GtoPdb ligands with structures and cross‑references.
- **Binding summaries**: a CSV aggregated by `(gtp_ligand_id, target_uniprot)` with Ki and pKi summaries.
- **Functional interactions**: a CSV of aggregated GtoPdb interaction evidence by `(gtp_ligand_id, target_uniprot)`.

- **Metrics report**: a timestamped `.txt` file containing counts for:
  - target matches
  - ligand matches (matched / ambiguous / unmatched)
  - eligible and successful functional joins
  - binding summary row count
  - example ambiguous and unmatched ligand records

## Caveats and practical implications

- Target mapping is generally robust; most misses will be symbol discrepancies or non‑human entries.
- Ligand mapping is the limiting factor because many PDSP entries lack structure or CAS identifiers.
- Name‑based matching is useful but risky and should be interpreted carefully.

These constraints are expected at this stage; the diagnostics report is intended to make the join quality transparent and actionable.
