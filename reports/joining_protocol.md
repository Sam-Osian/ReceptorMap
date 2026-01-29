# PDSP + GtoPdb Join Notes

This note explains the data sources, core domain concepts, and the join logic used by this repo's pipeline. The overall goal is to support rational within-class medication switching (starting with SSRIs) by making binding evidence and functional direction clear and comparable.

## Purpose and context

The tool aims to help clinical users compare drugs within the same class so they can avoid problematic targets while keeping desired similarity. The model separates:

- In vitro binding evidence (what receptors or targets a drug binds, and how strongly)
- Functional direction (whether the drug activates or blocks the target)

This repo uses PDSP and GtoPdb data to achieve this. Currently, it uses static downloads of these databases rather than a scraping pipeline.

## Core concepts (listing it here for convenience)

- Target: a protein that a drug can bind to (for example, SERT or a receptor). In datasets, targets are identified by gene symbols (HGNC) and UniProt accessions (the same target can have different names across different datasets).
- Ligand: a compound that binds a target. A drug is a ligand, but not all ligands are approved drugs.
- Binding affinity (Ki): a laboratory measurement of how tightly a ligand binds a target. Lower Ki means stronger binding.
- pKi: a transformed version of Ki used for comparability. It is computed as pKi = 9 - log10(Ki) if Ki is in nanomolar (nM). Higher pKi means stronger binding.
- Functional direction: what the ligand does to the target when it binds. Examples include agonist (activates), antagonist (blocks), inhibitor, or modulator.
- In vitro: experiments performed in laboratory systems, not in patients. These are informative but not always clinically relevant without exposure data.
- UniProt accession: a stable identifier for a protein, used to join across datasets.

## Data sources in this repo

- PDSP Ki Database (data/KiDatabase.csv)
  - Long format binding measurements (ligand x target x Ki).
  - Targets are given as gene symbols in the Unigene column.
  - Ligand identifiers are inconsistent; SMILES and CAS are often missing.

- GtoPdb interactions (data/interactions.csv)
  - Curated ligand target interactions with functional actions (agonist, antagonist, etc), affinity summaries, and PubMed IDs.
  - Includes Target UniProt ID and Ligand ID as join keys.
  - The first line is a metadata row and must be skipped when parsing.

- GtoPdb targets and families (data/targets_and_families.csv)
  - Target dictionary including HGNC symbols and human UniProt accessions (Human SwissProt).
  - The first line is a metadata row and must be skipped when parsing.

- GtoPdb ligands and ligand ID mapping (data/ligands.csv and data/ligand_id_mapping.csv)
  - Ligand dictionary, names, structures (SMILES, InChIKey), and cross references (CAS, DrugBank).
  - Each file begins with a metadata row and requires skiprows=1.

## Join strategy (implemented in scripts/join_data.py)

### Targets (PDSP to UniProt)
1. PDSP target symbols (Unigene) are normalised to uppercase.
2. They are matched to GtoPdb HGNC symbol.
3. The corresponding Human SwissProt entry is used as the UniProt normalised target ID.

This provides a route to UniProt for most PDSP targets and enables joins to GtoPdb interactions via Target UniProt ID.

### Ligands (PDSP to GtoPdb)
PDSP ligand identifiers are sparse, so matching uses a priority order. RDKit is used to canonicalise structures before matching, with two explicit decisions:

- Salt stripping is enabled to improve matches between free bases and salt forms.
- Stereochemistry is preserved (strict matching only, no stereo ignored fallback).

RDKit is a required dependency for canonical SMILES and InChIKey generation; if it is unavailable the pipeline falls back to string based matching only and logs a warning.

The priority order is:

1. InChIKey match (from canonicalised, salt stripped structures)
2. Canonical SMILES match (salt stripped, stereochemistry preserved)
3. Raw SMILES match (string equality)
4. CAS match (PDSP CAS to GtoPdb CAS via ligand mapping)
5. Name or synonym match (PDSP Ligand Name to GtoPdb Name, INN, Synonyms)

If multiple GtoPdb ligands match a single PDSP record, the match is labelled ambiguous and no ligand ID is assigned. Repeated IDs are de duplicated before deciding whether a match is truly ambiguous. Ambiguities are reported for manual review.

### Functional annotations (GtoPdb interactions)
After target and ligand normalisation, PDSP records are left joined to an aggregated GtoPdb interactions table keyed on:

- gtp_ligand_id
- target_uniprot

Interaction records are aggregated to preserve provenance while avoiding row explosions:

- Unique sets of Action, Type, Selectivity, and PubMed ID are collapsed into pipe delimited strings.
- gtp_interaction_count indicates how many interaction rows were collapsed.

### Binding summaries (PDSP)
The pipeline produces an aggregated binding table grouped by (gtp_ligand_id, target_uniprot), with:

- ki_count, ki_median, ki_min, ki_max
- pKi_median, pKi_min, pKi_max

pKi is computed as 9 - log10(Ki) under the assumption that Ki values are in nM. This assumption is documented in the diagnostics report.

## Outputs

- Joined dataset: a wide CSV with PDSP binding rows enriched by UniProt targets, GtoPdb ligand metadata, and aggregated functional annotations.
- Normalised targets: a CSV of GtoPdb targets with UniProt IDs.
- Normalised ligands: a CSV of GtoPdb ligands with structures and cross references.
- Binding summaries: a CSV aggregated by (gtp_ligand_id, target_uniprot) with Ki and pKi summaries.
- Functional interactions: a CSV of aggregated GtoPdb evidence by (gtp_ligand_id, target_uniprot).
- Metrics report: a timestamped .txt file with target match counts, ligand match counts, functional join coverage, and example ambiguous or unmatched records.

## Caveats and practical implications

- Target mapping is generally robust; most misses are symbol discrepancies or non human entries.
- Ligand mapping is the limiting factor because many PDSP entries lack structures or CAS identifiers.
- Name based matching is useful but risky and should be interpreted carefully.
- Functional direction is attached only where a ligand and target join successfully.
- Canonicalisation decisions (salt stripping and stereochemistry preserved) are chosen to balance match rate and correctness.

See the diagnostics report for specific figures for the above.
