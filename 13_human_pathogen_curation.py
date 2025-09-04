##### 13_Human_Pathogen_Curation
##### By Carter Clinton



#!/usr/bin/env python3
# human_pathogen_curation.py — Curate & flag human-associated/pathogenic taxa from QIIME taxonomy outputs
# Maps original file: Hu_microbiome_pathogen_list
# Analysis Stage: Human-Associated / Pathogen Curation
# Language: Python

"""
Overview
--------
Flags taxa in your QIIME taxonomy output that are documented as human-associated
or pathogenic, by cross-referencing one or more **curated reference lists** you provide.

Key design choices
- **No online calls, no API keys.** All inputs are local files you control.
- Flexible inputs: supports `taxonomy.tsv` (from `qiime tools export`) or species-only TSVs
  like those emitted by `species_taxonomy_export.sh`.
- Robust normalization of names: strips rank prefixes (e.g., `g__`, `s__`), common
  brackets/qualifiers, and compares at **species and/or genus** level.
- Transparent outputs: per-feature matches + species summaries + an audit report.
"""

from __future__ import annotations
import argparse, csv, json, re, sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

try:
    import pandas as pd
except ImportError as e:
    sys.stderr.write("ERROR: pandas is required. Try: pip install pandas\n")
    raise

# ---------------------------
# Normalization helpers
# ---------------------------

RANK_PREFIX_RE = re.compile(r"(^|\s)([dkpcofgs]__)")
SQUARE_BRACKETS_RE = re.compile(r"[\[\]]")
UNCULTURED_RE = re.compile(r"\buncultured\b|\bmetagenome\b|\bsp\.\b|\bincertae\s+sedis\b", re.I)

def normalize_taxon_name(name: str, strip_uncultured: bool = True) -> str:
    """Strip QIIME-style rank prefixes, square brackets, and common qualifiers."""
    if name is None:
        return ""
    n = str(name).strip()
    n = SQUARE_BRACKETS_RE.sub("", n)
    # keep the raw for decision; then strip rank tokens like 'g__', 's__'
    n = RANK_PREFIX_RE.sub(" ", n)
    n = re.sub(r"\s+", " ", n).strip()
    if strip_uncultured:
        n = UNCULTURED_RE.sub("", n)
        n = re.sub(r"\s+", " ", n).strip()
    return n

def extract_rank_from_taxon_string(taxon: str, rank_letter: str) -> Optional[str]:
    """
    From a hierarchical string like 'k__Bacteria; p__Firmicutes; g__Bacillus; s__B. subtilis',
    extract the piece that starts with '<rank_letter>__'. Returns content after '__'.
    """
    if not taxon:
        return None
    parts = [p.strip() for p in re.split(r";\s*", taxon)]
    for p in parts[::-1]:  # prefer the *last* occurrence at a rank
        if p.lower().startswith(f"{rank_letter.lower()}__"):
            return p.split("__", 1)[1].strip()
    return None

# ---------------------------
# I/O helpers
# ---------------------------

def read_taxonomy_table(path: Path) -> pd.DataFrame:
    """
    Accepts either:
      (A) taxonomy.tsv exported by QIIME (`Feature ID`, `Taxon`, `Confidence`)
      (B) species-only TSV from our export script (`FeatureID`, `Species`, `Taxon`, `Confidence`)
    """
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    cols = {c.lower().strip(): c for c in df.columns}
    # Standardize column keys
    if "feature id" in cols and "taxon" in cols:
        df = df.rename(columns={cols["feature id"]:"FeatureID", cols["taxon"]:"Taxon"})
        confcol = next((cols[k] for k in cols if k.startswith("confidence")), None)
        if confcol: df = df.rename(columns={confcol:"Confidence"})
        df["Species"] = df["Taxon"].apply(lambda t: extract_rank_from_taxon_string(t, "s") or "")
    elif "featureid" in cols and "taxon" in cols:
        # already species export format (or close)
        pass
    else:
        raise ValueError(f"Unrecognized taxonomy columns in {path}. "
                         "Expect `Feature ID`/`Taxon` or `FeatureID`/`Taxon`.")
    # Normalize convenience columns
    df["Species_norm"] = df["Species"].map(normalize_taxon_name)
    df["Genus"] = df["Taxon"].apply(lambda t: extract_rank_from_taxon_string(t, "g") or "")
    df["Genus_norm"] = df["Genus"].map(normalize_taxon_name)
    return df

def _read_ref_one(path: Path) -> pd.DataFrame:
    """
    Read a curated reference list (CSV/TSV). Expected columns (case-insensitive):
      - taxon  (required)  e.g., 'Staphylococcus aureus' OR 'Staphylococcus'
      - rank   (optional)  one of {'species','genus'} (default inferred by whitespace)
      - risk|risk_level (optional)  free text or a short code
      - source|evidence_source|reference (optional)  free text
      - notes  (optional)
    """
    # Infer delimiter
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    lower = {c.lower(): c for c in df.columns}
    if "taxon" not in lower:
        raise ValueError(f"{path} must contain a 'taxon' column.")
    df = df.rename(columns={lower["taxon"]: "taxon"})
    if "rank" not in lower:
        # Heuristic: has a space -> species; else genus
        df["rank"] = df["taxon"].map(lambda x: "species" if " " in str(x).strip() else "genus")
    else:
        df = df.rename(columns={lower["rank"]:"rank"})
        df["rank"] = df["rank"].str.lower().str.strip().replace({"sp":"species"})
    # optional columns
    for k in ("risk_level", "risk", "source", "evidence_source", "reference", "notes"):
        if k in lower:
            df = df.rename(columns={lower[k]: k})
    # Normalize
    df["taxon_norm"] = df["taxon"].map(normalize_taxon_name)
    df["rank"] = df["rank"].str.lower().str.strip()
    return df

def read_reference_lists(paths: List[Path]) -> pd.DataFrame:
    refs = []
    for p in paths:
        refs.append(_read_ref_one(p))
    if not refs:
        raise ValueError("At least one curated reference list is required.")
    ref = pd.con



