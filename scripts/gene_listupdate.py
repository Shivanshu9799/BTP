#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

# -------------------------
# Paths
# -------------------------
GENE_LIST = Path("/mnt/wwn-0x5000c500e64c5cc0/shivanhu/Data/hest1k_datasets/BREAST/processed_data/used_selected_gene_list.txt")

HGNC = Path("/mnt/wwn-0x5000c500e64c5cc0/shivanhu/HINGE/CellFM-main/csv/updated_hgcn.tsv")

OUTPUT = GENE_LIST.parent / "used_selected_gene_list_updated.txt"

# -------------------------
# Load HGNC
# -------------------------
hgnc = pd.read_csv(HGNC, sep="\t", dtype=str).fillna("")

print("Columns:")
print(hgnc.columns.tolist())

# -------- automatically detect approved symbol column --------
approved_col = None

for c in hgnc.columns:
    lc = c.lower()
    if "approved symbol" in lc or lc == "symbol":
        approved_col = c
        break

if approved_col is None:
    raise RuntimeError("Could not detect approved symbol column.")

# -------- build alias dictionary --------

alias = {}

for _, row in hgnc.iterrows():

    approved = row[approved_col].strip()

    for col in hgnc.columns:

        name = col.lower()

        if "previous" in name or "alias" in name:

            val = row[col]

            if val:

                for g in str(val).split(","):
                    g = g.strip()

                    if g:
                        alias[g] = approved

print(f"Loaded {len(alias):,} aliases")

# -------------------------
# Convert gene list
# -------------------------

genes = [g.strip() for g in open(GENE_LIST)]

updated = []

changed = []

for g in genes:

    ng = alias.get(g, g)

    updated.append(ng)

    if ng != g:
        changed.append((g, ng))

# remove duplicates while preserving order

seen = set()
final = []

for g in updated:
    if g not in seen:
        final.append(g)
        seen.add(g)

with open(OUTPUT, "w") as f:
    for g in final:
        f.write(g + "\n")

print(f"\nOriginal genes : {len(genes)}")
print(f"Updated genes  : {len(final)}")
print(f"Renamed genes  : {len(changed)}")

print("\nFirst 20 changes:")

for old, new in changed[:20]:
    print(f"{old:20s} -> {new}")

print(f"\nSaved to:\n{OUTPUT}")
