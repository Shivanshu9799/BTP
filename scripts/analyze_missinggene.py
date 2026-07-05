import pandas as pd
from pathlib import Path

BASE = Path("/mnt/wwn-0x5000c500e64c5cc0/shivanhu")

gene_info = pd.read_csv(
    BASE / "HINGE/CellFM-main/csv/gene_info_hinge.csv"
)

cellfm = set(gene_info["gene_name"].astype(str))

selected = [
    g.strip()
    for g in open(
        BASE / "Data/hest1k_datasets/BREAST/processed_data/used_selected_gene_list.txt"
    )
]

missing = sorted([g for g in selected if g not in cellfm])

print(f"Missing genes: {len(missing)}")

# --------------------------
# HGNC alias database
# --------------------------

hgnc = pd.read_csv(
    BASE / "HINGE/CellFM-main/csv/updated_hgcn.tsv",
    sep="\t",
    dtype=str
).fillna("")

approved_col = next(
    c for c in hgnc.columns
    if "approved symbol" in c.lower() or c.lower() == "symbol"
)

alias = {}

for _, row in hgnc.iterrows():

    approved = row[approved_col]

    for col in hgnc.columns:

        if "previous" in col.lower() or "alias" in col.lower():

            vals = row[col]

            if vals:

                for x in vals.split(","):
                    alias[x.strip()] = approved

alias_hits = []
true_missing = []

for g in missing:

    if g in alias:
        alias_hits.append((g, alias[g]))
    else:
        true_missing.append(g)

print("\n========================")
print("Alias candidates")
print("========================")

for old, new in alias_hits:
    print(f"{old:25s} -> {new}")

print("\n========================")
print("Truly missing")
print("========================")

for g in true_missing:
    print(g)

print("\nSummary")
print("-------")
print("Missing:", len(missing))
print("Alias :", len(alias_hits))
print("True  :", len(true_missing))
