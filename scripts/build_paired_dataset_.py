"""
Build a PAIRED dataset manifest: for every TCGA sample that has BOTH
a WSI slide and Spatial Transcriptomics data, list out the file(s) on
each side so you can directly use this as your training manifest.

Usage:
    python build_paired_dataset.py \
        /path/to/WSI/TCGA_COAD/FF \
        /path/to/ST/TCGA_COAD/FF \
        --out paired_manifest.csv

Output CSV columns:
    barcode          - sample-level TCGA barcode (e.g. TCGA-A6-2674-01A)
    wsi_files         - semicolon-separated list of all WSI filenames for this sample
    wsi_count         - how many WSI files exist for this sample
    st_files          - semicolon-separated list of all ST filenames for this sample
    st_count          - how many ST files exist for this sample
    wsi_primary       - the single WSI file picked as "primary" (see --prefer)
    st_primary        - the single ST file picked as "primary" (usually only 1)

Only samples present in BOTH sides are written to the paired manifest.
Unpaired samples (one side only) are written to a separate _unpaired.csv
for visibility, not silently dropped.
"""

import os
import re
import csv
import argparse

BARCODE_RE = re.compile(r'^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}(?:-\d{2}[A-Z](?:-\d{2}(?:-[A-Za-z0-9]+)?)?)?)')


def extract_barcode(filename, level="sample"):
    m = BARCODE_RE.match(filename)
    if not m:
        return None
    parts = m.group(1).split('-')
    if level == "case":
        return '-'.join(parts[0:3])
    elif level == "sample":
        return '-'.join(parts[0:4])
    else:
        return m.group(1)


def collect_barcodes(folder, level="sample"):
    barcodes = {}
    for root, dirs, files in os.walk(folder):
        for f in files:
            if f.startswith('.') or f.endswith('.log') or f.endswith('.partial') or f.endswith('.parcel'):
                continue
            if f in ('annotations.txt',) or f.startswith('gdc_manifest'):
                continue
            barcode = extract_barcode(f, level=level)
            if barcode is None:
                continue
            barcodes.setdefault(barcode, []).append(os.path.join(root, f))
    return barcodes


def pick_primary(files, prefer="TS"):
    """
    When a sample has multiple WSI files (e.g. BS1, BS2, TS1), pick one as
    'primary' for training. Default preference: prefer 'TS' (top slide) over
    'BS' (bottom slide) since TS is generally the standard diagnostic slide.
    Falls back to first file alphabetically if no preference match.
    """
    preferred = [f for f in files if f"-{prefer}" in os.path.basename(f).upper()]
    if preferred:
        return sorted(preferred)[0]
    return sorted(files)[0]


def main():
    parser = argparse.ArgumentParser(description="Build paired WSI+ST manifest.")
    parser.add_argument("wsi_folder", help="Folder containing WSI (.svs) files")
    parser.add_argument("st_folder", help="Folder containing ST (.h5ad.gz) files")
    parser.add_argument("--level", choices=["slide", "sample", "case"], default="sample")
    parser.add_argument("--prefer", default="TS", help="Preferred slide-type suffix when multiple WSI files exist per sample (default TS)")
    parser.add_argument("--out", default="paired_manifest.csv")
    parser.add_argument("--unpaired-out", default="unpaired_manifest.csv")
    args = parser.parse_args()

    print(f"Scanning WSI folder: {args.wsi_folder}")
    wsi = collect_barcodes(args.wsi_folder, level=args.level)
    print(f"  -> {sum(len(v) for v in wsi.values())} files -> {len(wsi)} unique samples")

    print(f"Scanning ST folder: {args.st_folder}")
    st = collect_barcodes(args.st_folder, level=args.level)
    print(f"  -> {sum(len(v) for v in st.values())} files -> {len(st)} unique samples")

    wsi_set, st_set = set(wsi.keys()), set(st.keys())
    paired = sorted(wsi_set & st_set)
    only_wsi = sorted(wsi_set - st_set)
    only_st = sorted(st_set - wsi_set)

    print(f"\nPaired samples (both WSI + ST): {len(paired)}")
    print(f"Unpaired (WSI only): {len(only_wsi)}")
    print(f"Unpaired (ST only):  {len(only_st)}")

    # --- write paired manifest: ONE ROW PER RAW WSI FILE (no dedup to one-per-sample) ---
    total_rows = 0
    with open(args.out, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["barcode", "wsi_file", "st_file"])
        for bc in paired:
            wsi_files = sorted(wsi[bc])
            st_files = sorted(st[bc])
            # normally 1 ST file per sample; if a sample has multiple WSI files,
            # every one of them gets paired with that same ST file (since the
            # expression profile belongs to the sample, not to a specific slide cut)
            st_primary = st_files[0]
            if len(st_files) > 1:
                print(f"  NOTE: {bc} has {len(st_files)} ST files, using {os.path.basename(st_primary)} for all its WSI rows")
            for wsi_f in wsi_files:
                writer.writerow([bc, os.path.basename(wsi_f), os.path.basename(st_primary)])
                total_rows += 1
    print(f"\nPaired manifest written: {args.out}  ({total_rows} rows -- one row per WSI file)")

    # --- write unpaired manifest (for visibility, not silently dropped) ---
    with open(args.unpaired_out, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["barcode", "missing_side", "available_files"])
        for bc in only_wsi:
            writer.writerow([bc, "ST", ';'.join(os.path.basename(x) for x in wsi[bc])])
        for bc in only_st:
            writer.writerow([bc, "WSI", ';'.join(os.path.basename(x) for x in st[bc])])
    print(f"Unpaired manifest written: {args.unpaired_out}  ({len(only_wsi) + len(only_st)} rows)")

    if only_wsi:
        print(f"\nSamples with WSI but no ST ({len(only_wsi)}):")
        for bc in only_wsi:
            print(f"  - {bc}")


if __name__ == "__main__":
    main()
