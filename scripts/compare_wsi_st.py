"""
Compare WSI (.svs) and Spatial Transcriptomics (.h5ad.gz) folders for the
SAME TCGA-COAD cohort by matching on TCGA BARCODE, not literal filename.

Why this matters:
  WSI file:  TCGA-A6-2671-01A-01-TS1.79a30fac-....svs
  ST file:   TCGA-A6-2671-01A-01-TS1.9c415218-....h5ad.gz
  -> different extension AND different UUID, but SAME sample.

This script strips the extension + UUID + parcel suffix and matches on the
TCGA barcode portion (e.g. "TCGA-A6-2671-01A-01-TS1"), so it correctly
identifies which samples have BOTH a slide and ST data, and which are
missing one side.

Usage:
    python compare_wsi_st.py /path/to/WSI/TCGA_COAD/FF /path/to/ST/TCGA_COAD/FF

Matching level (choose with --level):
    slide   (default) -> TCGA-A6-2671-01A-01-TS1   (exact slide/portion match)
    sample             -> TCGA-A6-2671-01A          (sample-vial level — use this
                           if one ST file can correspond to multiple slide files,
                           which is common: BS1 + TS1 slides but one ST profile)
    case               -> TCGA-A6-2671              (patient-level, coarsest)
"""

import os
import re
import sys
import argparse

BARCODE_RE = re.compile(r'^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}(?:-\d{2}[A-Z](?:-\d{2}(?:-[A-Za-z0-9]+)?)?)?)')


def extract_barcode(filename, level="slide"):
    """
    Pull the TCGA barcode out of a filename regardless of extension/UUID.
    filename examples handled:
      TCGA-A6-2671-01A-01-TS1.79a30fac-81ad-4c1f-8a1f-5e649afe275d.h5ad.gz
      TCGA-A6-2671-01A-01-BS1.aab0019c-....svs
      TCGA-A6-2671-01A-01-BS1.aab0019c-....svs.parcel
    """
    m = BARCODE_RE.match(filename)
    if not m:
        return None
    full_barcode = m.group(1)
    parts = full_barcode.split('-')

    if level == "case":
        return '-'.join(parts[0:3])          # TCGA-A6-2671
    elif level == "sample":
        return '-'.join(parts[0:4])          # TCGA-A6-2671-01A
    else:  # slide (most specific)
        return full_barcode                  # TCGA-A6-2671-01A-01-TS1


def collect_barcodes(folder, level="slide"):
    """Walk a folder, extract barcodes from every relevant file, skip .parcel/log duplicates."""
    barcodes = {}  # barcode -> list of matching filenames (for traceability)
    skipped = []

    for root, dirs, files in os.walk(folder):
        for f in files:
            if f.startswith('.') or f.endswith('.log') or f.endswith('.partial'):
                continue
            if f.endswith('.parcel'):
                continue  # .parcel is a duplicate marker of the .svs file, skip to avoid double count
            if f in ('annotations.txt',) or f.startswith('gdc_manifest'):
                continue  # metadata/manifest files sitting in the folder, not actual data

            barcode = extract_barcode(f, level=level)
            if barcode is None:
                skipped.append(f)
                continue
            barcodes.setdefault(barcode, []).append(f)

    return barcodes, skipped


def main():
    parser = argparse.ArgumentParser(description="Compare WSI and ST folders by TCGA barcode.")
    parser.add_argument("folder_a", help="First folder (e.g. WSI slides, .svs)")
    parser.add_argument("folder_b", help="Second folder (e.g. Spatial Transcriptomics, .h5ad.gz)")
    parser.add_argument("--label-a", default="A", help="Short label for folder_a in the report (e.g. WSI)")
    parser.add_argument("--label-b", default="B", help="Short label for folder_b in the report (e.g. ST)")
    parser.add_argument("--level", choices=["slide", "sample", "case"], default="sample",
                         help="Barcode matching granularity (default: sample)")
    parser.add_argument("--out", default="barcode_missing_report.txt", help="Output report file")
    args = parser.parse_args()

    print(f"Scanning {args.label_a}: {args.folder_a}")
    barcodes_a, skipped_a = collect_barcodes(args.folder_a, level=args.level)
    print(f"  -> {len(barcodes_a)} unique barcodes ({args.level}-level)")
    if skipped_a:
        print(f"  -> WARNING: {len(skipped_a)} files didn't match TCGA barcode pattern, skipped")

    print(f"Scanning {args.label_b}: {args.folder_b}")
    barcodes_b, skipped_b = collect_barcodes(args.folder_b, level=args.level)
    print(f"  -> {len(barcodes_b)} unique barcodes ({args.level}-level)")
    if skipped_b:
        print(f"  -> WARNING: {len(skipped_b)} files didn't match TCGA barcode pattern, skipped")

    set_a = set(barcodes_a.keys())
    set_b = set(barcodes_b.keys())

    only_in_a = sorted(set_a - set_b)   # has WSI, missing ST (or vice versa depending on arg order)
    only_in_b = sorted(set_b - set_a)
    matched = sorted(set_a & set_b)

    print(f"\n{'='*60}")
    print(f"MATCHED (present in both)        : {len(matched)}")
    print(f"Only in {args.label_a:<10}              : {len(only_in_a)}")
    print(f"Only in {args.label_b:<10}              : {len(only_in_b)}")
    print(f"{'='*60}\n")

    with open(args.out, 'w') as f:
        f.write(f"{args.label_a}: {args.folder_a}\n")
        f.write(f"{args.label_b}: {args.folder_b}\n")
        f.write(f"Matching level: {args.level}\n")
        f.write(f"Matched: {len(matched)} | Only-{args.label_a}: {len(only_in_a)} | Only-{args.label_b}: {len(only_in_b)}\n\n")

        f.write(f"--- ONLY IN {args.label_a} (missing from {args.label_b}) ({len(only_in_a)}) ---\n")
        for bc in only_in_a:
            f.write(f"{bc}\t<- {', '.join(barcodes_a[bc])}\n")

        f.write(f"\n--- ONLY IN {args.label_b} (missing from {args.label_a}) ({len(only_in_b)}) ---\n")
        for bc in only_in_b:
            f.write(f"{bc}\t<- {', '.join(barcodes_b[bc])}\n")

        f.write(f"\n--- MATCHED BARCODES ({len(matched)}) ---\n")
        for bc in matched:
            f.write(f"{bc}\n")

    print(f"Full report: {args.out}")

    if only_in_a:
        print(f"\nFirst 10 barcodes only in {args.label_a}:")
        for bc in only_in_a[:10]:
            print(f"  - {bc}")

    if only_in_b:
        print(f"\nFirst 10 barcodes only in {args.label_b}:")
        for bc in only_in_b[:10]:
            print(f"  - {bc}")


if __name__ == "__main__":
    main()
