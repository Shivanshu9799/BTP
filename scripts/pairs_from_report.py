"""
Builds the svs_path / h5ad_path / barcode pairing CSV using the ALREADY-VERIFIED
matched barcode list from barcode_missing_report.txt, rather than re-deriving
matches independently. This report is ground truth (962 matched, confirmed by
your earlier pairing work) -- this script only resolves each barcode to its
actual file path on disk.

Usage:
    python pairs_from_report.py \
        --report /path/to/barcode_missing_report.txt \
        --wsi_root "/mnt/.../WSI/TCGA_COAD/FF" \
        --st_root "/mnt/.../TCGA_data/data/TCGA_COAD/FF" \
        --out_csv /home/shivanshu/coad_ff_pairs.csv
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_matched_barcodes(report_path: str) -> list[str]:
    barcodes = []
    in_matched_section = False
    with open(report_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("--- MATCHED BARCODES"):
                in_matched_section = True
                continue
            if in_matched_section:
                if not line or line.startswith("---"):
                    break
                barcodes.append(line)
    return barcodes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True)
    parser.add_argument("--wsi_root", required=True)
    parser.add_argument("--st_root", required=True)
    parser.add_argument("--out_csv", required=True)
    args = parser.parse_args()

    wsi_root = Path(args.wsi_root)
    st_root = Path(args.st_root)

    barcodes = parse_matched_barcodes(args.report)
    print(f"Parsed {len(barcodes)} matched barcodes from report")
    if len(barcodes) != 962:
        print(f"[warn] Expected 962, got {len(barcodes)} -- check report format "
              f"hasn't changed / parsing didn't stop early.")

    rows = []
    missing_svs = []
    missing_h5ad = []
    ambiguous = []

    for bc in barcodes:
        svs_matches = list(wsi_root.rglob(f"{bc}.*.svs"))
        h5ad_matches = list(st_root.glob(f"{bc}.*.h5ad.gz")) or \
                       list(st_root.glob(f"{bc}.*.gz")) or \
                       list(st_root.glob(f"{bc}.*.h5ad"))

        if not svs_matches:
            missing_svs.append(bc)
            continue
        if not h5ad_matches:
            missing_h5ad.append(bc)
            continue
        if len(svs_matches) > 1 or len(h5ad_matches) > 1:
            ambiguous.append((bc, svs_matches, h5ad_matches))

        rows.append({
            "barcode": bc,
            "svs_path": str(svs_matches[0]),
            "h5ad_path": str(h5ad_matches[0]),
        })

    print(f"\nResolved: {len(rows)}")
    print(f"Missing on disk despite being in 'matched' list -- no svs found: {len(missing_svs)}")
    print(f"Missing on disk despite being in 'matched' list -- no h5ad found: {len(missing_h5ad)}")
    print(f"Ambiguous (multiple files matched one barcode): {len(ambiguous)}")

    if missing_svs:
        print(f"\n[!] These barcodes are in the report as matched, but no .svs "
              f"found on disk right now -- files may have moved/been deleted "
              f"since the report was generated:")
        for bc in missing_svs[:10]:
            print(f"    {bc}")

    if missing_h5ad:
        print(f"\n[!] These barcodes are in the report as matched, but no h5ad "
              f"found on disk right now:")
        for bc in missing_h5ad[:10]:
            print(f"    {bc}")

    if ambiguous:
        print(f"\n[!] These barcodes matched multiple files -- using the first, "
              f"but you should check these manually:")
        for bc, svs_m, h5ad_m in ambiguous[:10]:
            print(f"    {bc}: {len(svs_m)} svs, {len(h5ad_m)} h5ad")

    df = pd.DataFrame(rows)
    df.to_csv(args.out_csv, index=False)
    print(f"\nWrote {len(df)} pairs to {args.out_csv}")


if __name__ == "__main__":
    main()
