"""
Batch-runs extract_patches_for_slide() over every row of your pairing manifest.

Expects a CSV with at minimum:
    svs_path, h5ad_path   (one row per WSI-ST pair)
Adjust COL_SVS / COL_H5AD below if your pairing file uses different column names.

Usage:
    python batch_extract_patches.py \
        --pairing_csv /path/to/your_962_pairs.csv \
        --out_root /home/shivanshu/patches \
        --fov_um 112 --output_px 224
"""

import argparse
from pathlib import Path

import pandas as pd

from extract_st_matched_patches import get_fullres_spot_coords, extract_patches_for_slide

COL_SVS = "svs_path"
COL_H5AD = "h5ad_path"
COL_BARCODE = "barcode"  # sample-level barcode, used to name the output subfolder


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairing_csv", required=True)
    parser.add_argument("--out_root", required=True)
    parser.add_argument("--fov_um", type=float, default=112.0)
    parser.add_argument("--output_px", type=int, default=224)
    parser.add_argument("--scalefactors", default=None,
                         help="Set this only if verify_alignment showed you need it")
    parser.add_argument("--limit", type=int, default=None,
                         help="Process only first N pairs, for a quick test run")
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairing_csv)
    if args.limit:
        pairs = pairs.head(args.limit)

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    failures = []
    for i, row in pairs.iterrows():
        svs_path = row[COL_SVS]
        h5ad_path = row[COL_H5AD]
        sample_id = row[COL_BARCODE] if COL_BARCODE in pairs.columns else Path(svs_path).stem

        print(f"\n[{i+1}/{len(pairs)}] {sample_id}")
        try:
            coords = get_fullres_spot_coords(h5ad_path, args.scalefactors)
            out_dir = out_root / sample_id
            extract_patches_for_slide(
                svs_path, coords, str(out_dir),
                fov_um=args.fov_um, output_px=args.output_px,
            )
        except Exception as e:
            print(f"  [FAILED] {sample_id}: {e}")
            failures.append({"sample_id": sample_id, "svs_path": svs_path,
                              "h5ad_path": h5ad_path, "error": str(e)})

    if failures:
        fail_df = pd.DataFrame(failures)
        fail_path = out_root / "extraction_failures.csv"
        fail_df.to_csv(fail_path, index=False)
        print(f"\n{len(failures)} slides failed. See {fail_path}")
    else:
        print("\nAll slides processed successfully.")


if __name__ == "__main__":
    main()
