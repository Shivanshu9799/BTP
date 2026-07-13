"""
Exact DeepSpot-M tutorial patching (predict_tcga_skcm.ipynb):
  - Plain non-overlapping grid, stride = tile size = 224x224
  - Background filter: mean pixel value > white_mean (220.0) -> discard
  - No spot alignment, no coordinate matching
"""

import argparse
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
from PIL import Image

try:
    import openslide
    HAVE_OPENSLIDE = True
except ImportError:
    HAVE_OPENSLIDE = False

PATCH_SIZE = 224
WHITE_MEAN_THRESHOLD = 220.0


def passes_background_filter(patch_img, white_mean=WHITE_MEAN_THRESHOLD):
    arr = np.asarray(patch_img.convert("L"), dtype=np.float32)
    return arr.mean() <= white_mean


def extract_slide(svs_path, out_dir, patch_size=PATCH_SIZE,
                   white_mean=WHITE_MEAN_THRESHOLD, verify_only=False, limit=None):
    if not HAVE_OPENSLIDE:
        raise RuntimeError("openslide not available")

    slide = openslide.OpenSlide(svs_path)
    W, H = slide.dimensions
    print(f"[diag] slide dims: {W}x{H}")

    os.makedirs(out_dir, exist_ok=True)
    patches_dir = os.path.join(out_dir, "patches")
    os.makedirs(patches_dir, exist_ok=True)

    manifest_rows = []
    count = 0

    for top in range(0, H - patch_size + 1, patch_size):
        for left in range(0, W - patch_size + 1, patch_size):
            if limit is not None and count >= limit:
                break

            region = slide.read_region((left, top), 0, (patch_size, patch_size)).convert("RGB")

            if not passes_background_filter(region, white_mean):
                continue

            fname = f"x{left}_y{top}.png"
            if not verify_only:
                region.save(os.path.join(patches_dir, fname))
                manifest_rows.append({
                    "x": left, "y": top, "patch_size": patch_size, "file": fname,
                })
            count += 1
        if limit is not None and count >= limit:
            break

    if not verify_only:
        pd.DataFrame(manifest_rows).to_csv(os.path.join(out_dir, "manifest.csv"), index=False)
    print(f"[diag] kept {count} patches after background filter")

    slide.close()


def _worker(args_tuple):
    slide_id, svs_path, out_dir, patch_size, white_mean, verify_only, limit = args_tuple
    try:
        extract_slide(
            svs_path, out_dir,
            patch_size=patch_size, white_mean=white_mean,
            verify_only=verify_only, limit=limit,
        )
        return (slide_id, None)
    except Exception as e:
        return (slide_id, str(e))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report_csv", required=True, help="coad_ff_pairs.csv")
    ap.add_argument("--out_root", required=True)
    ap.add_argument("--patch_size", type=int, default=PATCH_SIZE)
    ap.add_argument("--white_mean", type=float, default=WHITE_MEAN_THRESHOLD)
    ap.add_argument("--limit", type=int, default=None, help="limit patches per slide (dry-run)")
    ap.add_argument("--slide_limit", type=int, default=None, help="limit number of slides (dry-run)")
    ap.add_argument("--verify_only", action="store_true")
    ap.add_argument("--workers", type=int, default=1, help="parallel CPU workers (processes)")
    args = ap.parse_args()

    pairs = pd.read_csv(args.report_csv)
    if args.slide_limit:
        pairs = pairs.head(args.slide_limit)

    jobs = []
    for _, row in pairs.iterrows():
        slide_id = row["barcode"]
        out_dir = os.path.join(args.out_root, slide_id)
        jobs.append((slide_id, row["svs_path"], out_dir,
                     args.patch_size, args.white_mean, args.verify_only, args.limit))

    if args.workers <= 1:
        for job in jobs:
            print(f"[diag] extracting {job[0]}")
            slide_id, err = _worker(job)
            if err:
                print(f"[error] {slide_id}: {err}")
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(_worker, job): job[0] for job in jobs}
            for fut in as_completed(futures):
                slide_id, err = fut.result()
                if err:
                    print(f"[error] {slide_id}: {err}")
                else:
                    print(f"[diag] done {slide_id}")


if __name__ == "__main__":
    main()
