"""
Extract WSI patches centered exactly at ST spot coordinates.

Critical correctness invariants (do not relax without re-verifying):
1. Spot centers are pulled from full-resolution pixel space (pxl_row_in_fullres /
   pxl_col_in_fullres, or obsm['spatial'] IF AND ONLY IF it is already full-res).
   If your coords are hires-scaled, they get divided by tissue_hires_scalef here.
2. Patch size is computed per-slide from that slide's own mpp/objective-power,
   so every patch covers the SAME physical field of view (default: 112 um) even
   though TCGA slides mix 20x and 40x scans.
3. Patches are read from openslide level 0 (native resolution) using read_region,
   then resized down to a fixed pixel size (default 224x224) for the encoder.

Before running on your real data: run verify_alignment() on 3-5 slides and
visually confirm the crops line up on tissue, not background/whitespace.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import openslide
except ImportError:
    openslide = None

try:
    import anndata as ad
except ImportError:
    ad = None


# ---------------------------------------------------------------------------
# Step 1: get spot centers in FULL-RES pixel space, unambiguously
# ---------------------------------------------------------------------------

def get_fullres_spot_coords(h5ad_path: str, scalefactors_path: str | None = None):
    """
    Returns a DataFrame [barcode, x_fullres, y_fullres] in WSI level-0 pixel space.

    Handles two common cases:
      (a) adata.obsm['spatial'] is already full-res (Space Ranger default when
          the h5ad was built directly from spaceranger's tissue_positions.csv)
      (b) adata.obsm['spatial'] is hires-scaled (some ST atlases store coords
          matched to tissue_hires_image.png) -> must divide by tissue_hires_scalef

    You MUST know which case applies to your specific data source. Check
    adata.uns['spatial'][sample_id]['scalefactors'] if present -- if it exists,
    trust it over guessing.
    """
    adata = ad.read_h5ad(h5ad_path)

    coords = adata.obsm["spatial"]  # (N, 2) -> typically [x, y] i.e. [col, row]
    barcodes = adata.obs_names.to_numpy()

    scale = 1.0
    source = "assumed full-res (no scalefactors found)"

    # Case: scalefactors embedded in the h5ad itself (typical scanpy/squidpy layout)
    if "spatial" in adata.uns:
        sample_keys = list(adata.uns["spatial"].keys())
        if sample_keys:
            sf = adata.uns["spatial"][sample_keys[0]].get("scalefactors", {})
            if "tissue_hires_scalef" in sf:
                # If coords were stored against the hires image, undo that scaling.
                # IMPORTANT: only do this if you've confirmed obsm['spatial'] is
                # actually in hires space -- some pipelines already store full-res
                # coords AND scalefactors side by side without coords needing scaling.
                # Default here: do NOT auto-scale; just surface the factor so you
                # can decide explicitly. See verify_alignment().
                source = f"scalefactors found: {sf}. NOT auto-applied -- verify first."

    # Case: separate scalefactors_json.json supplied explicitly
    if scalefactors_path:
        with open(scalefactors_path) as f:
            sf = json.load(f)
        scale = 1.0 / sf["tissue_hires_scalef"]
        source = f"external scalefactors_json.json, applied scale={scale:.4f}"

    df = pd.DataFrame({
        "barcode": barcodes,
        "x_fullres": coords[:, 0] * scale,
        "y_fullres": coords[:, 1] * scale,
    })
    print(f"[coords] {len(df)} spots loaded. Scale source: {source}")
    return df


# ---------------------------------------------------------------------------
# Step 2: compute per-slide patch size so physical FOV is constant
# ---------------------------------------------------------------------------

def get_patch_size_px(slide: "openslide.OpenSlide", fov_um: float = 112.0) -> int:
    """
    Returns the level-0 pixel size such that the extracted square patch covers
    `fov_um` microns per side, regardless of whether this slide was scanned at
    20x or 40x. This is what makes patches comparable ACROSS slides.
    """
    props = slide.properties
    mpp_x = props.get(openslide.PROPERTY_NAME_MPP_X)
    mpp_y = props.get(openslide.PROPERTY_NAME_MPP_Y)

    if mpp_x is None or mpp_y is None:
        # Fallback: infer from objective power (20x ~ 0.5 mpp, 40x ~ 0.25 mpp)
        obj_power = props.get(openslide.PROPERTY_NAME_OBJECTIVE_POWER)
        if obj_power is None:
            raise ValueError(
                "Slide has neither mpp nor objective-power metadata -- "
                "cannot safely compute physical patch size. Inspect manually."
            )
        obj_power = float(obj_power)
        mpp = 0.5 * (40.0 / obj_power)  # scanner-standard approximation
        print(f"  [warn] no mpp in metadata, inferred mpp={mpp:.4f} from "
              f"objective_power={obj_power}")
    else:
        mpp = (float(mpp_x) + float(mpp_y)) / 2.0

    patch_px = int(round(fov_um / mpp))
    return patch_px


# ---------------------------------------------------------------------------
# Step 3: extract patches
# ---------------------------------------------------------------------------

def extract_patches_for_slide(
    svs_path: str,
    coords_df: pd.DataFrame,
    out_dir: str,
    fov_um: float = 112.0,
    output_px: int = 224,
):
    """
    coords_df: DataFrame with columns [barcode, x_fullres, y_fullres] for THIS slide only
    Saves one .png per spot named <barcode>.png, plus a manifest.csv
    """
    slide = openslide.OpenSlide(svs_path)
    patch_px_native = get_patch_size_px(slide, fov_um=fov_um)
    half = patch_px_native // 2

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []

    w0, h0 = slide.level_dimensions[0]

    for _, row in coords_df.iterrows():
        x, y = int(row["x_fullres"]), int(row["y_fullres"])
        top_left = (x - half, y - half)

        # Guard against spots too close to the slide edge
        if top_left[0] < 0 or top_left[1] < 0 or \
           top_left[0] + patch_px_native > w0 or top_left[1] + patch_px_native > h0:
            manifest_rows.append({**row, "status": "skipped_edge"})
            continue

        region = slide.read_region(top_left, level=0, size=(patch_px_native, patch_px_native))
        region = region.convert("RGB").resize((output_px, output_px))

        out_path = out_dir / f"{row['barcode']}.png"
        region.save(out_path)
        manifest_rows.append({**row, "status": "ok", "patch_path": str(out_path),
                               "native_patch_px": patch_px_native})

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(out_dir / "manifest.csv", index=False)
    n_ok = (manifest["status"] == "ok").sum()
    print(f"[extract] {svs_path}: {n_ok}/{len(coords_df)} patches saved "
          f"(native_px={patch_px_native}, output_px={output_px})")
    return manifest


# ---------------------------------------------------------------------------
# Step 4: MANDATORY sanity check before trusting any of this
# ---------------------------------------------------------------------------

def verify_alignment(svs_path: str, coords_df: pd.DataFrame, out_path: str,
                      n_check: int = 12, fov_um: float = 112.0):
    """
    Draws the first n_check spot centers as circles on a downscaled thumbnail
    of the WSI, so you can visually confirm they land on tissue (not background,
    not systematically offset). RUN THIS FIRST on every new data source before
    extracting patches for the full dataset.
    """
    import PIL.Image
    import PIL.ImageDraw

    slide = openslide.OpenSlide(svs_path)
    w0, h0 = slide.level_dimensions[0]
    thumb_w = 1600
    scale = thumb_w / w0
    thumb = slide.get_thumbnail((thumb_w, int(h0 * scale)))
    draw = PIL.ImageDraw.Draw(thumb)

    for _, row in coords_df.head(n_check).iterrows():
        x, y = row["x_fullres"] * scale, row["y_fullres"] * scale
        r = 5
        draw.ellipse([x - r, y - r, x + r, y + r], outline="red", width=2)

    thumb.save(out_path)
    print(f"[verify] Saved alignment check to {out_path} -- open it and confirm "
          f"red dots sit on tissue, in a pattern matching the Visium grid.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--svs", required=True)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--scalefactors", default=None)
    parser.add_argument("--out", required=True)
    parser.add_argument("--fov_um", type=float, default=112.0)
    parser.add_argument("--output_px", type=int, default=224)
    parser.add_argument("--verify_only", action="store_true")
    args = parser.parse_args()

    coords = get_fullres_spot_coords(args.h5ad, args.scalefactors)

    if args.verify_only:
        verify_alignment(args.svs, coords, out_path=str(Path(args.out) / "alignment_check.png"),
                          fov_um=args.fov_um)
    else:
        extract_patches_for_slide(args.svs, coords, args.out, fov_um=args.fov_um,
                                   output_px=args.output_px)
