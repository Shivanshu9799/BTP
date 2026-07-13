#!/usr/bin/env python3
"""
check_spot_spacing.py

Quick standalone check: for one h5ad (ST) file, report the median nearest-
neighbor spot spacing in pixels (and in microns if the paired WSI is given),
so you can see exactly how much overlap a given --patch_px will produce
before running the full batch.

Usage:
  python check_spot_spacing.py --h5ad path/to/slide.h5ad.gz
  python check_spot_spacing.py --h5ad path/to/slide.h5ad.gz --svs path/to/slide.svs --patch_px 224
"""
import argparse
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def is_hdf5(path):
    with open(path, "rb") as f:
        return f.read(8) == b"\x89HDF\r\n\x1a\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--svs", default=None, help="optional, to also report mpp and micron spacing")
    ap.add_argument("--patch_px", type=int, default=None, help="optional, report overlap % for this patch size")
    args = ap.parse_args()

    if not is_hdf5(args.h5ad):
        raise ValueError(f"{args.h5ad} is not valid HDF5 despite its name")

    import anndata as ad
    adata = ad.read_h5ad(args.h5ad)
    obs = adata.obs

    if {"x_pixel_fullres", "y_pixel_fullres"}.issubset(obs.columns):
        x = obs["x_pixel_fullres"].to_numpy(dtype=float)
        y = obs["y_pixel_fullres"].to_numpy(dtype=float)
        source = "obs.x/y_pixel_fullres"
    elif "spatial" in adata.obsm:
        sp = np.asarray(adata.obsm["spatial"], dtype=float)
        x, y = sp[:, 0], sp[:, 1]
        source = "obsm['spatial'] (raw, no scalefactor applied here)"
    else:
        raise ValueError("no usable spatial coordinates found in this h5ad")

    print(f"n_spots        : {len(x)}")
    print(f"coord source   : {source}")

    xy = np.stack([x, y], axis=1)
    tree = cKDTree(xy)
    dists, _ = tree.query(xy, k=2)
    nn = dists[:, 1]

    print(f"NN spacing px  : median={np.median(nn):.1f}  mean={nn.mean():.1f}  "
          f"min={nn.min():.1f}  max={nn.max():.1f}  std={nn.std():.1f}")

    if args.svs:
        import openslide
        slide = openslide.OpenSlide(args.svs)
        mpp = slide.properties.get("openslide.mpp-x")
        if mpp:
            mpp = float(mpp)
            print(f"slide mpp-x    : {mpp:.4f}")
            print(f"NN spacing um  : median={np.median(nn) * mpp:.1f}um")
        else:
            print("slide has no openslide.mpp-x property")

    if args.patch_px:
        med = np.median(nn)
        if args.patch_px > med:
            pct = 100 * (args.patch_px - med) / args.patch_px
            print(f"--patch_px {args.patch_px} -> ~{pct:.1f}% overlap between adjacent patches")
        else:
            pct = 100 * (med - args.patch_px) / med
            print(f"--patch_px {args.patch_px} -> no overlap, ~{pct:.1f}% gap between adjacent patches")


if __name__ == "__main__":
    main()
