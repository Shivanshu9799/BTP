#!/usr/bin/env python3
"""
Prepare HEST breast slides for HINGE with robust logging.

What this script does
---------------------
For each slide:
1) Load HEST patch H5: barcode, coords, img
2) Load raw slide h5ad
3) Treat patch H5 as the patch-eligible subset of spots
4) Subset/reorder adata to patch barcode order
5) Save filtered/aligned h5ad into HINGE st/
6) Extract UNI features
7) Extract CONCH features
8) Save .pt feature files in HINGE processed_data/

Output layout
-------------
<out_root>/
├── st/
│   └── <slide>.h5ad
└── processed_data/
    ├── all_slide_lst.txt
    ├── failed_slide_lst.txt
    ├── used_selected_gene_list.txt
    ├── 1spot_uni_ebd/
    │   └── <slide>_uni.pt
    └── 1spot_conch_ebd/
        └── <slide>_conch.pt
"""

import os
import sys
import json
import logging
import argparse
from pathlib import Path

import h5py
import numpy as np
import scanpy as sc
import torch
from PIL import Image
from tqdm import tqdm

logger = logging.getLogger("HEST_Prep")


def setup_logging(out_root):
    os.makedirs(out_root, exist_ok=True)
    log_file = os.path.join(out_root, "prepare_hinge.log")

    log_format = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    root = logging.getLogger()
    root.setLevel(logging.INFO)

    # Avoid duplicate handlers if re-run in same process
    if root.handlers:
        for h in list(root.handlers):
            root.removeHandler(h)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_format)
    root.addHandler(console_handler)

    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_format)
    root.addHandler(file_handler)

    logger.info(f"Logging initialized. Writing logs to: {log_file}")


def decode_barcode_array(arr):
    """Decode HEST barcode array to flat np.array[str]."""
    if hasattr(arr, "ndim") and arr.ndim > 1:
        arr = np.squeeze(arr)

    out = []
    if not isinstance(arr, (np.ndarray, list)):
        arr = [arr]

    for x in arr:
        if isinstance(x, np.ndarray):
            if x.size == 0:
                out.append("")
                continue
            x = x[0]
        if isinstance(x, bytes):
            x = x.decode("utf-8").strip("\x00")
        else:
            x = str(x).strip()
        out.append(x)
    return np.array(out, dtype=object)


def load_hest_patch_file(h5_path):
    """Load HEST patch H5 and return patch_barcodes, coords, imgs."""
    try:
        with h5py.File(h5_path, "r") as f:
            keys = list(f.keys())
            expected = {"barcode", "coords", "img"}
            missing = expected - set(keys)
            if missing:
                raise KeyError(f"Missing keys: {missing}. Found keys: {keys}")

            patch_barcodes = decode_barcode_array(f["barcode"][:])
            coords = f["coords"][:]
            imgs = f["img"][:]

        if len(patch_barcodes) != len(coords) or len(patch_barcodes) != len(imgs):
            raise ValueError(
                f"Inconsistent H5 lengths in {h5_path}: "
                f"barcodes={len(patch_barcodes)}, coords={len(coords)}, imgs={len(imgs)}"
            )

        return patch_barcodes, coords, imgs
    except Exception as e:
        logger.exception(f"Failed parsing H5 file at {h5_path}")
        raise e


def align_adata_to_patch_subset(adata, patch_barcodes, strict=True):
    """
    Treat patch_barcodes as source of truth for valid spots.
    Returns adata[patch_barcodes].copy()
    """
    adata_names = adata.obs_names.astype(str).values
    adata_set = set(adata_names)

    missing_in_adata = [bc for bc in patch_barcodes if bc not in adata_set]
    if len(missing_in_adata) > 0:
        preview = missing_in_adata[:10]
        raise ValueError(
            f"{len(missing_in_adata)} patch barcodes not found in adata.obs_names. "
            f"Examples: {preview}"
        )

    adata_aligned = adata[patch_barcodes].copy()

    if strict:
        aligned_names = adata_aligned.obs_names.astype(str).values
        assert len(aligned_names) == len(patch_barcodes)
        assert np.all(aligned_names == patch_barcodes), "adata order != patch barcode order"

    return adata_aligned


def infer_slide_ids_from_patch_dir(patch_dir):
    return sorted(p.stem for p in Path(patch_dir).glob("*.h5") if p.is_file())


def safe_torch_loadable_tensor(x):
    if isinstance(x, (tuple, list)):
        x = x[0]
    elif isinstance(x, dict):
        if "features" in x:
            x = x["features"]
        elif "image_features" in x:
            x = x["image_features"]
        else:
            raise ValueError(f"Unknown dict output keys: {list(x.keys())}")

    if not isinstance(x, torch.Tensor):
        raise TypeError(f"Expected torch.Tensor output, got {type(x)}")

    return x.detach().cpu()


def load_uni_model(device):
    """Load MahmoodLab UNI from Hugging Face using timm."""
    import timm
    from timm.data import resolve_data_config
    from timm.data.transforms_factory import create_transform

    logger.info("Initializing UNI ViT model via timm...")
    model = timm.create_model(
        "hf-hub:MahmoodLab/UNI",
        pretrained=True,
        init_values=1e-5,
        dynamic_img_size=True,
    )
    transform = create_transform(**resolve_data_config(model.pretrained_cfg, model=model))
    model = model.to(device).eval()
    return model, transform


def load_conch_model(device, hf_token=None, conch_ckpt=None):
    """
    Load MahmoodLab CONCH image encoder.

    Requires:
      pip install git+https://github.com/Mahmoodlab/CONCH.git
    """
    try:
        from conch.open_clip_custom import create_model_from_pretrained
    except ImportError as e:
        logger.error(
            "CONCH library missing. Install with:\n"
            "  pip install git+https://github.com/Mahmoodlab/CONCH.git"
        )
        raise ImportError("CONCH library missing.") from e

    logger.info("Initializing CONCH OpenCLIP model...")
    if conch_ckpt is not None:
        model, preprocess = create_model_from_pretrained("conch_ViT-B-16", conch_ckpt)
    else:
        if hf_token is None:
            logger.warning(
                "No explicit Hugging Face token provided for CONCH. "
                "Will try cached auth / huggingface-cli login."
            )
        model, preprocess = create_model_from_pretrained(
            "conch_ViT-B-16",
            "hf_hub:MahmoodLab/CONCH",
            hf_auth_token=hf_token,
        )

    model = model.to(device).eval()
    return model, preprocess


@torch.inference_mode()
def extract_uni_features(imgs, model, transform, device, batch_size=64):
    feats = []
    for start in tqdm(range(0, len(imgs), batch_size), desc="Extracting UNI Features", leave=False):
        batch_imgs = imgs[start:start + batch_size]
        tensors = [transform(Image.fromarray(arr)) for arr in batch_imgs]
        batch = torch.stack(tensors, dim=0).to(device, non_blocking=True)
        out = safe_torch_loadable_tensor(model(batch))
        feats.append(out)
    return torch.cat(feats, dim=0)


@torch.inference_mode()
def extract_conch_features(imgs, model, preprocess, device, batch_size=64):
    feats = []
    for start in tqdm(range(0, len(imgs), batch_size), desc="Extracting CONCH Features", leave=False):
        batch_imgs = imgs[start:start + batch_size]
        tensors = [preprocess(Image.fromarray(arr)) for arr in batch_imgs]
        batch = torch.stack(tensors, dim=0).to(device, non_blocking=True)
        out = model.encode_image(batch, proj_contrast=False, normalize=False)
        feats.append(safe_torch_loadable_tensor(out))
    return torch.cat(feats, dim=0)


def write_gene_list(gene_list, out_path):
    with open(out_path, "w") as f:
        for g in gene_list:
            f.write(str(g) + "\n")


def derive_gene_list_from_first_slide(h5ad_path):
    adata = sc.read_h5ad(h5ad_path)
    return adata.var_names.astype(str).tolist()


def process_one_slide(
    slide,
    raw_st_dir,
    patch_dir,
    out_st_dir,
    uni_out_dir,
    conch_out_dir,
    uni_model,
    uni_transform,
    conch_model,
    conch_preprocess,
    device,
    batch_size=64,
    overwrite=False,
):
    patch_path = os.path.join(patch_dir, f"{slide}.h5")
    raw_h5ad_path = os.path.join(raw_st_dir, f"{slide}.h5ad")

    if not os.path.exists(patch_path) or not os.path.exists(raw_h5ad_path):
        logger.error(f"Missing patch/h5ad for slide: {slide}")
        raise FileNotFoundError(f"Missing patch/h5ad for slide {slide}")

    out_h5ad_path = os.path.join(out_st_dir, f"{slide}.h5ad")
    uni_out_path = os.path.join(uni_out_dir, f"{slide}_uni.pt")
    conch_out_path = os.path.join(conch_out_dir, f"{slide}_conch.pt")

    if (
        (not overwrite)
        and os.path.exists(out_h5ad_path)
        and os.path.exists(uni_out_path)
        and os.path.exists(conch_out_path)
    ):
        logger.info(f"[SKIP] Cache hit for slide: {slide}")
        try:
            old_adata = sc.read_h5ad(out_h5ad_path, backed="r")
            n_spots = old_adata.n_obs
        except Exception:
            n_spots = None
        return {
            "slide": slide,
            "n_spots": n_spots,
            "uni_shape": "Cached",
            "conch_shape": "Cached",
            "skipped": True,
        }

    logger.info(f"Loading patch data: {slide}")
    patch_barcodes, coords, imgs = load_hest_patch_file(patch_path)

    logger.info(f"Loading raw h5ad: {slide}")
    adata = sc.read_h5ad(raw_h5ad_path)

    logger.info(f"Aligning adata to patch subset: {slide}")
    adata_aligned = align_adata_to_patch_subset(adata, patch_barcodes, strict=True)

    assert adata_aligned.n_obs == len(patch_barcodes)
    assert imgs.shape[0] == adata_aligned.n_obs
    assert np.all(adata_aligned.obs_names.astype(str).values == patch_barcodes)

    logger.info(f"Saving aligned h5ad: {out_h5ad_path}")
    adata_aligned.write_h5ad(out_h5ad_path)

    logger.info(f"Running UNI + CONCH feature extraction: {slide}")
    uni_feats = extract_uni_features(imgs, uni_model, uni_transform, device, batch_size)
    conch_feats = extract_conch_features(imgs, conch_model, conch_preprocess, device, batch_size)

    assert uni_feats.shape[0] == adata_aligned.n_obs
    assert conch_feats.shape[0] == adata_aligned.n_obs

    logger.info(
        f"Saving features for {slide} | spots={adata_aligned.n_obs} | "
        f"UNI={tuple(uni_feats.shape)} | CONCH={tuple(conch_feats.shape)}"
    )
    torch.save(uni_feats, uni_out_path)
    torch.save(conch_feats, conch_out_path)

    meta = {
        "slide": slide,
        "n_spots": int(adata_aligned.n_obs),
        "uni_shape": list(uni_feats.shape),
        "conch_shape": list(conch_feats.shape),
    }
    meta_path = os.path.join(out_st_dir, f"{slide}_prep_meta.json")
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    return {
        "slide": slide,
        "n_spots": int(adata_aligned.n_obs),
        "uni_shape": tuple(uni_feats.shape),
        "conch_shape": tuple(conch_feats.shape),
        "skipped": False,
    }


def parse_args():
    p = argparse.ArgumentParser(description="Prepare HEST breast slides for HINGE with UNI + CONCH features.")
    p.add_argument("--patch_dir", type=str, required=True, help="Directory containing HEST patch .h5 files")
    p.add_argument("--raw_st_dir", type=str, required=True, help="Directory containing raw source .h5ad files")
    p.add_argument("--out_root", type=str, required=True, help="HINGE dataset root, e.g. Data/hest1k_datasets/BREAST")
    p.add_argument("--slides", type=str, nargs="*", default=None, help="Optional explicit slide IDs")
    p.add_argument("--batch_size", type=int, default=64)
    p.add_argument("--device", type=str, default=None, help="cuda / cuda:0 / cpu")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--gene_list_path", type=str, default=None, help="Optional one-gene-per-line file")
    p.add_argument("--hf_token", type=str, default=None, help="HF token for gated UNI/CONCH if needed")
    p.add_argument("--conch_ckpt", type=str, default=None, help="Optional local CONCH checkpoint path")
    return p.parse_args()


def main():
    args = parse_args()
    setup_logging(args.out_root)

    out_st_dir = os.path.join(args.out_root, "st")
    processed_dir = os.path.join(args.out_root, "processed_data")
    uni_out_dir = os.path.join(processed_dir, "1spot_uni_ebd")
    conch_out_dir = os.path.join(processed_dir, "1spot_conch_ebd")

    os.makedirs(out_st_dir, exist_ok=True)
    os.makedirs(processed_dir, exist_ok=True)
    os.makedirs(uni_out_dir, exist_ok=True)
    os.makedirs(conch_out_dir, exist_ok=True)

    device = args.device if args.device is not None else ("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"Using device: {device}")

    slides = list(args.slides) if args.slides else infer_slide_ids_from_patch_dir(args.patch_dir)
    if len(slides) == 0:
        logger.error("No slides found in patch_dir.")
        raise ValueError("No slides found to process.")

    logger.info(f"Resolved {len(slides)} slide(s): {slides}")

    logger.info("Loading UNI model...")
    uni_model, uni_transform = load_uni_model(device=device)

    logger.info("Loading CONCH model...")
    conch_model, conch_preprocess = load_conch_model(
        device=device,
        hf_token=args.hf_token,
        conch_ckpt=args.conch_ckpt,
    )

    summary = []
    for slide in slides:
        logger.info(f"--- Processing slide: {slide} ---")
        try:
            result = process_one_slide(
                slide=slide,
                raw_st_dir=args.raw_st_dir,
                patch_dir=args.patch_dir,
                out_st_dir=out_st_dir,
                uni_out_dir=uni_out_dir,
                conch_out_dir=conch_out_dir,
                uni_model=uni_model,
                uni_transform=uni_transform,
                conch_model=conch_model,
                conch_preprocess=conch_preprocess,
                device=device,
                batch_size=args.batch_size,
                overwrite=args.overwrite,
            )
            summary.append(result)
        except Exception:
            logger.exception(f"Unhandled error while processing slide: {slide}")

    processed_slides = [x["slide"] for x in summary]
    failed_slides = [s for s in slides if s not in processed_slides]

    all_slide_path = os.path.join(processed_dir, "all_slide_lst.txt")
    with open(all_slide_path, "w") as f:
        for s in processed_slides:
            f.write(s + "\n")
    logger.info(f"Wrote processed slide list: {all_slide_path}")

    failed_path = os.path.join(processed_dir, "failed_slide_lst.txt")
    with open(failed_path, "w") as f:
        for s in failed_slides:
            f.write(s + "\n")
    logger.info(f"Wrote failed slide list: {failed_path}")

    gene_out = os.path.join(processed_dir, "used_selected_gene_list.txt")
    if args.gene_list_path is not None:
        with open(args.gene_list_path, "r") as f:
            genes = [line.strip() for line in f if line.strip()]
        write_gene_list(genes, gene_out)
        logger.info(f"Wrote provided gene list: {gene_out} ({len(genes)} genes)")
    else:
        if processed_slides:
            first_slide = processed_slides[0]
            first_h5ad = os.path.join(out_st_dir, f"{first_slide}.h5ad")
            if os.path.exists(first_h5ad):
                genes = derive_gene_list_from_first_slide(first_h5ad)
                write_gene_list(genes, gene_out)
                logger.info(f"Wrote gene list from first slide var_names: {gene_out} ({len(genes)} genes)")
            else:
                logger.error(f"First processed h5ad missing: {first_h5ad}")
        else:
            logger.warning("No slides processed successfully; skipping gene list derivation.")

    logger.info(f"Processed slides: {len(processed_slides)}")
    logger.info(f"Failed slides: {len(failed_slides)}")
    if failed_slides:
        logger.warning(f"Failed slide IDs: {failed_slides}")

    logger.info("=== Final Summary ===")
    for row in summary:
        logger.info(row)

    logger.info("Done.")


if __name__ == "__main__":
    main()
