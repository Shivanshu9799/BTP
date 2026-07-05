import h5py
import scanpy as sc
import numpy as np
import sys

slide = sys.argv[1]   # e.g. SPA71
patch_path = f"HEST_BREAST/patches/{slide}.h5"
h5ad_path = f"Data/hest1k_datasets/BREAST/st/{slide}.h5ad"

with h5py.File(patch_path, "r") as f:
    patch_barcodes = f["barcode"][:]
    patch_coords = f["coords"][:]

# decode barcode
decoded = []
for x in patch_barcodes:
    v = x[0] if isinstance(x, np.ndarray) else x
    if isinstance(v, bytes):
        v = v.decode("utf-8")
    else:
        v = str(v)
    decoded.append(v)

patch_barcodes = np.array(decoded)

adata = sc.read_h5ad(h5ad_path)

print("slide:", slide)
print("patch spots:", len(patch_barcodes))
print("adata spots:", adata.n_obs)

print("\nFirst 10 patch barcodes:")
print(patch_barcodes[:10])

print("\nFirst 10 adata.obs_names:")
print(adata.obs_names[:10].tolist())

common = set(patch_barcodes).intersection(set(adata.obs_names.astype(str)))
print("\nCommon barcodes with adata.obs_names:", len(common))

print("\nadata.obs columns:")
print(list(adata.obs.columns))
