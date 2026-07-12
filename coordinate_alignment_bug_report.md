# ST Spot Coordinate Alignment Bug — Diagnosis & Fix

**Dataset:** TCGA-COAD Fresh Frozen (FF), WSI–ST matched pairs
**Script affected:** `extract_st_matched_patches.py`
**Status:** Fixed and verified on 3 pairs

---

## 1. The Issue

Patch extraction centers a crop on each ST spot's pixel coordinates. The
original coordinate-loading logic pulled coordinates from `adata.obsm['spatial']`,
assuming this was already in full-resolution WSI pixel space (the standard
scanpy/Space Ranger convention).

Running `--verify_only` on `TCGA-3L-AA1B` showed spot centers clustered in the
top-left corner of the slide, nowhere near the tissue:

![Before: broken alignment](images/01_before_broken_topleft.png)

This meant `obsm['spatial']` in this dataset was **not** full-resolution pixel
space — it was a much smaller-scale coordinate system.

---

## 2. Diagnosis

A diagnostic script inspected `adata.obs`, `adata.uns`, and cross-referenced
against the actual `.svs` dimensions via OpenSlide. Key findings:

- `adata.obs` contained explicit columns: `x_array`, `y_array`, `x_pixel`,
  `y_pixel`, `x_pixel_fullres`, `y_pixel_fullres`.
- `adata.uns` had **no** `spatial` key (i.e. no standard scanpy scalefactors),
  which is why the original code silently fell back to "assumed full-res."
- `obsm['spatial']` values turned out to be exactly `x_pixel_fullres / 10` and
  `y_pixel_fullres / 10` — a downsampled preview coordinate set, not the
  target scale.

### First fix attempt (partial)

Switching to `adata.obs['x_pixel_fullres']` / `adata.obs['y_pixel_fullres']`
moved the spot cloud dramatically — but it still didn't land on tissue. It
sat consistently offset to the left:

![After partial fix: still offset](images/02_after_partial_fix_offset.png)

### Root cause: axis swap

Adding a diagnostic comparing coordinate ranges against actual WSI dimensions
(`129480 × 33445` for this slide) revealed the real bug:

| Column | Range | Fits inside... |
|---|---|---|
| `x_pixel_fullres` | 7,161 – 30,041 (span ≈ 22,880) | slide **height** (33,445), not width |
| `y_pixel_fullres` | 10,681 – 118,921 (span ≈ 108,240) | slide **width** (129,480), not height |

This is the classic Space Ranger `pxl_row_in_fullres` / `pxl_col_in_fullres`
naming collapsed into generic `x_pixel_fullres` / `y_pixel_fullres` column
names by this atlas's export — but the underlying values kept row (→height)
and col (→width) semantics. So `x_pixel_fullres` was actually the **row**
(vertical) coordinate, and `y_pixel_fullres` was actually the **column**
(horizontal) coordinate — the opposite of what the names suggest.

---

## 3. The Fix

In `get_fullres_spot_coords()`, when `x_pixel_fullres`/`y_pixel_fullres` are
present, the columns are now swapped on load:

```python
df = pd.DataFrame({
    "barcode": barcodes,
    "x_fullres": adata.obs["y_pixel_fullres"].to_numpy(dtype=float),  # col -> x
    "y_fullres": adata.obs["x_pixel_fullres"].to_numpy(dtype=float),  # row -> y
})
```

Two supporting fixes were made alongside this:

1. **`verify_alignment()` now samples randomly**, not `coords_df.head(n)`.
   Taking the first N barcodes alphabetically is not a spatial sample — it
   can cluster near tissue edges and produce a misleading "still broken"
   read even when the fix is correct. It also now plots 300 points by
   default instead of 12, so the true Visium grid pattern is visible.
2. **A `[diag]` block** prints coordinate min/max/mean against actual WSI
   dimensions on every `--verify_only` run, so any future scale/axis issue
   is caught from console output alone, without needing a visual re-check.

---

## 4. Verification

### TCGA-3L-AA1B (primary fix target)

Diagnostic output after the fix:

```
[diag] WSI dims (level 0): 129480 x 33445
[diag] x_fullres: min=10681.0, max=118921.0, mean=67987.3
[diag] y_fullres: min=7161.0, max=30041.0, mean=15896.6
[diag] x range as fraction of WSI width: 0.8360
[diag] y range as fraction of WSI height: 0.6841
[diag] implied x scale factor (w0/x_range): 1.1962
[diag] implied y scale factor (h0/y_range): 1.4618
```

300 randomly-sampled spots now land cleanly on tissue, tracing the tumor
mass outline:

![Fixed: TCGA-3L-AA1B](images/03_fixed_TCGA-3L-AA1B.png)

### Generalization check — 2 additional pairs

The fix was re-verified on two more pairs chosen to differ from the first in
tissue type, portion, and serial-section position:

- **TCGA-A6-2672-11A-01-TS1** — normal tissue (sample type `11`, not tumor),
  Top Slide serial section.
- **TCGA-A6-2677-01B-01-BS1** — primary tumor, second portion (`01B`),
  Bottom Slide serial section, two-fragment specimen.

Both show spot centers landing correctly on tissue, including tracing
individual glandular/mucosal folds in the normal-tissue slide:

![Verified: TCGA-A6-2677-01B](images/04_verified_TCGA-A6-2677-01B.png)

![Verified: TCGA-A6-2672-11A](images/05_verified_TCGA-A6-2672-11A.png)

---

## 5. Conclusion

The axis-swap fix generalizes across tumor vs. normal tissue, single- vs.
multi-fragment slides, and different serial-section positions (Top/Bottom).
Fix confirmed correct on 3/3 tested pairs.

**Outstanding item before full batch extraction:** confirm whether
`batch_extract_patches.py` imports `get_fullres_spot_coords()` from
`extract_st_matched_patches.py` directly, or has a separate/duplicated copy
of the coordinate-loading logic that would also need this patch.
