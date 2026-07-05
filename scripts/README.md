# HEST Spatial Transcriptomics Scripts

This directory contains Python scripts for processing and analyzing spatial transcriptomics data from the HEST (Histology-driven Explainable Spatial Transcriptomics) dataset, with a focus on breast tissue samples.

## Overview

These scripts handle the complete pipeline for preparing HEST breast slides for the HINGE framework, including:
- Data download and alignment
- Spatial transcriptomics preprocessing
- Feature extraction using deep learning models (UNI, CONCH)
- Gene list management and validation

---

## Scripts

### 1. **`data_preprocessing.py`** ⭐ (Main Pipeline)
The primary preprocessing script that orchestrates the entire data preparation workflow.

**What it does:**
- Loads HEST patch H5 files (barcodes, coordinates, images)
- Aligns raw slide h5ad files with patch-eligible spots
- Saves filtered/aligned h5ad files
- Extracts UNI vision transformer features
- Extracts CONCH image encoder features
- Generates processed output and metadata

**Output structure:**
```
<out_root>/
├── st/
│   └── <slide>.h5ad              # Aligned spatial transcriptomics data
├── processed_data/
│   ├── all_slide_lst.txt         # List of successfully processed slides
│   ├── failed_slide_lst.txt      # List of failed slides
│   ├── used_selected_gene_list.txt
│   ├── 1spot_uni_ebd/            # UNI features
│   │   └── <slide>_uni.pt
│   └── 1spot_conch_ebd/          # CONCH features
│       └── <slide>_conch.pt
```

**Usage:**
```bash
python data_preprocessing.py \
  --patch_dir HEST_BREAST/patches \
  --raw_st_dir Data/hest1k_datasets/BREAST/raw \
  --out_root Data/hest1k_datasets/BREAST \
  --batch_size 64 \
  --device cuda
```

---

### 2. **`dowload_hest_data.py`**
Downloads HEST breast tissue samples from the Hugging Face Hub.

**What it does:**
- Loads HEST metadata (v1.3.0)
- Filters samples by organ type (Breast)
- Downloads breast-specific data using Hugging Face Hub API

**Requirements:**
- Hugging Face account credentials (via `huggingface-cli login`)
- ~500GB+ storage for full breast dataset

**Usage:**
```bash
python dowload_hest_data.py
```

---

### 3. **`analyze_missinggene.py`**
Validates gene lists against reference databases and identifies missing genes.

**What it does:**
- Compares selected genes with CellFM gene database
- Cross-references with HGNC (HUGO Gene Nomenclature Committee) alias database
- Identifies gene name aliases for standardization
- Reports missing, aliased, and truly missing genes

**Output:**
- Gene categorization: Alias candidates vs. truly missing genes
- Summary statistics

**Usage:**
```bash
python analyze_missinggene.py
```

---

### 4. **`gene_listupdate.py`**
Updates and standardizes gene names using HGNC database.

**What it does:**
- Loads gene list from preprocessing output
- Maps outdated/alternate gene names to current approved symbols
- Removes duplicate entries while preserving order
- Outputs standardized gene list

**Usage:**
```bash
python gene_listupdate.py
```

---

### 5. **`check_hest_alignment.py`**
Validates alignment between patch H5 files and h5ad spatial data.

**What it does:**
- Loads patch barcodes from H5 file
- Loads spatial transcriptomics h5ad file
- Compares barcode matching between patch and h5ad
- Reports alignment statistics and metadata

**Usage:**
```bash
python check_hest_alignment.py <slide_id>  # e.g., SPA71
```

---

### 6. **`inspection.py`** (Utility)
Quick diagnostic script to inspect H5 file structure.

**What it does:**
- Reads H5 file keys and metadata
- Prints dataset shapes and data types

**Usage:**
```bash
python inspection.py
```

---

### 7. **`uni.py`** (Model Loader)
Standalone script to load the UNI vision transformer model.

**What it does:**
- Initializes MahmoodLab UNI pretrained model via TIMM
- Configures image preprocessing transforms

**Note:** This is typically called as part of `data_preprocessing.py`

**Usage:**
```bash
python uni.py
```

---

## Key Dependencies

```
pandas
numpy
scanpy
torch
h5py
PIL (Pillow)
timm >= 0.9.0  # UNI via Hugging Face Hub
conch  # Install: pip install git+https://github.com/Mahmoodlab/CONCH.git
huggingface-hub
```

## Requirements

### System
- CUDA GPU recommended (CPU mode available but slow)
- ~500GB for full HEST breast dataset
- ~100GB for processed outputs

### Credentials
- Hugging Face Hub access for model downloads
- Proper file permissions for input/output directories

---

## Typical Workflow

1. **Download data:**
   ```bash
   python dowload_hest_data.py
   ```

2. **Validate data:**
   ```bash
   python check_hest_alignment.py <slide_id>
   ```

3. **Preprocess and extract features:**
   ```bash
   python data_preprocessing.py --patch_dir ... --raw_st_dir ... --out_root ...
   ```

4. **Analyze gene coverage:**
   ```bash
   python analyze_missinggene.py
   ```

5. **Update gene list:**
   ```bash
   python gene_listupdate.py
   ```

---

## Logging

- **`data_preprocessing.py`** writes comprehensive logs to `<out_root>/prepare_hinge.log`
- **`dowload_hest_data.py`** writes logs to `<output_dir>/download_hest_breast.log`

---

## Notes

- All scripts use absolute paths configured internally—modify paths as needed for your environment
- The pipeline is optimized for breast tissue; organ types can be modified in download script
- Feature extraction uses batch processing for memory efficiency
- Supports cache skipping with `--overwrite` flag in preprocessing
