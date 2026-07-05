import os
import logging
import pandas as pd
#from hest import download_hest  # official helper
from huggingface_hub import snapshot_download, hf_hub_download
OUT_DIR = "/mnt/wwn-0x5000c500e64c5cc0/shivanhu/HEST_BREAST"
os.makedirs(OUT_DIR, exist_ok=True)

LOG_FILE = os.path.join(OUT_DIR, "download_hest_breast.log")
logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

# Step 1: Load metadata
logging.info("Loading HEST metadata...")
meta_df = pd.read_csv("hf://datasets/MahmoodLab/hest/HEST_v1_3_0.csv")
logging.info(f"Total samples: {len(meta_df)}")

# Step 2: Filter breast
breast_df = meta_df[
    (meta_df["organ"] == "Breast") ].copy()
logging.info(f"Breast samples found: {len(breast_df)}")

# Step 3: Save filtered metadata
breast_df.to_csv(
    os.path.join(OUT_DIR, "HEST_v1_3_0_BREAST_METADATA.csv"),
    index=False,
)

# Step 4: Build official patterns
ids_to_query = breast_df["id"].values
list_patterns = [f"*{sid}[_.]**" for sid in ids_to_query]
logging.info(f"Total patterns: {len(list_patterns)}")
logging.info(f"Sample patterns: {list_patterns[:3]}")

# Step 5: Download
try:
    snapshot_download(
        repo_id="MahmoodLab/hest",
        repo_type="dataset",
        local_dir=OUT_DIR,
        allow_patterns=list_patterns,
        resume_download=True,
    )
    logging.info("Download completed successfully.")
except Exception as e:
    logging.error(f"Download failed: {e}")
    raise
