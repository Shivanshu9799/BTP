import h5py

# Replace with your actual path
h5_path = "patches/SPA71.h5"

with h5py.File(h5_path, 'r') as f:
    print("Keys in file:", list(f.keys()))
    for key in f.keys():
        print(f"  {key}: shape={f[key].shape}, dtype={f[key].dtype}")
