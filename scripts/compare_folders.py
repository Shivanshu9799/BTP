"""
Compare two folders and report which files are missing from each side.

Usage:
    python compare_folders.py /path/to/folder1 /path/to/folder2

Typical use case here:
    folder1 = expected files (e.g. extracted from GDC manifest.txt or metadata.json)
    folder2 = actually downloaded WSI folder (gdc-client output dir)

Handles two common GDC download layouts:
  1. Flat folder:      folder/<file_id>_<file_name>.svs
  2. gdc-client layout: folder/<file_id>/<file_name>.svs   (one subfolder per file_id)

Matching is done by FILE NAME (not file_id), since that's what you'd
usually eyeball. If you want to match by file_id instead, use --by-id.
"""

import os
import sys
import argparse
import json


def collect_filenames(folder, by_id=False):
    """
    Walk a folder and return a dict {filename: full_path}.
    Handles gdc-client's <file_id>/<filename> subfolder structure automatically.
    """
    result = {}
    if not os.path.isdir(folder):
        print(f"ERROR: '{folder}' is not a valid directory.")
        sys.exit(1)

    for root, dirs, files in os.walk(folder):
        for f in files:
            # skip hidden/system/log files gdc-client sometimes drops
            if f.startswith('.') or f.endswith('.log') or f.endswith('.partial'):
                continue
            full_path = os.path.join(root, f)
            key = f if not by_id else os.path.basename(root)  # subfolder name = file_id in gdc-client layout
            if key in result:
                # duplicate name across subfolders — keep track of both
                if isinstance(result[key], list):
                    result[key].append(full_path)
                else:
                    result[key] = [result[key], full_path]
            else:
                result[key] = full_path
    return result


def collect_from_manifest(manifest_path):
    """
    If folder1 is actually a GDC manifest.txt (tab-separated: id, filename, md5, size, state),
    parse expected filenames from it instead of walking a directory.
    """
    expected = {}
    with open(manifest_path) as f:
        header = f.readline().strip().split('\t')
        try:
            fname_idx = header.index('filename')
        except ValueError:
            fname_idx = 1  # fallback: GDC manifest column order is id, filename, md5, size, state
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > fname_idx:
                expected[parts[fname_idx]] = line.strip()
    return expected


def collect_from_metadata_json(json_path):
    """
    If folder1 is actually a GDC metadata.repository.json export,
    parse expected filenames from it instead of walking a directory.
    """
    with open(json_path) as f:
        data = json.load(f)
    expected = {}
    for entry in data:
        fname = entry.get('file_name')
        if fname:
            expected[fname] = entry.get('file_id', '')
    return expected


def load_side(path, by_id=False):
    """Decide whether 'path' is a folder, a manifest.txt, or a metadata.json — and load accordingly."""
    if os.path.isdir(path):
        return collect_filenames(path, by_id=by_id)
    elif path.endswith('.json'):
        return collect_from_metadata_json(path)
    elif path.endswith('.txt') or path.endswith('.tsv'):
        return collect_from_manifest(path)
    else:
        print(f"ERROR: Don't know how to read '{path}' (expected folder, .json, or .txt manifest).")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Compare two folders/manifests and report missing files.")
    parser.add_argument("source", help="Expected file list: folder, manifest.txt, or metadata.json")
    parser.add_argument("downloaded", help="Actual downloaded folder to check against")
    parser.add_argument("--by-id", action="store_true", help="Match by file_id (subfolder name) instead of filename")
    parser.add_argument("--out", default="missing_report.txt", help="Output report file")
    args = parser.parse_args()

    print(f"Reading expected file list from: {args.source}")
    expected = load_side(args.source, by_id=args.by_id)
    print(f"  -> {len(expected)} expected entries")

    print(f"Reading downloaded folder: {args.downloaded}")
    downloaded = collect_filenames(args.downloaded, by_id=args.by_id)
    print(f"  -> {len(downloaded)} files found on disk")

    expected_keys = set(expected.keys())
    downloaded_keys = set(downloaded.keys())

    missing_from_downloaded = sorted(expected_keys - downloaded_keys)   # in source, not on disk
    extra_in_downloaded = sorted(downloaded_keys - expected_keys)       # on disk, not in source

    print(f"\n{'='*60}")
    print(f"MISSING from downloaded folder : {len(missing_from_downloaded)}")
    print(f"EXTRA in downloaded (unexpected): {len(extra_in_downloaded)}")
    print(f"MATCHED                         : {len(expected_keys & downloaded_keys)}")
    print(f"{'='*60}\n")

    with open(args.out, 'w') as f:
        f.write(f"Source: {args.source}\n")
        f.write(f"Downloaded folder: {args.downloaded}\n")
        f.write(f"Expected: {len(expected)} | On disk: {len(downloaded)}\n\n")

        f.write(f"--- MISSING FROM DOWNLOADED FOLDER ({len(missing_from_downloaded)}) ---\n")
        for name in missing_from_downloaded:
            f.write(f"{name}\n")

        f.write(f"\n--- EXTRA FILES IN DOWNLOADED FOLDER (not in source list) ({len(extra_in_downloaded)}) ---\n")
        for name in extra_in_downloaded:
            f.write(f"{name}\n")

    print(f"Full report written to: {args.out}")

    if missing_from_downloaded:
        print("\nFirst 10 missing files:")
        for name in missing_from_downloaded[:10]:
            print(f"  - {name}")
        if len(missing_from_downloaded) > 10:
            print(f"  ... and {len(missing_from_downloaded) - 10} more (see {args.out})")


if __name__ == "__main__":
    main()
