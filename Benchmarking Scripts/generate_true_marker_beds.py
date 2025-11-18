#!/usr/bin/env python3
import os
import pandas as pd
import time
from pathlib import Path

base_dir = "/home/mahe/Downloads/PhD/phd_scripts/stan"

REQUIRED_COLS = {"genome_id", "marker_id", "start", "end"}

def backup_if_exists(path: Path):
    if path.exists():
        ts = time.strftime("%Y%m%dT%H%M%S")
        bak = path.with_name(f"{path.name}.bak.{ts}")
        path.replace(bak)
        print(f"  ↳ Existing {path.name} backed up → {bak.name}")

def safe_write_df_to_bed(df: pd.DataFrame, bed_path: Path):
    """Write df to bed_path safely: temp file then atomic replace. Assumes df has genome_id,start,end columns."""
    tmp = bed_path.with_suffix(bed_path.suffix + ".tmp")
    # ensure order & columns
    out_df = df[["genome_id", "start", "end"]].copy()
    out_df.to_csv(tmp, sep="\t", header=False, index=False)
    # backup original if exists then replace
    if bed_path.exists():
        backup_if_exists(bed_path)
    tmp.replace(bed_path)
    print(f"  ✅ Wrote {bed_path} ({len(out_df)} rows)")

def validate_and_clean_df(df: pd.DataFrame, src_path: Path):
    """
    Validate required columns, coerce numeric starts/ends, remove bad rows,
    remove .fasta suffix from genome_id, drop duplicates, return cleaned df and a report dict.
    """
    report = {"read_rows": len(df), "dropped_non_numeric": 0, "dropped_invalid_coords": 0, "duplicates": 0}

    # Check required columns
    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {src_path}: {', '.join(missing)}")

    working = df.copy()

    # Strip whitespace from string columns (helpful)
    working["genome_id"] = working["genome_id"].astype(str).str.strip()
    working["marker_id"] = working["marker_id"].astype(str).str.strip()

    # Remove .fasta suffix from genome_id if present
    working["genome_id"] = working["genome_id"].str.replace(r"\.fasta$", "", regex=True)

    # Coerce start/end to numeric (integers). If cannot, drop those rows.
    # Use errors='coerce' then dropna
    working["start"] = pd.to_numeric(working["start"], errors="coerce").astype("Int64")
    working["end"] = pd.to_numeric(working["end"], errors="coerce").astype("Int64")
    non_numeric_mask = working["start"].isna() | working["end"].isna()
    if non_numeric_mask.any():
        report["dropped_non_numeric"] = int(non_numeric_mask.sum())
        working = working.loc[~non_numeric_mask].copy()

    # Enforce start >= 1 and start < end
    invalid_mask = (working["start"] < 1) | (working["start"] >= working["end"])
    if invalid_mask.any():
        report["dropped_invalid_coords"] = int(invalid_mask.sum())
        working = working.loc[~invalid_mask].copy()

    # Now convert start to 0-based for BED (remember original is 1-based)
    # but keep end as-is
    working["start"] = (working["start"] - 1).astype(int)
    working["end"] = working["end"].astype(int)

    # Sort and drop exact duplicate genome/start/end triplets (keep first)
    before_dups = len(working)
    working = working.sort_values(["genome_id", "start", "end"])
    working_before = working.copy()
    working = working.drop_duplicates(subset=["genome_id", "start", "end"])
    report["duplicates"] = before_dups - len(working)

    report["written_rows"] = len(working)
    return working, report

def to_bed(df, bed_filename):
    """Wrapper kept for compatibility but uses safe_write_df_to_bed."""
    bed_path = Path(bed_filename)
    safe_write_df_to_bed(df, bed_path)

def process_marker_file(marker_file: str):
    marker_path = Path(marker_file)
    print(f"\n🧬 Processing {marker_path} ...")
    try:
        df = pd.read_csv(marker_path, sep="\t", dtype=str, comment="#")
    except Exception as e:
        print(f"  ❌ Failed to read {marker_file}: {e}")
        return

    if df.empty:
        print(f"  ⚠️ Empty file: {marker_file}")
        return

    try:
        cleaned, report = validate_and_clean_df(df, marker_path)
    except ValueError as e:
        print(f"  ❌ Validation error: {e}")
        return

    print(f"  Rows read: {report['read_rows']}")
    if report["dropped_non_numeric"]:
        print(f"  ⚠ Dropped non-numeric start/end rows: {report['dropped_non_numeric']}")
    if report["dropped_invalid_coords"]:
        print(f"  ⚠ Dropped invalid coordinate rows (start<1 or start>=end): {report['dropped_invalid_coords']}")
    if report["duplicates"]:
        print(f"  ℹ Dropped duplicate rows: {report['duplicates']}")

    # Prepare outputs
    outdir = marker_path.parent / "ground_truth"
    outdir.mkdir(exist_ok=True)

    unique = cleaned[cleaned["marker_id"].str.startswith("unique_")]
    shared = cleaned[cleaned["marker_id"].str.startswith("shared_")]

    unique_bed = outdir / "true_unique_marker_all.bed"
    shared_bed = outdir / "true_shared_marker_all.bed"

    # Write only if there is content; otherwise remove existing target (after backup) or skip
    if not unique.empty:
        safe_write_df_to_bed(unique, unique_bed)
    else:
        print("  ℹ No unique markers found; skipping unique bed write.")

    if not shared.empty:
        safe_write_df_to_bed(shared, shared_bed)
    else:
        print("  ℹ No shared markers found; skipping shared bed write.")

    print(f"  ✅ Finished processing {marker_path.name}: wrote unique={len(unique)} shared={len(shared)}")

# Main traversal
for size_dir in sorted(os.listdir(base_dir)):
    size_path = os.path.join(base_dir, size_dir)
    if not os.path.isdir(size_path):
        continue

    for sim_dir in sorted(os.listdir(size_path)):
        sim_path = os.path.join(size_path, sim_dir)
        if not os.path.isdir(sim_path):
            continue

        marker_file = os.path.join(sim_path, "marker_positions_all.txt")
        if not os.path.exists(marker_file):
            continue

        process_marker_file(marker_file)

