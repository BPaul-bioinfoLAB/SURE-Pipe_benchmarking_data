#!/usr/bin/env python3
"""
Robust marker evaluation for FUR, SURE-Pipe, and KEC
---------------------------------------------------
✔ Compatible with new folder structure
✔ Auto-discovers marker FASTAs
✔ Parallel BLAST evaluation
✔ Computes TP/FP/TN/FN → Sn, Sp, Acc, MCC
✔ Exports:
   - Global summary (all tools)
   - Mean summary (by size/tool/type)
   - Tool-wise summary files
"""

import os
import math
import numpy as np
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess

# ============================================================
# 🧭 CONFIGURATION
# ============================================================
BASE_DIR = "/home/mahe/Downloads/PhD/phd_scripts/SURE-Pipe_benchmarking_datas"
GLOBAL_SUMMARY = os.path.join(BASE_DIR, "marker_eval_summary_all.tsv")
MEAN_SUMMARY = os.path.join(BASE_DIR, "marker_eval_mean_by_size.tsv")
MAX_THREADS = 4  # parallel BLASTs per simulation

# ============================================================
# 🧰 HELPER FUNCTIONS
# ============================================================
def run_cmd(args):
    """Run a shell command safely with error handling."""
    print(f"🔹 Executing:\n   {' '.join(args)}")
    try:
        subprocess.run(args, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ Command failed: {' '.join(e.cmd)}")
        print(f"Exit code: {e.returncode}")
        raise

def genome_length_sum(genomes):
    total = 0
    for g in genomes:
        for record in SeqIO.parse(g, "fasta"):
            total += len(record.seq)
    return total

def bed_total_length(bed):
    total = 0
    if not os.path.exists(bed):
        return 0
    with open(bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) >= 3:
                total += abs(int(cols[2]) - int(cols[1]))
    return total

def compute_metrics(TP, FP, TN, FN):
    TP, FP, TN, FN = map(float, (TP, FP, TN, FN))
    Sn = TP / (TP + FN) if (TP + FN) else 0
    Sp = TN / (TN + FP) if (TN + FP) else 0
    Acc = (TP + TN) / (TP + FP + FN + TN) if (TP + FP + FN + TN) else 0
    MCC = "NA"
    denom = (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
    if denom > 0:
        MCC = (TP * TN - FP * FN) / math.sqrt(denom)
    return Sn, Sp, Acc, MCC

def blast_to_bed(blast_out, bed_out):
    """Convert BLAST tab output → unique, sorted, merged BED file."""
    intervals_by_chrom = {}

    with open(blast_out) as fin:
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            chrom = parts[1]
            try:
                s1, s2 = int(parts[8]), int(parts[9])
            except (ValueError, IndexError):
                continue
            start, end = sorted([s1, s2])
            if chrom not in intervals_by_chrom:
                intervals_by_chrom[chrom] = []
            intervals_by_chrom[chrom].append((start, end))

    merged_intervals = []
    for chrom, ivals in intervals_by_chrom.items():
        # Sort by start coordinate
        ivals.sort(key=lambda x: x[0])
        merged = []
        for start, end in ivals:
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                # Merge overlapping or adjacent intervals
                merged[-1][1] = max(merged[-1][1], end)
        # Add to final list
        for s, e in merged:
            merged_intervals.append((chrom, s, e))

    # Sort globally for deterministic output
    merged_intervals.sort(key=lambda x: (x[0], x[1], x[2]))

    with open(bed_out, "w") as fout:
        for chrom, start, end in merged_intervals:
            fout.write(f"{chrom}\t{start}\t{end}\n")

def compare_bed_files(true_bed_file, pred_bed_file, fasta_genomes):
    """Compare predicted vs true BED masks."""
    genome_lengths = {}
    for g in fasta_genomes:
        for record in SeqIO.parse(g, "fasta"):
            genome_lengths[record.id] = len(record.seq)

    true_masks = {g: np.zeros(genome_lengths[g], dtype=bool) for g in genome_lengths}
    pred_masks = {g: np.zeros(genome_lengths[g], dtype=bool) for g in genome_lengths}

    with open(true_bed_file) as tb:
        for line in tb:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, s, e = parts[0], int(parts[1]), int(parts[2])
            if chrom not in true_masks:
                continue
            end = min(e, genome_lengths[chrom])
            true_masks[chrom][s:end] = True

    with open(pred_bed_file) as pb:
        for line in pb:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, s, e = parts[0], int(parts[1]), int(parts[2])
            if chrom not in pred_masks:
                continue
            end = min(e, genome_lengths[chrom])
            pred_masks[chrom][s:end] = True

    TP = FP = TN = FN = 0
    for chrom in genome_lengths:
        tmask = true_masks[chrom]
        pmask = pred_masks[chrom]
        TP += np.sum(tmask & pmask)
        FP += np.sum(~tmask & pmask)
        FN += np.sum(tmask & ~pmask)
        TN += np.sum(~tmask & ~pmask)
    return TP, FP, TN, FN

# ============================================================
# 🚀 MAIN EXECUTION
# ============================================================
results = []

for size_dir in sorted(os.listdir(BASE_DIR)):
    size_path = os.path.join(BASE_DIR, size_dir)
    if not os.path.isdir(size_path):
        continue
    print(f"\n🔬 Genome size: {size_dir}")

    for sim_dir in sorted(os.listdir(size_path)):
        sim_path = os.path.join(size_path, sim_dir)
        if not os.path.isdir(sim_path):
            continue
        print(f"   ▶ Simulation: {sim_dir}")

        target_genomes = sorted(str(p) for p in Path(sim_path, "target").glob("*.fasta"))
        non_target_genomes = sorted(str(p) for p in Path(sim_path, "non_target").glob("*.fasta"))
        all_genomes = target_genomes + non_target_genomes
        if not all_genomes:
            print(f"      ⚠️ No genomes found for {sim_dir}")
            continue

        genome_len_total = genome_length_sum(all_genomes)
        print(f"      📏 Total genome length: {genome_len_total:,} bp")

        # === BLAST DB CREATION ===
        db_dir = Path(sim_path, "blast_db")
        db_dir.mkdir(exist_ok=True)
        db_name = db_dir / "all_genomes_db"
        combined_fa = db_dir / "all_genomes.fasta"

        if not Path(str(db_name) + ".nhr").exists():
            print(f"      🧱 Creating BLAST DB ...")
            with open(combined_fa, "w") as out:
                for g in all_genomes:
                    with open(g) as f:
                        out.write(f.read())
            run_cmd(["makeblastdb", "-in", str(combined_fa), "-dbtype", "nucl", "-out", str(db_name)])
        else:
            print(f"      ✅ Using existing BLAST DB")

        # === Ground truth ===
        true_dir = Path(sim_path, "ground_truth")
        true_unique = true_dir / "true_unique_marker_all.bed"
        true_shared = true_dir / "true_shared_marker_all.bed"

        # === Marker sources (auto-discovery) ===
        tool_markers = {
            "FUR": {
                "unique": list(Path(sim_path, "fur_results").glob("*_unique.fasta"))
            },
            "KEC": {
                "unique": list(Path(sim_path, "kec_results").glob("KEC_unique*.fasta"))
            },
            "SURE-Pipe": {
                "unique": list(Path(sim_path, "surepipe_results/shared").rglob("unique_genomic_regions.fasta")),
                "shared": list(Path(sim_path, "surepipe_results/shared").rglob("shared_genomic_regions.fasta"))
            }
        }

        outdir = Path(sim_path, "blast_eval_results")
        outdir.mkdir(exist_ok=True)

        # === Run BLASTs in parallel ===
        futures = []
        with ThreadPoolExecutor(max_workers=MAX_THREADS) as executor:
            for tool, subdict in tool_markers.items():
                for kind, fasta_list in subdict.items():
                    for fa_path in fasta_list:
                        if not fa_path.exists():
                            print(f"      ⚠️ Missing {tool} ({kind}) FASTA: {fa_path}")
                            continue
                        blast_out = outdir / f"{tool}_{kind}_{size_dir}_{sim_dir}_blast.tsv"
                        args = [
                            "blastn",
                            "-query", str(fa_path),
                            "-db", str(db_name),
                            "-outfmt", "6",
                            "-num_threads", "2",
                            "-out", str(blast_out)
                        ]
                        futures.append(executor.submit(run_cmd, args))
            for future in as_completed(futures):
                future.result()

        # === Convert and compare ===
        for tool, subdict in tool_markers.items():
            for kind, fasta_list in subdict.items():
                for fa_path in fasta_list:
                    if not fa_path.exists(): 
                        continue
                    true_bed = true_shared if kind == "shared" else true_unique
                    if not true_bed.exists():
                        continue

                    blast_out = outdir / f"{tool}_{kind}_{size_dir}_{sim_dir}_blast.tsv"
                    bed_out = outdir / f"{tool}_{kind}_{size_dir}_{sim_dir}_predicted.bed"
                    if not blast_out.exists():
                        continue

                    blast_to_bed(blast_out, bed_out)
                    TP, FP, TN, FN = compare_bed_files(true_bed, bed_out, all_genomes)
                    Sn, Sp, Acc, MCC = compute_metrics(TP, FP, TN, FN)

                    results.append([
                        size_dir, sim_dir, tool, kind, genome_len_total,
                        bed_total_length(true_bed), TP, FP, TN, FN, Sn, Sp, Acc, MCC
                    ])
                    print(f"      ✅ {tool} ({kind}) → Sn={Sn:.3f}, Sp={Sp:.3f}, Acc={Acc:.3f}")

# ============================================================
# 📊 SAVE RESULTS
# ============================================================
# Create main DataFrame
df = pd.DataFrame(results, columns=[
    "GenomeSize", "Simulation", "Tool", "MarkerType",
    "TotalGenomeLength", "TrueRegionBP", "TP", "FP", "TN", "FN",
    "Sensitivity", "Specificity", "Accuracy", "MCC"
])

# --- Save raw (unmodified) results ---
RAW_FILE = os.path.join(BASE_DIR, "marker_eval_raw_results.tsv")
df.to_csv(RAW_FILE, sep="\t", index=False)
print(f"💾 Raw results (unmodified) → {RAW_FILE}")

# --- Clean and convert numeric fields ---
df["MCC"] = pd.to_numeric(df["MCC"], errors="coerce")
for col in ["Sensitivity", "Specificity", "Accuracy", "TP", "FP", "TN", "FN"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# --- Save cleaned summary ---
df.to_csv(GLOBAL_SUMMARY, sep="\t", index=False)
print(f"💾 Cleaned summary → {GLOBAL_SUMMARY}")

# --- TOOL-WISE RESULTS (no averaging) ---
toolwise_file = os.path.join(BASE_DIR, "marker_eval_by_tool.tsv")
df_toolwise = df.sort_values(["Tool", "MarkerType", "GenomeSize", "Simulation"])
df_toolwise.to_csv(toolwise_file, sep="\t", index=False)
print(f"🧩 Tool-wise detailed results (no averaging) → {toolwise_file}")

# --- Aggregate mean and sums ---
agg = {
    "TP": "sum", "FP": "sum", "TN": "sum", "FN": "sum",
    "Sensitivity": "mean", "Specificity": "mean",
    "Accuracy": "mean", "MCC": "mean"
}
mean_df = (
    df.groupby(["GenomeSize", "Tool", "MarkerType"], as_index=False)
      .agg(agg)
)

# --- Save mean summary ---
mean_df.to_csv(MEAN_SUMMARY, sep="\t", index=False)
print(f"📊 Mean summary (with TP/FP/TN/FN) → {MEAN_SUMMARY}")

# --- Optional print summary ---
print("\n🧠 Mean results by genome size:")
print(mean_df[["GenomeSize", "Tool", "MarkerType", "Sensitivity", "Specificity", "Accuracy", "MCC"]]
      .round(3)
      .to_string(index=False))


