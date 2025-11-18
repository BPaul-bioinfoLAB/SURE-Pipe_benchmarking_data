#!/usr/bin/env python3
"""
Stan + Marker Genome Simulator (Unified & Enhanced)
---------------------------------------------------
🧬 Uses Stan to generate genomes (target + non-target)
💥 Adds unique + shared markers (random orientation, optional duplication)
📄 Saves:
  - marker_positions_all.txt
  - simulation_info.txt
  - <size_label>_tree.nwk
"""

import os
import random
import argparse
import subprocess
import shutil
from datetime import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# -----------------------------
# Helper functions
# -----------------------------
def random_dna(length):
    """Generate a random DNA sequence of given length."""
    return ''.join(random.choices("ATGC", k=length))

def introduce_mutations(seq, mutation_rate):
    """Introduce random point mutations in a DNA sequence."""
    seq = list(seq)
    for i in range(len(seq)):
        if random.random() < mutation_rate:
            seq[i] = random.choice([b for b in "ATGC" if b != seq[i]])
    return ''.join(seq)

def random_orientation(seq):
    """Randomly flip sequence orientation."""
    if random.random() < 0.5:
        return seq, "+"
    else:
        return str(Seq(seq).reverse_complement()), "-"

def insert_marker(genome_list, marker_seq, size, used_positions, allow_duplication=False):
    """Insert marker randomly into genome sequence."""
    marker_len = len(marker_seq)
    inserted_positions = []
    copies = 2 if allow_duplication and random.random() < 0.25 else 1

    for _ in range(copies):
        for _ in range(10000):
            pos = random.randint(1000, size - marker_len - 1000)
            if not any(abs(pos - p) < marker_len + 500 for p in used_positions):
                genome_list[pos:pos + marker_len] = marker_seq
                used_positions.append(pos)
                inserted_positions.append((pos, pos + marker_len))
                break
    return inserted_positions

def rename_fasta_files(folder, prefix):
    """Rename FASTA files in a folder to prefix_#.fasta."""
    if not os.path.exists(folder):
        print(f"⚠️ Missing folder: {folder}")
        return
    files = sorted([f for f in os.listdir(folder) if f.endswith((".fa", ".fasta"))])
    for i, fname in enumerate(files, 1):
        new_name = f"{prefix}_{i}.fasta"
        shutil.move(os.path.join(folder, fname), os.path.join(folder, new_name))
    print(f"📄 Renamed {len(files)} files → {prefix}_#.fasta in {folder}")

# -----------------------------
# Step 1: Run Stan
# -----------------------------
def run_stan(stan_bin, size, targets, nontargets, outdir, mutation_rate=0.01, marker_mut_rate=0, seed=42):
    """Run Stan simulator and return target/non-target directories + tree path."""
    os.makedirs(outdir, exist_ok=True)

    if size >= 1_000_000:
        label = f"{int(size / 1_000_000)}mb"
    elif size >= 1_000:
        label = f"{int(size / 1_000)}kb"
    else:
        label = f"{size}bp"

    tree_path = os.path.join(outdir, f"{label}_tree.nwk")

    cmd = [
        stan_bin, "-t", str(targets),
        "-n", str(nontargets),
        "-l", str(size),
        "-m", str(mutation_rate),
        "-M", str(marker_mut_rate),
        "-o", "-s", str(seed), "-c"
    ]

    print(f"🚀 Running Stan simulation:\n   {' '.join(cmd)}")
    with open(tree_path, "w") as tree_file:
        subprocess.run(cmd, cwd=outdir, check=True, stdout=tree_file)
    print(f"✅ Stan completed → {tree_path}")

    target_dir = os.path.join(outdir, "targets")
    non_target_dir = os.path.join(outdir, "neighbors")

    rename_fasta_files(target_dir, "target_genome")
    rename_fasta_files(non_target_dir, "non_target_genome")

    return target_dir, non_target_dir, tree_path

# -----------------------------
# Step 2: Marker Insertion
# -----------------------------
def add_markers(size, target_dir, non_target_dir, outdir, seed):
    print(f"\n🧬 [Seed={seed}] Adding markers → {outdir}")
    random.seed(seed)

    os.makedirs(f"{outdir}/target", exist_ok=True)
    os.makedirs(f"{outdir}/non_target", exist_ok=True)

    marker_sizes = [200, 500, 1000, 2000, 5000]
    unique_markers = {f"unique_marker{i+1}": random_dna(marker_sizes[i % len(marker_sizes)]) for i in range(5)}
    shared_markers = {f"shared_marker{i+1}": random_dna(marker_sizes[i % len(marker_sizes)]) for i in range(5)}

    all_marker_positions = []
    print(f"⚙️ Inserting 5 UNIQUE + 5 SHARED markers per genome")

    # ---- Targets ----
    for fname in sorted(os.listdir(target_dir)):
        if not fname.endswith(".fasta"): continue
        seq_record = SeqIO.read(os.path.join(target_dir, fname), "fasta")
        genome_list = list(str(seq_record.seq))
        used_positions = []

        # Unique markers
        for marker_id, seq in unique_markers.items():
            seq_oriented, strand = random_orientation(seq)
            placed = insert_marker(genome_list, seq_oriented, size, used_positions, True)
            for start, end in placed:
                all_marker_positions.append([fname, "target", marker_id, start, end, strand, end - start, len(placed) > 1])

        # Shared markers
        for marker_id, seq in shared_markers.items():
            seq_oriented, strand = random_orientation(seq)
            placed = insert_marker(genome_list, seq_oriented, size, used_positions, True)
            for start, end in placed:
                all_marker_positions.append([fname, "target", marker_id, start, end, strand, end - start, len(placed) > 1])

        SeqIO.write(
            SeqRecord(Seq(''.join(genome_list)), id=fname.split('.')[0], description="target genome with markers"),
            f"{outdir}/target/{fname}", "fasta"
        )

    # ---- Non-targets ----
    for fname in sorted(os.listdir(non_target_dir)):
        if not fname.endswith(".fasta"): continue
        seq_record = SeqIO.read(os.path.join(non_target_dir, fname), "fasta")
        genome_list = list(str(seq_record.seq))
        used_positions = []

        for marker_id, seq in shared_markers.items():
            mutated_seq = introduce_mutations(seq, 0.07)
            seq_oriented, strand = random_orientation(mutated_seq)
            placed = insert_marker(genome_list, seq_oriented, size, used_positions, True)
            for start, end in placed:
                all_marker_positions.append([fname, "non_target", marker_id, start, end, strand, end - start, len(placed) > 1])

        SeqIO.write(
            SeqRecord(Seq(''.join(genome_list)), id=fname.split('.')[0], description="non-target genome with markers"),
            f"{outdir}/non_target/{fname}", "fasta"
        )

    # ---- Save marker info ----
    out_table = f"{outdir}/marker_positions_all.txt"
    with open(out_table, "w") as f:
        f.write("genome_id\tgroup\tmarker_id\tstart\tend\torientation\tlength\tduplication\n")
        for row in all_marker_positions:
            f.write("\t".join(map(str, row)) + "\n")

    print(f"✅ Marker insertion complete → {out_table}")

# -----------------------------
# Main
# -----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="🧬 Stan + Marker Genome Simulator\n"
                    "Simulates target and non-target genomes using Stan, adds markers, "
                    "and saves outputs with simulation metadata.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "📘 Stan installation help:\n"
            "  • Repository: https://github.com/EvolBioInf/stan\n"
            "  • Clone and build:\n"
            "        git clone https://github.com/EvolBioInf/stan.git\n"
            "        cd stan\n"
            "        bash scripts/setup.sh\n"
            "        make\n"
            "  • Binary will be created in the 'bin/' directory.\n"
            "  • You can specify its path using --stan_bin.\n"
        )
    )

    parser.add_argument("--size", type=int, required=True, help="Genome size (bp)")
    parser.add_argument("--targets", type=int, default=5, help="Number of target genomes")
    parser.add_argument("--nontargets", type=int, default=5, help="Number of non-target genomes")
    parser.add_argument("--stan_bin", type=str, default="/usr/local/bin/stan", help="Path to Stan executable")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory")
    parser.add_argument("--seed", type=int, default=None, help="Random seed (optional)")
    args = parser.parse_args()

    # --- ✅ Check if Stan is installed / accessible
    stan_path = shutil.which(args.stan_bin) if os.path.sep not in args.stan_bin else args.stan_bin
    if not stan_path or not os.path.isfile(stan_path):
        print("\n❌ ERROR: Stan executable not found.")
        print("🔍 Tried path:", args.stan_bin)
        print("\n🧩 To install Stan, follow these steps:")
        print("    git clone https://github.com/EvolBioInf/stan.git")
        print("    cd stan")
        print("    bash scripts/setup.sh")
        print("    make")
        print("\nAfter building, you can point to it using:")
        print("    --stan_bin /path/to/stan/bin/stan")
        sys.exit(1)

    # --- Continue simulation if Stan is found
    seed = args.seed if args.seed else int(datetime.now().timestamp()) % 1_000_000
    stan_out = os.path.join(args.outdir, "stan_raw")

    # Run simulation
    target_dir, non_target_dir, tree_path = run_stan(
        stan_path, args.size, args.targets, args.nontargets, stan_out, seed=seed
    )
    add_markers(args.size, target_dir, non_target_dir, args.outdir, seed)

    # Save metadata
    with open(os.path.join(args.outdir, "simulation_info.txt"), "w") as f:
        f.write(f"Genome size: {args.size}\n")
        f.write(f"Targets: {args.targets}\n")
        f.write(f"Non-targets: {args.nontargets}\n")
        f.write(f"Random seed: {seed}\n")
        f.write(f"Coalescent tree: {tree_path}\n")
        f.write(f"Generated on: {datetime.now().isoformat()}\n")

    print(f"\n🎉 Simulation complete!")
    print(f"🧬 Coalescent tree: {tree_path}")
    print(f"📂 Output directory: {args.outdir}")

if __name__ == "__main__":
    main()


