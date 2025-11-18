#!/bin/bash
# ---------------------------------------------------------------------
# 🧬 Benchmark script for FUR, KEC, and SURE-Pipe across all genome sizes & simulations
# ---------------------------------------------------------------------

BASE_DIR="$HOME/Downloads/PhD/phd_scripts/stan"

echo "🔍 Starting FUR + KEC + SURE-Pipe benchmarking in $BASE_DIR"

for size_dir in "$BASE_DIR"/*; do
    [ -d "$size_dir" ] || continue
    size_label=$(basename "$size_dir")
    echo "📁 Genome size: $size_label"

    for sim_dir in "$size_dir"/sim_*; do
        [ -d "$sim_dir" ] || continue
        sim_label=$(basename "$sim_dir")
        echo "   🔬 Simulation: $sim_label"

        target_dir="$sim_dir/target"
        non_target_dir="$sim_dir/non_target"
        ref_genome="$target_dir/target_genome_1.fasta"

        if [ ! -f "$ref_genome" ]; then
            echo "      ⚠️ Skipping $sim_label — no reference genome found"
            continue
        fi
        if [ ! -d "$target_dir" ] || [ ! -d "$non_target_dir" ]; then
            echo "      ⚠️ Skipping $sim_label — missing target/ or non_target/"
            continue
        fi

        cd "$sim_dir" || continue

        # Create result directories
        mkdir -p fur_results kec_results surepipe_results/unique surepipe_results/shared

        # -------------------------------
        # 🧱 FUR Step 1: makeFurDb
        # -------------------------------
        echo "      🧱 [FUR-DB] Building database..."
        /usr/bin/time -v makeFurDb \
            -t "$target_dir" \
            -n "$non_target_dir" \
            -r "target_genome_1.fasta" \
            -d "fur_results/${size_label}_${sim_label}_fur_db" \
            -o \
            > "fur_results/stdout_fur_db_${size_label}_${sim_label}.log" \
            2> "fur_results/runtime_mem_fur_db_${size_label}_${sim_label}.log"

        # -------------------------------
        # ⚙️ FUR Step 2: fur run
        # -------------------------------
        echo "      ⚙️ [FUR-RUN] Running FUR..."
        /usr/bin/time -v fur \
            -d "fur_results/${size_label}_${sim_label}_fur_db" \
            -w "20" \
            > "fur_results/fur_run_${size_label}_${sim_label}_unique_raw.fasta" \
            2> "fur_results/runtime_mem_fur_run_${size_label}_${sim_label}.log"

        # 🧹 Clean FUR FASTA output (multi-record safe)
        echo "      🧹 [CLEAN] Formatting FUR FASTA output..."
        awk '/^>/ { if (NR>1) print ""; print $1; next } { printf "%s", $0 } END { print "" }' \
            "fur_results/fur_run_${size_label}_${sim_label}_unique_raw.fasta" \
            > "fur_results/fur_run_${size_label}_${sim_label}_unique.fasta"

        # -------------------------------
        # 🧬 KEC Step
        # -------------------------------
        echo "      🧬 [KEC] Running KEC unique region extraction..."
        cat "$target_dir"/*.fasta > "kec_results/master_target_${size_label}_${sim_label}.fasta"

        /usr/bin/time -v kec exclude \
            -t "kec_results/master_target_${size_label}_${sim_label}.fasta" \
            -n "$non_target_dir" \
            -o "kec_results/KEC_unique_${size_label}_${sim_label}.fasta" \
            -r \
            -min 100 \
            > "kec_results/stdout_kec_${size_label}_${sim_label}.log" \
            2> "kec_results/runtime_mem_kec_${size_label}_${sim_label}.log"

        # -------------------------------
        # 🔗 SURE-Pipe Step 2: Shared regions
        # -------------------------------
        echo "      🔗 [SURE-PIPE] Running (with shared regions)..."
        /usr/bin/time -v SURE-Pipe run \
            --reference_genome "$ref_genome" \
            --target_dir "$target_dir" \
            --neighbour_dir "$non_target_dir" \
            --output_dir "surepipe_results/shared/results_shared_${size_label}_${sim_label}" \
            --shared_region yes \
            --mode group \
            > "surepipe_results/shared/stdout_sure_pipe_${size_label}_${sim_label}_shared.log" \
            2> "surepipe_results/shared/runtime_mem_sure_pipe_${size_label}_${sim_label}_shared.log"

        echo "      ✅ Completed $sim_label"
        echo "---------------------------------------------------"

        cd "$BASE_DIR" || exit
    done
done

echo "🎯 All FUR + KEC + SURE-Pipe runs completed and organized successfully!"

