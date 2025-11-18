#!/bin/bash
# ============================================================
# 🧬 Parallel 4 Mb genome simulations (150 runs)
# ============================================================
# For each simulation:
#   1. Runs simulate_genomes_with_stan.py
#   2. Creates "ANI" folder
#   3. Copies FASTAs (target + non_target)
#   4. Runs ANIclustermap
#   5. Removes FASTA files in ANI/
# ============================================================

# ---- Config ----
SCRIPT="simulate_genomes_with_stan.py"
OUT_BASE="$HOME/Downloads/PhD/phd_scripts/stan/4mb"
STAN_BIN="/usr/local/bin/stan"
GENOME_SIZE=4000000
TARGETS=5
NONTARGETS=5
TOTAL=150
SEED_BASE=1000
THREADS=4   # how many simulations to run at once (tune for CPU cores)

# ---- Prepare output directory ----
mkdir -p "$OUT_BASE"

# ---- Define function for one simulation ----
run_simulation() {
    local i=$1
    local SIM_DIR="$OUT_BASE/sim_$(printf "%03d" $i)"
    local SEED=$((SEED_BASE + i))

    echo "🚀 [Sim $i] Starting..."
    mkdir -p "$SIM_DIR"

    python3 "$SCRIPT" \
        --size $GENOME_SIZE \
        --targets $TARGETS \
        --nontargets $NONTARGETS \
        --stan_bin "$STAN_BIN" \
        --outdir "$SIM_DIR" \
        --seed $SEED

    echo "✅ [Sim $i] Simulation complete. Preparing ANI step..."

    # ---- Create ANI folder and copy FASTAs ----
    mkdir -p "$SIM_DIR/ANI"
    cp "$SIM_DIR"/target/*.fasta "$SIM_DIR"/non_target/*.fasta "$SIM_DIR/ANI"/

    # ---- Run ANIclustermap ----
    cd "$SIM_DIR/ANI"
    ANIclustermap -i . -o results --annotation --annotation_fmt .1f --cmap_ranges 90,95,97,100

    # ---- Clean up FASTA files ----
    rm -f "$SIM_DIR/ANI"/*.fasta
    cd - >/dev/null

    echo "🧩 [Sim $i] ANI complete and FASTAs cleaned."
    echo "---------------------------------------------"
}

export -f run_simulation
export SCRIPT OUT_BASE STAN_BIN GENOME_SIZE TARGETS NONTARGETS TOTAL SEED_BASE

# ---- Run simulations in parallel ----
seq 1 $TOTAL | parallel -j $THREADS run_simulation {}

echo "🎉 All $TOTAL simulations completed successfully!"
echo "📂 Results available in: $OUT_BASE"

