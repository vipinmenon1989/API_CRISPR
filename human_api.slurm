#!/bin/bash
#SBATCH --job-name=Cas_HT_Scout
#SBATCH --output=scout_%j.out
#SBATCH --error=scout_%j.err
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=50      # Optimized for 50-worker parallel scout
#SBATCH --mem=80G               # Sufficient for API/JSON handling
#SBATCH --time=12:00:00         # Realistic wall-time for 10k paper scan

# --- 0. IMMUNITY: Prevent Nested Parallelism at the Shell Level ---
# These must be exported before launching Python
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# --- 1. Environment Setup ---
cd /local/projects-t3/lilab/vmenon/API_Cas12/
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate CRISPR

# --- 2. Execution ---
# Using the scout script to identify 'Big Boys' (Human, N>2000, HT)
python human_api_approach.py --output big_boys_manifest.csv