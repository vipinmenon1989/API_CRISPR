#!/bin/bash
#SBATCH --job-name=Cas12_Parallel_Miner
#SBATCH --output=logs_%j.out
#SBATCH --error=logs_%j.err
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G
#SBATCH --time=248:00:00

# --- ENVIRONMENT SETUP ---
cd /local/projects-t3/lilab/vmenon/API_Cas12/

# Load Anaconda
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate CRISPR

# --- CRITICAL: THREAD LOCKING ---
# We must force Python libraries (NumPy, PyTorch) to be single-threaded.
# If we don't, 50 workers x 50 threads = 2500 threads = NODE CRASH.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# --- EXECUTION ---
echo "------------------------------------------------------"
echo "🚀 Starting Parallel Cas12 Mining Job at $(date)"
echo "   - Allocated CPUs: $SLURM_CPUS_PER_TASK"
echo "   - Node: $(hostname)"
echo "------------------------------------------------------"

# Run the script
python dataloader_cas13.py

echo "------------------------------------------------------"
echo "✅ Job finished at $(date)"
echo "------------------------------------------------------"
