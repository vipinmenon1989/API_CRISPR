#!/bin/bash
#SBATCH --job-name=Cas13_VIP_Recovery
#SBATCH --output=recovery_%j.out
#SBATCH --error=recovery_%j.err
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --cpus-per-task=2           # Only 2 CPUs needed for single-threaded recovery
#SBATCH --mem=100G                   # Low memory is fine
#SBATCH --time=22:00:00             # 2 hours is plenty for a few VIPs

# --- Environment Setup ---
cd /local/projects-t3/lilab/vmenon/API_Cas12/
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate CRISPR

# --- Execution ---
echo "Starting VIP Recovery at $(date)"
python recovery_cas13.py
echo "Recovery finished at $(date)"
