#!/bin/bash
# Submit permutation chunks for one SEEG patient (all regions) as a SLURM array job (PCA greymatter version).
# Usage: sbatch submit_patient.sh DBSTRD001
#
# Array index = chunk_id (0 to N_CHUNKS-1).
# Chunk 0 also runs the true decoder.
# After all chunks finish, run: python seeg_pca_perm_jobs_gm/merge_seeg_pca_perm.py --patient_id <ID>

#SBATCH --array=0-19
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=4:00:00
#SBATCH --qos=default_tier
#SBATCH --output=logs/seeg_pca_perm_gm_%x_%A_%a.out
#SBATCH --error=logs/seeg_pca_perm_gm_%x_%A_%a.err

export XALT_EXECUTABLE_TRACKING=no
source ~/.bashrc
conda activate neural-analysis

PATIENT_ID=$1
N_CHUNKS=20
N_PERMUTATIONS=1000
N_JOBS=8   # match --cpus-per-task

python seeg_pca_perm_jobs_gm/perm_pca_seeg_chunk.py \
    --patient_id "${PATIENT_ID}" \
    --chunk_id "${SLURM_ARRAY_TASK_ID}" \
    --n_chunks "${N_CHUNKS}" \
    --n_permutations "${N_PERMUTATIONS}" \
    --n_jobs "${N_JOBS}"
