#!/bin/bash
# Submit permutation chunks for one SEEG patient + hemisphere as a SLURM array job.
# Usage: sbatch submit_patient.sh DBSTRD001 left
#
# After all chunks finish, merge with:
#   python seeg_pca_perm_jobs_gm_laterality/merge_seeg_pca_gm_laterality_perm.py \
#       --patient_id <ID> --hemisphere <HEMISPHERE>

#SBATCH --array=0-19
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=4:00:00
#SBATCH --qos=default_tier
#SBATCH --output=logs/seeg_pca_gm_lat_%x_%A_%a.out
#SBATCH --error=logs/seeg_pca_gm_lat_%x_%A_%a.err

export XALT_EXECUTABLE_TRACKING=no
source ~/.bashrc
conda activate neural-analysis

PATIENT_ID=$1
HEMISPHERE=$2
N_CHUNKS=20
N_PERMUTATIONS=1000
N_JOBS=8   # match --cpus-per-task

python seeg_pca_perm_jobs_gm_laterality/perm_pca_seeg_gm_laterality_chunk.py \
    --patient_id  "${PATIENT_ID}" \
    --hemisphere  "${HEMISPHERE}" \
    --chunk_id    "${SLURM_ARRAY_TASK_ID}" \
    --n_chunks    "${N_CHUNKS}" \
    --n_permutations "${N_PERMUTATIONS}" \
    --n_jobs      "${N_JOBS}"
