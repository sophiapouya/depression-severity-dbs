#!/bin/bash
# Submit permutation chunks for one DBS patient + probe type as a SLURM array job.
# Usage: sbatch submit_patient.sh <patient_id> <probe_type>
# Example: sbatch submit_patient.sh DBSTRD001 ALL

#SBATCH --array=0-19
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=4:00:00
#SBATCH --qos=default_tier
#SBATCH --output=logs/dbs_pca_perm_%x_%A_%a.out
#SBATCH --error=logs/dbs_pca_perm_%x_%A_%a.err

export XALT_EXECUTABLE_TRACKING=no
source ~/.bashrc
conda activate neural-analysis

PATIENT_ID=$1
PROBE_TYPE=$2
N_CHUNKS=20
N_PERMUTATIONS=1000
N_JOBS=8   # match --cpus-per-task

python dbs_pca_perm_jobs/perm_pca_dbs_chunk.py \
    --patient_id  "${PATIENT_ID}" \
    --probe_type  "${PROBE_TYPE}" \
    --chunk_id    "${SLURM_ARRAY_TASK_ID}" \
    --n_chunks    "${N_CHUNKS}" \
    --n_permutations "${N_PERMUTATIONS}" \
    --n_jobs      "${N_JOBS}"
