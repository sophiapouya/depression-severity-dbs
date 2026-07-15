#!/bin/bash
# Submit permutation chunks for one patient + condition as a SLURM array job.
#
# Single mode:
#   sbatch submit_patient.sh DBSTRD001 RVCVS single acc left
#
# Pair mode:
#   sbatch submit_patient.sh DBSTRD001 RVCVS pair acc amy
#
# After all chunks finish, merge with:
#   python dbs_seeg_pca_perm_jobs/merge_pca_dbs_seeg_perm.py \
#       --patient_id <ID> --condition_label <LABEL>

#SBATCH --array=0-19
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=4:00:00
#SBATCH --qos=default_tier
#SBATCH --output=logs/dbs_seeg_pca_%x_%A_%a.out
#SBATCH --error=logs/dbs_seeg_pca_%x_%A_%a.err

export XALT_EXECUTABLE_TRACKING=no
source ~/.bashrc
conda activate neural-analysis

PATIENT_ID=$1
DBS_REGION=$2
MODE=$3       # "single" or "pair"
ARG4=$4       # seeg_region (single) or seeg_left_region (pair)
ARG5=$5       # seeg_hemisphere (single) or seeg_right_region (pair)

N_CHUNKS=20
N_PERMUTATIONS=1000
N_JOBS=8

if [ "$MODE" = "pair" ]; then
    python dbs_seeg_pca_perm_jobs/perm_pca_dbs_seeg_chunk.py \
        --patient_id      "${PATIENT_ID}" \
        --dbs_region      "${DBS_REGION}" \
        --seeg_region     "${ARG4}" \
        --seeg_right_region "${ARG5}" \
        --chunk_id        "${SLURM_ARRAY_TASK_ID}" \
        --n_chunks        "${N_CHUNKS}" \
        --n_permutations  "${N_PERMUTATIONS}" \
        --n_jobs          "${N_JOBS}"
else
    python dbs_seeg_pca_perm_jobs/perm_pca_dbs_seeg_chunk.py \
        --patient_id      "${PATIENT_ID}" \
        --dbs_region      "${DBS_REGION}" \
        --seeg_region     "${ARG4}" \
        --seeg_hemisphere "${ARG5}" \
        --chunk_id        "${SLURM_ARRAY_TASK_ID}" \
        --n_chunks        "${N_CHUNKS}" \
        --n_permutations  "${N_PERMUTATIONS}" \
        --n_jobs          "${N_JOBS}"
fi
