#!/bin/bash
# Submit permutation jobs for all SEEG patients (PCA greymatter, all regions combined).
# Usage: bash seeg_pca_perm_jobs_gm/submit_all.sh
#
# After all jobs finish, merge each patient's results:
#   for P in "${PATIENTS[@]}"; do
#       python seeg_pca_perm_jobs_gm/merge_seeg_pca_perm.py --patient_id "$P" --n_chunks 20
#   done

PATIENTS=(
    "DBSTRD001"
    "DBSTRD002"
    "DBSTRD006"
    "DBSTRD008"
    "DBSTRD010"
    "DBSTRD011"
    "DBSTRD014"
)

# Resolve project root (parent of seeg_pca_perm_jobs_gm/) regardless of where this
# script is called from. --chdir tells SLURM to cd here before running the job.
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
echo "Project root: ${PROJECT_DIR}"

mkdir -p "${PROJECT_DIR}/logs"

for PATIENT in "${PATIENTS[@]}"; do
    echo "Submitting ${PATIENT}..."
    sbatch --job-name="seeg_pca_perm_gm_${PATIENT}" \
        --chdir="${PROJECT_DIR}" \
        "${PROJECT_DIR}/seeg_pca_perm_jobs_gm/submit_patient.sh" "${PATIENT}"
done

echo "All jobs submitted. Monitor with: squeue -u \$USER"
echo ""
echo "After all jobs finish, merge results with:"
for P in "${PATIENTS[@]}"; do
    echo "  python seeg_pca_perm_jobs_gm/merge_seeg_pca_perm.py --patient_id ${P} --n_chunks 20"
done
