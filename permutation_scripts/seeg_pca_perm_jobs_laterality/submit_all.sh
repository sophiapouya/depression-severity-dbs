#!/bin/bash
# Submit SEEG PCA laterality permutation jobs for all patients × all conditions.
#
# Conditions (region × hemisphere):
#   acc/all, acc/left, acc/right
#   amy/all, amy/left, amy/right
#   dlpfc/all, dlpfc/left, dlpfc/right
#   vmpfc/all, vmpfc/left, vmpfc/right
#   ofc/all, ofc/left, ofc/right
#   all/left, all/right
#
# Usage: bash seeg_pca_perm_jobs_laterality/submit_all.sh

PATIENTS=(
    "DBSTRD001"
    "DBSTRD002"
    "DBSTRD006"
    "DBSTRD008"
    "DBSTRD010"
    "DBSTRD011"
    "DBSTRD014"
)

# (region, hemisphere) pairs
CONDITIONS=(
    "acc all"
    "acc left"
    "acc right"
    "amy all"
    "amy left"
    "amy right"
    "dlpfc all"
    "dlpfc left"
    "dlpfc right"
    "vmpfc all"
    "vmpfc left"
    "vmpfc right"
    "ofc all"
    "ofc left"
    "ofc right"
    "all left"
    "all right"
)

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
echo "Project root: ${PROJECT_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

for PATIENT in "${PATIENTS[@]}"; do
    for COND in "${CONDITIONS[@]}"; do
        REGION=$(echo "$COND" | awk '{print $1}')
        HEMI=$(echo "$COND" | awk '{print $2}')
        JOB_NAME="seeg_lat_${PATIENT}_${REGION}_${HEMI}"
        echo "Submitting ${PATIENT} region=${REGION} hemisphere=${HEMI}..."
        sbatch --job-name="${JOB_NAME}" \
            --chdir="${PROJECT_DIR}" \
            "${PROJECT_DIR}/seeg_pca_perm_jobs_laterality/submit_patient.sh" \
            "${PATIENT}" "${REGION}" "${HEMI}"
    done
done

echo ""
echo "All jobs submitted ($(( ${#PATIENTS[@]} * ${#CONDITIONS[@]} )) total)."
echo "Monitor with: squeue -u \$USER"
echo ""
echo "After all jobs finish, merge with:"
for PATIENT in "${PATIENTS[@]}"; do
    for COND in "${CONDITIONS[@]}"; do
        REGION=$(echo "$COND" | awk '{print $1}')
        HEMI=$(echo "$COND" | awk '{print $2}')
        echo "  python seeg_pca_perm_jobs_laterality/merge_seeg_pca_laterality_perm.py --patient_id ${PATIENT} --region ${REGION} --hemisphere ${HEMI} --n_chunks 20"
    done
done
