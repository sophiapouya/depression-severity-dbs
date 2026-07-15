#!/bin/bash
# Submit SEEG PCA greymatter laterality permutation jobs for all patients × hemispheres.
# Hemispheres: all, left, right
# Usage: bash seeg_pca_perm_jobs_gm_laterality/submit_all.sh
#
# After all jobs finish, merge with:
#   for P in "${PATIENTS[@]}"; do
#     for H in "${HEMISPHERES[@]}"; do
#       python seeg_pca_perm_jobs_gm_laterality/merge_seeg_pca_gm_laterality_perm.py \
#           --patient_id "$P" --hemisphere "$H" --n_chunks 20
#     done
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

HEMISPHERES=("all" "left" "right")

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
echo "Project root: ${PROJECT_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

for PATIENT in "${PATIENTS[@]}"; do
    for HEMI in "${HEMISPHERES[@]}"; do
        echo "Submitting ${PATIENT} hemisphere=${HEMI}..."
        sbatch --job-name="seeg_gm_lat_${PATIENT}_${HEMI}" \
            --chdir="${PROJECT_DIR}" \
            "${PROJECT_DIR}/seeg_pca_perm_jobs_gm_laterality/submit_patient.sh" \
            "${PATIENT}" "${HEMI}"
    done
done

echo ""
echo "All jobs submitted ($(( ${#PATIENTS[@]} * ${#HEMISPHERES[@]} )) total)."
echo "Monitor with: squeue -u \$USER"
echo ""
echo "After all jobs finish, merge with:"
for P in "${PATIENTS[@]}"; do
    for H in "${HEMISPHERES[@]}"; do
        echo "  python seeg_pca_perm_jobs_gm_laterality/merge_seeg_pca_gm_laterality_perm.py --patient_id ${P} --hemisphere ${H} --n_chunks 20"
    done
done
