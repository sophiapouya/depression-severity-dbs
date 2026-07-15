#!/bin/bash
# Submit DBS PCA permutation jobs for all patients × all probe types.
# Usage: bash dbs_pca_perm_jobs/submit_all.sh

PATIENTS=(
    "DBSTRD001"
    "DBSTRD002"
    "DBSTRD006"
    "DBSTRD008"
    "DBSTRD010"
    "DBSTRD011"
    "DBSTRD014"
)

PROBE_TYPES=(
    "ALL"
    "LEFT"
    "RIGHT"
    "LSCC"
    "RSCC"
    "RVCVS"
    "LVCVS"
    "SCC"
    "VCVS"
)

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
echo "Project root: ${PROJECT_DIR}"

mkdir -p "${PROJECT_DIR}/logs"

for PATIENT in "${PATIENTS[@]}"; do
    for PROBE in "${PROBE_TYPES[@]}"; do
        echo "Submitting ${PATIENT} / ${PROBE}..."
        sbatch --job-name="dbs_pca_${PATIENT}_${PROBE}" \
            --chdir="${PROJECT_DIR}" \
            "${PROJECT_DIR}/dbs_pca_perm_jobs/submit_patient.sh" "${PATIENT}" "${PROBE}"
    done
done

echo ""
echo "All jobs submitted ($(( ${#PATIENTS[@]} * ${#PROBE_TYPES[@]} )) total)."
echo "Monitor with: squeue -u \$USER"
