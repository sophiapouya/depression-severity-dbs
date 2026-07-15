#!/bin/bash
# Submit DBS+SEEG combined PCA permutation jobs for all patients × all conditions.
#
# DBS regions (7): ALL, VCVS, LVCVS, RVCVS, SCC, LSCC, RSCC
#
# SEEG single conditions (10): 5 regions × left/right hemispheres
# SEEG pair conditions (25)  : all ordered (left_A + right_B) 5×5 region pairs
#
# Total: 7 patients × 7 DBS × 35 SEEG conditions = 1715 jobs
#
# Usage: bash dbs_seeg_pca_perm_jobs/submit_all.sh

PATIENTS=(
    "DBSTRD001"
    "DBSTRD002"
    "DBSTRD006"
    "DBSTRD008"
    "DBSTRD010"
    "DBSTRD011"
    "DBSTRD014"
)

DBS_REGIONS=("ALL" "VCVS" "LVCVS" "RVCVS" "SCC" "LSCC" "RSCC")
SEEG_REGIONS=("acc" "amy" "dlpfc" "vmpfc" "ofc")
SEEG_HEMISPHERES=("left" "right")

# 25 ordered (left_region, right_region) pairs — all 5×5 combinations
SEEG_PAIRS=(
    "acc acc"
    "acc amy"
    "acc dlpfc"
    "acc vmpfc"
    "acc ofc"
    "amy acc"
    "amy amy"
    "amy dlpfc"
    "amy vmpfc"
    "amy ofc"
    "dlpfc acc"
    "dlpfc amy"
    "dlpfc dlpfc"
    "dlpfc vmpfc"
    "dlpfc ofc"
    "vmpfc acc"
    "vmpfc amy"
    "vmpfc dlpfc"
    "vmpfc vmpfc"
    "vmpfc ofc"
    "ofc acc"
    "ofc amy"
    "ofc dlpfc"
    "ofc vmpfc"
    "ofc ofc"
)

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
echo "Project root: ${PROJECT_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

COUNT=0

for PATIENT in "${PATIENTS[@]}"; do
    for DBS in "${DBS_REGIONS[@]}"; do

        # ── single conditions ──────────────────────────────────────────────────
        for SEEG in "${SEEG_REGIONS[@]}"; do
            for HEMI in "${SEEG_HEMISPHERES[@]}"; do
                JOB_NAME="dbs_seeg_${PATIENT}_${DBS}_${SEEG}_${HEMI}"
                sbatch --job-name="${JOB_NAME}" \
                    --chdir="${PROJECT_DIR}" \
                    "${PROJECT_DIR}/dbs_seeg_pca_perm_jobs/submit_patient.sh" \
                    "${PATIENT}" "${DBS}" "single" "${SEEG}" "${HEMI}"
                COUNT=$((COUNT + 1))
            done
        done

        # ── pair conditions ────────────────────────────────────────────────────
        for PAIR in "${SEEG_PAIRS[@]}"; do
            LEFT=$(echo "$PAIR" | awk '{print $1}')
            RIGHT=$(echo "$PAIR" | awk '{print $2}')
            JOB_NAME="dbs_seeg_${PATIENT}_${DBS}_L${LEFT}_R${RIGHT}"
            sbatch --job-name="${JOB_NAME}" \
                --chdir="${PROJECT_DIR}" \
                "${PROJECT_DIR}/dbs_seeg_pca_perm_jobs/submit_patient.sh" \
                "${PATIENT}" "${DBS}" "pair" "${LEFT}" "${RIGHT}"
            COUNT=$((COUNT + 1))
        done

    done
done

echo ""
echo "Submitted ${COUNT} jobs. Monitor with: squeue -u \$USER"
echo ""
echo "After all jobs finish, merge with:"
for PATIENT in "${PATIENTS[@]}"; do
    for DBS in "${DBS_REGIONS[@]}"; do
        for SEEG in "${SEEG_REGIONS[@]}"; do
            for HEMI in "${SEEG_HEMISPHERES[@]}"; do
                LABEL="dbs_${DBS}_seeg_${SEEG}_${HEMI}"
                echo "  python dbs_seeg_pca_perm_jobs/merge_pca_dbs_seeg_perm.py --patient_id ${PATIENT} --condition_label ${LABEL} --n_chunks 20"
            done
        done
        for PAIR in "${SEEG_PAIRS[@]}"; do
            LEFT=$(echo "$PAIR" | awk '{print $1}')
            RIGHT=$(echo "$PAIR" | awk '{print $2}')
            LABEL="dbs_${DBS}_seeg_L${LEFT}_R${RIGHT}"
            echo "  python dbs_seeg_pca_perm_jobs/merge_pca_dbs_seeg_perm.py --patient_id ${PATIENT} --condition_label ${LABEL} --n_chunks 20"
        done
    done
done
