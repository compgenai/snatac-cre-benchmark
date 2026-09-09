#!/bin/bash
#SBATCH --job-name=ldsc_phase5_h2cts
#SBATCH --output=log_ldsc_phase5_h2cts_%j.out
#SBATCH --error=log_ldsc_phase5_h2cts_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

# ---------------------------------------------------------------------------
# S5 · Phase 5 — `--h2-cts` sensitivity run on the 20-annotation CTS file
# (10 cell types × 2 workflows).
# ---------------------------------------------------------------------------

set -e

BEDTOOLS_ENV_BIN="${BEDTOOLS_ENV_BIN:-/data/programs/conda/envs/snapatac2/bin}"
export PATH="${BEDTOOLS_ENV_BIN}:${PATH}"

LDSC_ROOT="${LDSC_ROOT:-/data/pipelines/ldsc}"
LDSC_ENV="${LDSC_ENV:-${LDSC_ROOT}/conda-env}"
LDSC_PY2="${LDSC_PY2:-${LDSC_ENV}/bin/python}"
LDSC_BIN="${LDSC_BIN:-${LDSC_ROOT}/ldsc/ldsc.py}"
LDSC_DATA="${LDSC_DATA:-${LDSC_ROOT}/data_ldsc}"
BASELINE_LD="${BASELINE_LD:-${LDSC_DATA}/1000G_Phase3_baselineLD_v2.2/baselineLD.}"
WEIGHTS="${WEIGHTS:-${LDSC_DATA}/weights_hm3_no_hla/weights.}"

LDSC_WORK="${LDSC_WORK:-./LDSC_work/v3}"
LDSC_SHARED="${LDSC_SHARED:-./LDSC_shared}"
SUMSTATS_DIR="${SUMSTATS_DIR:-${LDSC_SHARED}/munged_sumstats}"
CTS_FILE="${CTS_FILE:-${LDSC_WORK}/phase5_h2cts/cts_20_all.txt}"
OUT_DIR="${OUT_DIR:-${LDSC_WORK}/phase5_h2cts}"

TRAITS="Roselli_2025_AF_common Aragam_2022_CAD_primary Nauffal_2022_QT Ntalla_2020_PR_EUR"

T0=$(date +%s)
echo "=== Phase 5 start: $(date) ==="
echo ""

for trait in ${TRAITS}; do
    t0=$(date +%s)
    echo "----- ${trait} -----"
    "${LDSC_PY2}" "${LDSC_BIN}" \
        --h2-cts "${SUMSTATS_DIR}/${trait}.sumstats.gz" \
        --ref-ld-chr "${BASELINE_LD}" \
        --ref-ld-chr-cts "${CTS_FILE}" \
        --w-ld-chr "${WEIGHTS}" \
        --out "${OUT_DIR}/${trait}_h2cts" \
        > "${OUT_DIR}/${trait}_h2cts.console.log" 2>&1
    t1=$(date +%s)
    echo "[DONE] ${trait} ($((t1-t0))s)"
done

T1=$(date +%s)
echo ""
echo "=== Phase 5 end: $(date) (total $((T1-T0))s) ==="
