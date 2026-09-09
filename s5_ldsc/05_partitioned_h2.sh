#!/bin/bash
#SBATCH --job-name=ldsc_phase8_h2_multi
#SBATCH --output=log_ldsc_phase8_h2_multi_%A_%a.out
#SBATCH --error=log_ldsc_phase8_h2_multi_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --array=1-4

# ---------------------------------------------------------------------------
# S5 · Phase 8 — main analysis: `ldsc.py --h2 --overlap-annot --print-coefficients`
# Joint 10-annotation model per trait per workflow.
# ---------------------------------------------------------------------------

set -e

BEDTOOLS_ENV_BIN="${BEDTOOLS_ENV_BIN:-/data/programs/conda/envs/snapatac2/bin}"
export PATH="${BEDTOOLS_ENV_BIN}:${PATH}"

CONDA_EXE="${CONDA_EXE:-/data/programs/conda/bin/conda}"
LDSC_ROOT="${LDSC_ROOT:-/data/pipelines/ldsc}"
LDSC_ENV="${LDSC_ENV:-${LDSC_ROOT}/conda-env}"
LDSC_BIN="${LDSC_BIN:-${LDSC_ROOT}/ldsc/ldsc.py}"
LDSC_DATA="${LDSC_DATA:-${LDSC_ROOT}/data_ldsc}"
LDSC_WORK="${LDSC_WORK:-./LDSC_work/v3}"
LDSC_SHARED="${LDSC_SHARED:-./LDSC_shared}"

source "$(dirname ${CONDA_EXE})/activate" "${LDSC_ENV}"

# DS only (MM phase 8 results unchanged - symlinked from v2)
# 4 traits × 1 method (DS)
TASKS=(
    "_placeholder_"
    "Roselli_2025_AF_common,DS"
    "Aragam_2022_CAD_primary,DS"
    "Nauffal_2022_QT,DS"
    "Ntalla_2020_PR_EUR,DS"
)

TASK=${TASKS[${SLURM_ARRAY_TASK_ID}]}
TRAIT=$(echo "${TASK}" | cut -d, -f1)
METHOD=$(echo "${TASK}" | cut -d, -f2)
echo "===== Task ${SLURM_ARRAY_TASK_ID} : ${TRAIT} × ${METHOD} ====="
echo "Start: $(date)"

CELLS="Ventricular_Cardiomyocytes Atrial_Cardiomyocytes Fibroblasts Endothelial Smooth_Muscle Macrophages Nervous_Cells Myofibroblasts Primitive_Endoderm Trophectoderm"

REF_LD_CHR="${LDSC_DATA}/1000G_Phase3_baselineLD_v2.2/baselineLD."
for ct in ${CELLS}; do
    NAME="${METHOD}_${ct}"
    REF_LD_CHR="${REF_LD_CHR},${LDSC_WORK}/phase7_ldscore_overlap/${NAME}/${NAME}."
done

SUMSTATS="${LDSC_SHARED}/munged_sumstats/${TRAIT}.sumstats.gz"
OUT_DIR="${LDSC_WORK}/phase8_h2_overlap"
mkdir -p "${OUT_DIR}"

echo "REF_LD_CHR (first 200 chars): ${REF_LD_CHR:0:200}..."
echo ""

T0=$(date +%s)
python "${LDSC_BIN}" \
    --h2 "${SUMSTATS}" \
    --w-ld-chr "${LDSC_DATA}/weights_hm3_no_hla/weights." \
    --ref-ld-chr "${REF_LD_CHR}" \
    --overlap-annot \
    --frqfile-chr "${LDSC_DATA}/1000G_Phase3_frq/1000G.EUR.QC." \
    --print-coefficients \
    --out "${OUT_DIR}/${METHOD}_${TRAIT}"
T1=$(date +%s)

echo ""
echo "===== ${TRAIT} × ${METHOD} total: $((T1-T0))s ====="
echo "End: $(date)"
