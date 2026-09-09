#!/bin/bash
#SBATCH --job-name=ldsc_phase7_ldscore_overlap
#SBATCH --output=log_ldsc_phase7_ldscore_overlap_%A_%a.out
#SBATCH --error=log_ldsc_phase7_ldscore_overlap_%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-10

# ---------------------------------------------------------------------------
# S5 · Phase 7 — joint LD scores on the overlap annotation (no --thin-annot).
# ---------------------------------------------------------------------------

set -e

BEDTOOLS_ENV_BIN="${BEDTOOLS_ENV_BIN:-/data/programs/conda/envs/snapatac2/bin}"
export PATH="${BEDTOOLS_ENV_BIN}:${PATH}"

CONDA_EXE="${CONDA_EXE:-/data/programs/conda/bin/conda}"
LDSC_ROOT="${LDSC_ROOT:-/data/pipelines/ldsc}"
LDSC_ENV="${LDSC_ENV:-${LDSC_ROOT}/conda-env}"
LDSC_BIN="${LDSC_BIN:-${LDSC_ROOT}/ldsc/ldsc.py}"
LDSC_DATA="${LDSC_DATA:-${LDSC_ROOT}/data_ldsc}"
BFILE_DIR="${BFILE_DIR:-${LDSC_DATA}/1000G_EUR_Phase3_plink}"
SNPLIST_DIR="${SNPLIST_DIR:-${LDSC_DATA}/1000G_EUR_Phase3_baseline_snps}"

LDSC_WORK="${LDSC_WORK:-./LDSC_work/v3}"

source "$(dirname ${CONDA_EXE})/activate" "${LDSC_ENV}"

# DS only (MM unchanged from v2, symlinked)
TASKS=(
    "_placeholder_"
    "DS_Ventricular_Cardiomyocytes"
    "DS_Atrial_Cardiomyocytes"
    "DS_Fibroblasts"
    "DS_Endothelial"
    "DS_Smooth_Muscle"
    "DS_Macrophages"
    "DS_Nervous_Cells"
    "DS_Myofibroblasts"
    "DS_Primitive_Endoderm"
    "DS_Trophectoderm"
)

NAME=${TASKS[${SLURM_ARRAY_TASK_ID}]}
echo "===== Task ${SLURM_ARRAY_TASK_ID} : ${NAME} ====="
echo "Start: $(date)"

PHASE6="${LDSC_WORK}/phase6_annot_overlap/${NAME}"
PHASE7="${LDSC_WORK}/phase7_ldscore_overlap/${NAME}"
mkdir -p "${PHASE7}"

export OPENBLAS_NUM_THREADS=2

T0=$(date +%s)
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
    out="${PHASE7}/${NAME}.${chr}"
    if [ -f "${out}.l2.ldscore.gz" ] && [ -f "${out}.l2.M_5_50" ]; then
        echo "[SKIP] chr${chr}"
        continue
    fi
    t0=$(date +%s)
    python "${LDSC_BIN}" \
        --l2 \
        --bfile "${BFILE_DIR}/1000G.EUR.QC.${chr}" \
        --ld-wind-cm 1 \
        --annot "${PHASE6}/${NAME}.${chr}.annot.gz" \
        --out "${out}" \
        --print-snps "${SNPLIST_DIR}/baseline_snp.${chr}.snp" \
        > /dev/null 2>&1
    t1=$(date +%s)
    echo "[DONE] chr${chr} ($((t1-t0))s)"
done

T1=$(date +%s)
echo "===== ${NAME} total: $((T1-T0))s ====="
echo "End: $(date)"
