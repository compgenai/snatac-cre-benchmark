#!/bin/bash
#SBATCH --job-name=ldsc_phase4_ldscore
#SBATCH --output=log_ldsc_phase4_ldscore_%A_%a.out
#SBATCH --error=log_ldsc_phase4_ldscore_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-10

# ---------------------------------------------------------------------------
# S5 · Phase 4 — single-annotation LD scores (`ldsc.py --l2 --thin-annot`).
#
# One SLURM array task per cell type (10 total). DS only in the v3 rerun.
# ---------------------------------------------------------------------------

set -e

BEDTOOLS_ENV_BIN="${BEDTOOLS_ENV_BIN:-/data/programs/conda/envs/snapatac2/bin}"
export PATH="${BEDTOOLS_ENV_BIN}:${PATH}"

LDSC_ROOT="${LDSC_ROOT:-/data/pipelines/ldsc}"
LDSC_ENV="${LDSC_ENV:-${LDSC_ROOT}/conda-env}"
LDSC_PY2="${LDSC_PY2:-${LDSC_ENV}/bin/python}"
LDSC_BIN="${LDSC_BIN:-${LDSC_ROOT}/ldsc/ldsc.py}"
LDSC_DATA="${LDSC_DATA:-${LDSC_ROOT}/data_ldsc}"
BFILE_DIR="${BFILE_DIR:-${LDSC_DATA}/1000G_EUR_Phase3_plink}"
HM3_DIR="${HM3_DIR:-${LDSC_DATA}/hapmap3_snps}"

LDSC_WORK="${LDSC_WORK:-./LDSC_work/v3}"

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

ANNOT_DIR="${LDSC_WORK}/phase3_annot/${NAME}"
OUT_DIR="${LDSC_WORK}/phase4_ldscore/${NAME}"
mkdir -p "${OUT_DIR}"

T0=$(date +%s)
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
    out="${OUT_DIR}/${NAME}.${chr}"
    if [ -f "${out}.l2.ldscore.gz" ] && [ -f "${out}.l2.M_5_50" ]; then
        echo "[SKIP] chr${chr}"
        continue
    fi
    t0=$(date +%s)
    "${LDSC_PY2}" "${LDSC_BIN}" --l2 \
        --bfile "${BFILE_DIR}/1000G.EUR.QC.${chr}" \
        --ld-wind-cm 1 \
        --annot "${ANNOT_DIR}/${NAME}.${chr}.annot.gz" \
        --thin-annot \
        --out "${out}" \
        --print-snps "${HM3_DIR}/hm.${chr}.snp" \
        > /dev/null 2>&1
    t1=$(date +%s)
    echo "[DONE] chr${chr} ($((t1-t0))s)"
done

T1=$(date +%s)
echo "===== ${NAME} total: $((T1-T0))s ====="
echo "End: $(date)"
