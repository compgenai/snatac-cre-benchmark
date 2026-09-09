#!/bin/bash
#SBATCH --job-name=ldsc_phase3_annot
#SBATCH --output=log_ldsc_phase3_annot_%j.out
#SBATCH --error=log_ldsc_phase3_annot_%j.err
#SBATCH --time=120:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# ---------------------------------------------------------------------------
# S5 · Phase 3 — per-chr `make_annot.py` for the DS (Standard) annotation.
#
# MM (Suggested) annotations were prepared in an earlier LDSC v2 run and
# are re-used unchanged in v3 via symlinks (see s5_ldsc/README.md).
# ---------------------------------------------------------------------------

set -e

BEDTOOLS_ENV_BIN="${BEDTOOLS_ENV_BIN:-/data/programs/conda/envs/snapatac2/bin}"
export PATH="${BEDTOOLS_ENV_BIN}:${PATH}"

LDSC_ROOT="${LDSC_ROOT:-/data/pipelines/ldsc}"
LDSC_ENV="${LDSC_ENV:-${LDSC_ROOT}/conda-env}"
LDSC_PY2="${LDSC_PY2:-${LDSC_ENV}/bin/python}"
LDSC_DATA="${LDSC_DATA:-${LDSC_ROOT}/data_ldsc}"
MAKE_ANNOT="${MAKE_ANNOT:-${LDSC_ROOT}/ldsc/make_annot.py}"
BIM_DIR="${BIM_DIR:-${LDSC_DATA}/1000G_EUR_Phase3_plink}"

LDSC_WORK="${LDSC_WORK:-./LDSC_work/v3}"
IN="${IN:-${LDSC_WORK}/phase2_bed_hg19}"
OUT="${OUT:-${LDSC_WORK}/phase3_annot}"

CELLS="Ventricular_Cardiomyocytes Primitive_Endoderm Fibroblasts Atrial_Cardiomyocytes Endothelial Smooth_Muscle Macrophages Nervous_Cells Myofibroblasts Trophectoderm"

T0=$(date +%s)
echo "=== Phase 3 start: $(date) (DS only; MM symlinked from v2) ==="
echo ""

for ct in ${CELLS}; do
    name="DS_${ct}"
    bed="${IN}/${name}.hg19.bed"
    outdir="${OUT}/${name}"

    if [ -f "${outdir}/${name}.22.annot.gz" ] && [ -f "${outdir}/${name}.1.annot.gz" ]; then
        echo "[SKIP] ${name} (already 22 chr present)"
        continue
    fi

    mkdir -p "${outdir}"
    t0=$(date +%s)
    for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
        "${LDSC_PY2}" "${MAKE_ANNOT}" \
            --bed-file "${bed}" \
            --bimfile "${BIM_DIR}/1000G.EUR.QC.${chr}.bim" \
            --annot-file "${outdir}/${name}.${chr}.annot.gz" \
            > /dev/null 2>&1
    done
    t1=$(date +%s)
    echo "[DONE] ${name}  ($((t1-t0))s)"
done

T1=$(date +%s)
echo ""
echo "=== Phase 3 end: $(date) (total $((T1-T0))s) ==="
