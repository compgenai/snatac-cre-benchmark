#!/bin/bash
#SBATCH --job-name=ldsc_phase6_annot_overlap
#SBATCH --output=log_ldsc_phase6_annot_overlap_%j.out
#SBATCH --error=log_ldsc_phase6_annot_overlap_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# ---------------------------------------------------------------------------
# S5 · Phase 6 — build the overlap annotation for the joint LDSC model.
#
# Appends a per-cell-type binary column (peak-overlap indicator) to the
# 97-column baselineLD v2.2 annotation using `intersectBed -c`. DS only in
# the v3 rerun — MM inherited from v2 via symlink.
# ---------------------------------------------------------------------------

set -e

BEDTOOLS_ENV_BIN="${BEDTOOLS_ENV_BIN:-/data/programs/conda/envs/snapatac2/bin}"
export PATH="${BEDTOOLS_ENV_BIN}:${PATH}"

LDSC_ROOT="${LDSC_ROOT:-/data/pipelines/ldsc}"
LDSC_DATA="${LDSC_DATA:-${LDSC_ROOT}/data_ldsc}"
BASELINE_DIR="${BASELINE_DIR:-${LDSC_DATA}/1000G_Phase3_baselineLD_v2.2}"

LDSC_WORK="${LDSC_WORK:-./LDSC_work/v3}"
PHASE2="${PHASE2:-${LDSC_WORK}/phase2_bed_hg19}"
PHASE6="${PHASE6:-${LDSC_WORK}/phase6_annot_overlap}"

CELLS="Ventricular_Cardiomyocytes Primitive_Endoderm Fibroblasts Atrial_Cardiomyocytes Endothelial Smooth_Muscle Macrophages Nervous_Cells Myofibroblasts Trophectoderm"

T0=$(date +%s)
echo "=== Phase 6 start: $(date) (DS only; MM symlinked from v2) ==="

for ct in ${CELLS}; do
    EXPID="DS_${ct}"
    CELL_BED="${PHASE2}/${EXPID}.hg19.bed"
    OUT_DIR="${PHASE6}/${EXPID}"

    if [ -f "${OUT_DIR}/${EXPID}.22.annot.gz" ] && [ -f "${OUT_DIR}/${EXPID}.1.annot.gz" ]; then
        echo "[SKIP] ${EXPID}"
        continue
    fi

    mkdir -p "${OUT_DIR}"
    t0=$(date +%s)
    for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
        ANNOT_BL="${BASELINE_DIR}/baselineLD.${chr}.annot.gz"
        OUT="${OUT_DIR}/${EXPID}.${chr}.annot.gz"
        TMP="${OUT_DIR}/tmp.${chr}"

        zcat "${ANNOT_BL}" | awk -v OFS="\t" 'NR>1{print "chr"$1, $2-1, $2, $3}' \
            | intersectBed -c -a - -b "${CELL_BED}" \
            | awk -v EXPID="${EXPID}" 'BEGIN{print EXPID} {print $5}' > "${TMP}"

        zcat "${ANNOT_BL}" | cut -f 1-4 | paste - "${TMP}" | gzip > "${OUT}"
        rm -f "${TMP}"
    done
    t1=$(date +%s)
    echo "[DONE] ${EXPID} ($((t1-t0))s)"
done

T1=$(date +%s)
echo "=== Phase 6 end: $(date) (total $((T1-T0))s) ==="
