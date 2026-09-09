#!/bin/bash
#SBATCH --job-name=ldsc_phase2_liftover
#SBATCH --output=log_ldsc_phase2_liftover_%j.out
#SBATCH --error=log_ldsc_phase2_liftover_%j.err
#SBATCH --time=30:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=lab

# ---------------------------------------------------------------------------
# S5 · Phase 2 — liftOver hg38 → hg19 with CrossMap 0.7.0
#
# Runs CrossMap on the MM (Suggested) per-cell-type BEDs produced by
# s3_gkmqc/05_trim_max_above_for_ldsc.py, and mirrors the DS (Standard)
# hg19 BEDs from the corresponding LDSC work directory (they were lifted
# in the primary v3 run and did not need to be re-run).
#
# Override paths via env vars, e.g.
#   LDSC_WORK=/path/to/LDSC/v3 bash 02_liftover.sh
# ---------------------------------------------------------------------------

set -e

CROSSMAP="${CROSSMAP_BIN:-/data/programs/conda/envs/crossmap/bin/CrossMap}"
CHAIN="${LIFTOVER_CHAIN:-/data/pipelines/crossmap/chain/hg38ToHg19.over.chain.gz}"
LDSC_WORK="${LDSC_WORK:-./LDSC_work/v3}"

IN="${IN:-${LDSC_WORK}/phase1_bed_hg38}"
OUT="${OUT:-${LDSC_WORK}/phase2_bed_hg19}"
DS_SRC="${DS_SRC:-${LDSC_WORK}/phase2_bed_hg19}"   # DS hg19 BEDs (already lifted in primary v3 run)

CELLS="Atrial_Cardiomyocytes Endothelial Fibroblasts Macrophages Myofibroblasts Nervous_Cells Primitive_Endoderm Smooth_Muscle Trophectoderm Ventricular_Cardiomyocytes"

echo "=== phase2: MM CrossMap start $(date) ==="
mkdir -p "${OUT}/unmapped"
for ct in ${CELLS}; do
  in_bed="${IN}/MM_${ct}.bed"
  out_bed="${OUT}/MM_${ct}.hg19.bed"
  t0=$(date +%s)
  "${CROSSMAP}" bed "${CHAIN}" "${in_bed}" "${out_bed}" > "${OUT}/unmapped/MM_${ct}.crossmap.log" 2>&1
  # CrossMap writes unmapped to <out>.unmap
  [ -f "${out_bed}.unmap" ] && mv "${out_bed}.unmap" "${OUT}/unmapped/MM_${ct}.unmap"
  sort -k1,1 -k2,2n "${out_bed}" -o "${out_bed}"
  n_in=$(wc -l < "${in_bed}"); n_out=$(wc -l < "${out_bed}"); t1=$(date +%s)
  echo "[DONE] MM_${ct}  hg38=${n_in}  hg19=${n_out}  ($((t1-t0))s)"
done

echo ""
echo "=== phase2: DS mirror from primary v3 (no re-run) ==="
for ct in ${CELLS}; do
  ln -sfn "${DS_SRC}/DS_${ct}.hg19.bed" "${OUT}/DS_${ct}.hg19.bed"
  echo "[LINK] DS_${ct}.hg19.bed"
done
echo "=== phase2 end $(date) ==="
