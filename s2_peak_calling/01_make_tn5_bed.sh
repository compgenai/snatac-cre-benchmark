#!/bin/bash
#SBATCH --job-name=tn5_extraction
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=tn5_extraction_%j.log
#SBATCH --error=tn5_extraction_%j.err

# ---------------------------------------------------------------------------
# S2 · Tn5 insertion extraction
#
# Reads per-(sample × cell type) fragment BEDs (`mcluster<CT>.<SAMPLE>_fragments.bed`)
# and writes 1-bp insertion events (both fragment ends, no +4/-5 shift).
#
# Override paths via env vars, e.g.
#   DATA_ROOT=/path/to/output_so sbatch 01_make_tn5_bed.sh
# ---------------------------------------------------------------------------

DATA_ROOT="${DATA_ROOT:-./output_so}"
INPUT_DIR="${INPUT_DIR:-${DATA_ROOT}/fragments}"
OUT_DIR="${OUT_DIR:-${DATA_ROOT}/tn5_insertions}"

mkdir -p "${OUT_DIR}"

echo "[START] Tn5 insertion site extraction"
echo "INPUT_DIR: ${INPUT_DIR}"
echo "OUT_DIR:   ${OUT_DIR}"

for f in "${INPUT_DIR}"/mcluster*_fragments.bed; do
  fname=$(basename "${f}")
  # mcluster<CELLTYPE>.<SAMPLE>_fragments.bed → mcluster<CELLTYPE>.<SAMPLE>
  prefix=${fname%_fragments.bed}

  TN5_BED="${OUT_DIR}/${prefix}_tn5_ins.bed"

  awk 'BEGIN{OFS="\t"}{
    # fragment start insertion
    print $1, $2, $2+1, $4, $5
    # fragment end insertion
    print $1, $3-1, $3, $4, $5
  }' "${f}" > "${TN5_BED}"

  echo "[OK] Generated: ${TN5_BED}"
done

echo "[COMPLETE] All files processed"
