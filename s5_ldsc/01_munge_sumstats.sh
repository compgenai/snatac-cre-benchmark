#!/bin/bash
# ---------------------------------------------------------------------------
# S5 · Munge GWAS summary statistics with LDSC's munge_sumstats.py
#
# Reconstructed verbatim from the four munge logs at:
#   ~/LDSC_shared/munged_sumstats/{Aragam_2022_CAD_primary,
#                                  Nauffal_2022_QT,
#                                  Ntalla_2020_PR_EUR,
#                                  Roselli_2025_AF_common}.log
# All four traits were munged with identical flag sets — only --sumstats /
# --out paths differ per trait.
#
# Chunksize 500,000 matches the log-quoted setting; --N-col N assumes the
# preprocessed .tsv exposes a per-SNP N column.
# ---------------------------------------------------------------------------

set -euo pipefail

# LDSC install (Python 2.7.18 env)
CONDA_EXE="${CONDA_EXE:-/data/programs/conda/bin/conda}"
LDSC_ROOT="${LDSC_ROOT:-/data/pipelines/ldsc}"
LDSC_ENV="${LDSC_ENV:-${LDSC_ROOT}/conda-env}"
MUNGE="${MUNGE:-${LDSC_ROOT}/munge_sumstats.py}"     # LDSC's munge script
HAPMAP3="${HAPMAP3:-${LDSC_ROOT}/data_ldsc/w_hm3.snplist}"   # HapMap3 allele list (n=1,217,311)

LDSC_SHARED="${LDSC_SHARED:-./LDSC_shared}"
RAW_DIR="${RAW_DIR:-${LDSC_SHARED}/raw_sumstats}"
OUT_DIR="${OUT_DIR:-${LDSC_SHARED}/munged_sumstats}"
mkdir -p "${OUT_DIR}"

eval "$(${CONDA_EXE} shell.bash hook)"
conda activate "${LDSC_ENV}"

# TRAIT_KEY = <output stem>
# INPUT_TSV = preprocessed sumstats fed to --sumstats
TRAITS=(
  "Aragam_2022_CAD_primary    ${RAW_DIR}/Aragam_2022_CAD/Aragam_2022_CAD_primary.preproc.tsv"
  "Nauffal_2022_QT             ${RAW_DIR}/Nauffal_2022_QT/Nauffal_2022_QT.preproc.tsv"
  "Ntalla_2020_PR_EUR          ${RAW_DIR}/Ntalla_2020_PR_EUR/Ntalla_2020_PR_EUR.preproc.tsv"
  "Roselli_2025_AF_common      ${RAW_DIR}/Roselli_2025_AF/Roselli_2025_AF_common.preproc.tsv"
)

for row in "${TRAITS[@]}"; do
  read -r TRAIT INPUT <<<"${row}"
  echo "[INFO] $(date '+%F %T')  munging ${TRAIT}"
  python "${MUNGE}" \
    --sumstats       "${INPUT}" \
    --out            "${OUT_DIR}/${TRAIT}" \
    --merge-alleles  "${HAPMAP3}" \
    --signed-sumstats BETA,0 \
    --N-col N \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p  P \
    --chunksize 500000
done
