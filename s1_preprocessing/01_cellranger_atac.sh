#!/bin/bash
# ---------------------------------------------------------------------------
# S1 · Cell Ranger ATAC — per-sample count
#
# Reconstructed from the 9 per-sample `_cmdline` files under
#   output_so/cell_ranger_output_heart/sample<N>/_cmdline
# All 9 samples used byte-identical command lines that differed only in the
# `sample<N>` label passed to --id and --fastqs.
#
# Binary:    cellranger-atac 2.1.0  (ATAC-only pipeline)
# Chemistry: ARC-v1                 (multiome ATAC library format)
# Reference: refdata-cellranger-arc-GRCh38-2020-A-2.0.0  (ARC reference)
#
# Override binary and reference locations via env vars:
#   CELLRANGER_BIN=/opt/cellranger-atac-2.1.0/cellranger-atac \
#   CELLRANGER_REF=/refs/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
#   FASTQ_ROOT=./data \
#   bash 01_cellranger_atac.sh
# ---------------------------------------------------------------------------

set -euo pipefail

CELLRANGER_BIN="${CELLRANGER_BIN:-cellranger-atac}"
CELLRANGER_REF="${CELLRANGER_REF:?export CELLRANGER_REF=/path/to/refdata-cellranger-arc-GRCh38-2020-A-2.0.0}"
FASTQ_ROOT="${FASTQ_ROOT:-./data}"

SAMPLES=(sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8 sample9)

for S in "${SAMPLES[@]}"; do
  echo "[INFO] $(date '+%F %T')  cellranger-atac count --id ${S}"
  "${CELLRANGER_BIN}" count \
    --id "${S}" \
    --reference "${CELLRANGER_REF}" \
    --fastqs "${FASTQ_ROOT}/${S}" \
    --chemistry=ARC-v1
done
