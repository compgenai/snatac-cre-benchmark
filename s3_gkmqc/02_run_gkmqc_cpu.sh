#!/bin/bash
#SBATCH --job-name=gkmqc_cpu
#SBATCH --output=log_gkmqc_cpu_%x_%j.out
#SBATCH --error=log_gkmqc_cpu_%x_%j.err
#SBATCH --time=80:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --account=lab
#SBATCH --export=NONE

# ---------------------------------------------------------------------------
# gkmQC (CPU) evaluate — Branch A (feeds Fig 4a / Table S2)
#
# Consolidated from 10 per-cell-type SLURM scripts originally at:
#   output_seung/gkmQC/Sug_dedup_50pct/job_gkmqc_Improved_<CELLTYPE>.dedup_50pct.sh
#
# The originals differed only in the CELLTYPE literal in 5 lines
# (job-name, --output, --error, echo line, -i, -n).  This template accepts
# CELLTYPE as $1; SBATCH directives are otherwise byte-identical.
#
# Usage (submit once per cell type):
#   for CT in Atrial_Cardiomyocytes Endothelial Fibroblasts Macrophages \
#             Myofibroblasts Nervous_Cells Primitive_Endoderm \
#             Smooth_Muscle Trophectoderm Ventricular_Cardiomyocytes; do
#     sbatch --job-name=gkmqc_cpu_${CT} 02_run_gkmqc_cpu.sh ${CT}
#   done
#
# Working directory must contain the dedup_50pct narrowPeak inputs, since -i
# is relative ("./Improved_<CT>.dedup_50pct.narrowPeak") — matches the
# original job-script behavior (all 10 originals used the same relative path
# and were submitted from output_seung/gkmQC/Sug_dedup_50pct/).
# ---------------------------------------------------------------------------

CELLTYPE="${1:?Usage: sbatch $0 <CELLTYPE>}"

CONDA_EXE="${CONDA_EXE:-/data/programs/conda/bin/conda}"
GKMQC_ENV="${GKMQC_ENV:-/data/pipelines/gkmqc/conda-env}"
GKMQC_BIN="${GKMQC_BIN:-/data/pipelines/gkmqc/bin/gkmqc.py}"

eval "$(${CONDA_EXE} shell.bash hook)"
conda activate "${GKMQC_ENV}"

echo "[INFO] Running gkmqc (CPU) for Improved_${CELLTYPE}.dedup_50pct"
python "${GKMQC_BIN}" evaluate \
  -i "./Improved_${CELLTYPE}.dedup_50pct.narrowPeak" \
  -g "hg38" \
  -n "Improved_${CELLTYPE}.dedup_50pct" \
  -@ "100" \
  -o "8" \
  -re "40"
