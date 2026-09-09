#!/bin/bash
#SBATCH --job-name=gkmqc_cpu_ldsc
#SBATCH --output=log_gkmqc_cpu_ldsc_%x_%j.out
#SBATCH --error=log_gkmqc_cpu_ldsc_%x_%j.err
#SBATCH --time=80:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --account=lab
#SBATCH --export=NONE

# ---------------------------------------------------------------------------
# gkmQC (CPU) evaluate — Branch B (feeds LDSC Suggested annotation)
#
# Runs the SAME CPU gkmQC build as Branch A (`02_run_gkmqc_cpu.sh`) but with
# two differences that define Branch B:
#   • Input   : raw MM 3D-centroid narrowPeak (BEFORE 50% reciprocal dedup)
#   • -re 350 : evaluate up to 350 subsets (default is 40 for Branch A)
#
# Historical note: on the machine that produced the manuscript's LDSC v3
# results, this branch was originally run with the GPU build
# (`gkmqc-gpu.py -P slurm_gpu`) for wall-clock reasons. Repo standardizes on
# the CPU build for reproducibility (no GPU dependency). The two builds
# produce equivalent .eval.out AUC profiles; only wall time differs.
#
# Usage (submit once per cell type):
#   for CT in Atrial_Cardiomyocytes Endothelial Fibroblasts Macrophages \
#             Myofibroblasts Nervous_Cells Primitive_Endoderm \
#             Smooth_Muscle Trophectoderm Ventricular_Cardiomyocytes; do
#     sbatch --job-name=gkmqc_cpu_ldsc_${CT} 04_run_gkmqc_cpu_for_ldsc.sh ${CT}
#   done
#
# Working directory must contain the raw MM centroid narrowPeak files
# (Improved_<CT>.narrowPeak).
# ---------------------------------------------------------------------------

CELLTYPE="${1:?Usage: sbatch $0 <CELLTYPE>}"

CONDA_EXE="${CONDA_EXE:-/data/programs/conda/bin/conda}"
GKMQC_ENV="${GKMQC_ENV:-/data/pipelines/gkmqc/conda-env}"
GKMQC_BIN="${GKMQC_BIN:-/data/pipelines/gkmqc/bin/gkmqc.py}"

eval "$(${CONDA_EXE} shell.bash hook)"
conda activate "${GKMQC_ENV}"

echo "[INFO] Running gkmqc (CPU) for Improved_${CELLTYPE} (raw MM centroid, LDSC branch)"
python "${GKMQC_BIN}" evaluate \
  -i "./Improved_${CELLTYPE}.narrowPeak" \
  -g "hg38" \
  -n "Improved_${CELLTYPE}" \
  -@ "100" \
  -o "8" \
  -re "350"
