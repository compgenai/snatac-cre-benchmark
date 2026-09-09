#!/bin/bash
# ---------------------------------------------------------------------------
# S2 · MACS2 peak calling — Suggested workflow (permissive `-p 0.01`)
#
# Identical to `02_macs2_standard.sh` except for the additional `-p 0.01`
# flag and a different OUTPUT_BASE default. Both scripts read the same Tn5
# insertion BEDs.
# ---------------------------------------------------------------------------

DATA_ROOT="${DATA_ROOT:-./output_so}"
INPUT_DIR="${INPUT_DIR:-${DATA_ROOT}/tn5_insertions}"
OUTPUT_BASE="${OUTPUT_BASE:-${DATA_ROOT}/peakcalling_p0.01}"
LOG_BASE="${OUTPUT_BASE}/_logs_local"

CONDA_EXE="${CONDA_EXE:-/data/programs/conda/bin/conda}"
MACS2_ENV="${MACS2_ENV:-/data/pipelines/atac_bulk/conda-env}"
SLURM_PARTITION="${SLURM_PARTITION:-cpu}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-lab}"
THREADS="${THREADS:-8}"
PREV_JOBID=""

mkdir -p "${LOG_BASE}"

FILES=( "${INPUT_DIR}"/*.bed )

for BED in "${FILES[@]}"; do
  [ -f "${BED}" ] || continue

  BASE=$(basename "${BED}" .bed)
  JOB_SCRIPT="job_macs2_${BASE}.sh"

  echo "[INFO] Preparing MACS2 job for ${BASE}"

  cat <<EOF > "${JOB_SCRIPT}"
#!/bin/bash
#SBATCH --job-name=macs2_${BASE}
#SBATCH --output=${LOG_BASE}/slurm_${BASE}.out
#SBATCH --error=${LOG_BASE}/slurm_${BASE}.err
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${THREADS}
#SBATCH --mem=16G
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --export=NONE

eval "\$(${CONDA_EXE} shell.bash hook)"
conda activate ${MACS2_ENV}

OUTDIR="${OUTPUT_BASE}/${BASE}/macs2_${BASE}"
mkdir -p "\$OUTDIR"

macs2 callpeak \\
  -n "${BASE}" \\
  -g hs \\
  -p 0.01 \\
  --nomodel \\
  --shift -50 --extsize 100 \\
  --keep-dup all \\
  -f BED \\
  -t "${BED}" \\
  -B --SPMR \\
  --buffer-size 1000000 \\
  --outdir "\$OUTDIR"
EOF

  if [ -z "${PREV_JOBID}" ]; then
    JOBID=$(sbatch "${JOB_SCRIPT}" | awk '{print $4}')
  else
    JOBID=$(sbatch --dependency=afterok:${PREV_JOBID} "${JOB_SCRIPT}" | awk '{print $4}')
  fi

  echo "[INFO] Submitted MACS2 job → JobID=${JOBID}"
  PREV_JOBID=${JOBID}
done

echo "[DONE] All MACS2 jobs submitted."
