# envs/

Conda environment specifications for the four (four+one) environments the
pipeline runs in. Exported with `conda env export --no-builds`.

| File | Purpose | Original conda location |
|------|---------|-------------------------|
| `snapatac2.yml` | S1 preprocessing, S2 peak calling (Python parts), S3 dedup, S4 bedtools | `~/.conda/envs/sohyeong-snapatac2` (Python 3.9.21; snapatac2 2.8.0, harmonypy 0.0.10, scanpy 1.10.3, bedtools 2.31.1) |
| `macs2.yml`     | S2 MACS2 peak calling | `/data/pipelines/atac_bulk/conda-env` (MACS2 2.2.6) |
| `gkmqc.yml`     | S3 gkmQC evaluate (CPU and GPU wrappers activate the same env) | `/data/pipelines/gkmqc/conda-env` (Python 3.7.12; numpy 1.21.6, scikit-learn 1.0.2, scipy 1.7.3, pyfasta 0.5.2; gkmQC v1.0.0) |
| `ldsc.yml`      | S5 munge, LDSC phases 3–8 | `/data/pipelines/ldsc/conda-env` (Python 2.7.18; LDSC v1.0.1, pybedtools 0.8.1, bedtools 2.30.0) |

## CrossMap

CrossMap 0.7.0 lives in a separate env at
`/data/programs/conda/envs/crossmap`. Not yet exported here — export
similarly with `conda env export --prefix /data/programs/conda/envs/crossmap
--no-builds > envs/crossmap.yml` if needed to reproduce S5 · phase 2.

## gkmQC install

The gkmQC binaries at `/data/pipelines/gkmqc/bin/{gkmqc.py, gkmqc-gpu.py}`
are not distributed with a git tag or commit hash (source dir is not a git
checkout). Install artifacts are dated 2024-12-24 per `CODE_MANUAL.md` in
the install tree. Reproducing this branch requires the same install; see
the manuscript's Software table for the version string reported by
`gkmqc.py evaluate --version`.
