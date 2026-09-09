# S5 · Partitioned heritability (S-LDSC)

Runs LDSC v1.0.1 partitioned heritability jointly against baselineLD v2.2,
per trait × workflow (Standard vs Suggested), for four cardiovascular GWAS
(Aragam 2022 CAD, Nauffal 2022 QT, Ntalla 2020 PR, Roselli 2025 AF).

Phase numbering below mirrors the on-disk layout at
`~/LDSC_shared/partitioned_h2_v3/` (phases 3–8 native to v3; phase 1/2
scripts reconstructed as noted).

## Scripts (run in order per workflow)

| Order | File | Purpose |
|---|---|---|
| 00 | `00_phase1_bed_ds.py`      | Assemble per-cell-type hg38 BEDs for the **Standard** annotation (Suggested is produced by `s3_gkmqc/05_trim_max_above_for_ldsc.py`) |
| 01 | `01_munge_sumstats.sh`     | Reconstructed from munge logs — 4-trait loop, `--merge-alleles w_hm3.snplist`, `--signed-sumstats BETA,0`, `--N-col N`, `--chunksize 500000` |
| 02 | `02_liftover.sh`           | hg38 → hg19 via **CrossMap 0.7.0**; ±1 kb extension is applied upstream inside phase-1 (in pandas, no `bedtools merge`) |
| 03 | `03_make_annot.sh`         | `make_annot.py` per chromosome, cell type, workflow (feeds phase 4) |
| 03b | `03b_annot_overlap.sh`    | Custom `intersectBed -c` to append cell-type binary columns to baselineLD v2.2 (feeds phase 7) |
| 04 | `04_ldscore.sh`            | `ldsc.py --l2 --thin-annot --ld-wind-cm 1 --print-snps hapmap3_snps/...` (single-annot LD scores) |
| 04b | `04b_ldscore_overlap.sh`  | `ldsc.py --l2 --ld-wind-cm 1 --print-snps 1000G_EUR_Phase3_baseline_snps/...` on the phase-6 overlap annots (no `--thin-annot`) |
| 05 | `05_partitioned_h2.sh`     | Main analysis — `ldsc.py --h2 --overlap-annot --print-coefficients`, joint 10-annotation model per workflow × trait |
| 06 | `06_h2cts_sensitivity.sh`  | `ldsc.py --h2-cts` sensitivity run on all 20 (10 CT × 2 workflow) annotations |

## Notes

- **No `bedtools merge` step**: peaks are extended ±1 kb in pandas
  (`start - 1000` clipped to 0; `end + 1000`) inside `00_phase1_bed_ds.py`
  (Standard) and `s3_gkmqc/05_trim_max_above_for_ldsc.py` (Suggested). This
  contradicts an earlier draft of the supplement.
- **Liftover tool = CrossMap 0.7.0** at
  `/data/programs/conda/envs/crossmap/bin/CrossMap`, chain
  `/data/pipelines/crossmap/chain/hg38ToHg19.over.chain.gz`. `02_liftover.sh`
  originally lived at `partitioned_h2_v3_ruleA/run_phase2_crossmap.sh` — the
  v3 main variant has the equivalent hg19 outputs on disk but no separate
  driver script survives in the v3 folder itself.
- **baselineLD v2.2** carries 97 annotation columns (101 total in the
  `.annot.gz`; first 4 are `CHR BP SNP CM`).
- **MM (Suggested) LDSC results in v3** were inherited unchanged from
  `partitioned_h2_v2` via symlinks — only phase 5 (`--h2-cts`) was re-run in
  v3 with the joint 20-annotation CTS file. Standard (DS) phase 3, 4, 6, 7,
  8 were re-run natively in v3.

## Inputs

- Preprocessed sumstats at `~/LDSC_shared/raw_sumstats/<Trait>/*.preproc.tsv`
  (4 traits)
- HapMap3 allele list `/data/pipelines/ldsc/data_ldsc/w_hm3.snplist`
- baselineLD v2.2, 1000G EUR Phase 3 BIM / weights / frq files (installed
  under `/data/pipelines/ldsc/data_ldsc/`)
- Per-cell-type per-workflow hg38 BEDs (Standard from
  `00_phase1_bed_ds.py`; Suggested from
  `s3_gkmqc/05_trim_max_above_for_ldsc.py`)

## Outputs (not tracked)

- Munged sumstats `~/LDSC_shared/munged_sumstats/<Trait>.sumstats.gz`
- Per-chromosome annot / LD-score files
- Per-trait × per-workflow `.results` (LDSC enrichment + `Coefficient_z-score`)

## Environments

- Steps 01, 03, 04, 05, 06 (LDSC / munge): `envs/ldsc.yml` (Python 2.7.18,
  LDSC v1.0.1)
- Step 02 (CrossMap): `envs/crossmap.yml` — **not yet exported**; the
  binary lives at `/data/programs/conda/envs/crossmap/bin/CrossMap` and
  reports `CrossMap 0.7.0`
- Step 03b (`intersectBed`), step 00 / MM phase-1 (`pandas`):
  `envs/snapatac2.yml` (bedtools 2.31.1)
