# Inputs

Nothing in this directory is tracked by git.

## Expected layout

```
data/
├── encode/            snATAC-seq fragment files, one per sample (Table S1 accessions)
├── reference/         refdata-cellranger-arc-GRCh38-2020-A-2.0.0
├── gwas/              GWAS summary statistics (Table S3)
└── ldsc/              baselineLD v2.2, 1000G EUR Phase 3, HapMap3 allele list, weights
```

## How to obtain

- **ENCODE** — download by accession from https://www.encodeproject.org (Table S1).
- **GWAS** — obtain from the sources cited in Table S3.
- **LDSC reference files** — from the LDSC distribution (baselineLD v2.2, 1000G EUR
  Phase 3 frq files, weights excluding HLA).
