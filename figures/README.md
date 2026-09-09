# figures/

Figure-producing scripts for the manuscript.

## Scripts

| File | Produces | Original path |
|------|----------|---------------|
| `fig4_peak_overlap.py`     | Figure 4a bar chart (Standard vs Suggested peak counts and asymmetric concordance) | `/data/shared/sohyeong_data/figures/figure_benchmark/make_figure_v3_split_A_B.py` |
| `fig3_ldsc_enrichment.py`  | Figure 3 heritability enrichment heatmap by phenotype, ordered by count | `~/LDSC_shared/partitioned_h2_v3/make_enrichment_heatmap_by_phenotype_ordered_by_count.py` |
| `fig3_tau_zscore.py`       | Supplementary τ z-score (`Coefficient_z-score`) heatmap | `~/LDSC_shared/partitioned_h2_v3/make_coef_zscore_heatmap_by_phenotype.py` |
| `make_master_table.py`     | Consolidates all `phase8_h2_overlap/*.results` into `LDSC_master_table.csv` (input to the two heatmap scripts above) | `~/LDSC_shared/partitioned_h2_v3/make_master_table.py` |

## Inputs

- Fig 4a: per-cell-type overlap TSV produced by `s4_peak_comparison/01_concordance.py`
- Fig 3 / τ heatmap: `LDSC_master_table.csv` produced by `make_master_table.py`, which itself reads `s5_ldsc/05_partitioned_h2.sh` outputs

## Environment

Figures use matplotlib and pandas. Any of `envs/snapatac2.yml` or a plain
Python 3 environment with `matplotlib`, `pandas`, `numpy` suffices for the
plotting scripts. `make_master_table.py` requires the LDSC output tree to
be present but does not invoke LDSC itself.
