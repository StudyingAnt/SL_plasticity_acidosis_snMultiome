# Sphingolipid axes: regulatory pressure, dominant axis and diversity

## Ceramide clearance axes

The four axes were manually curated from the KEGG sphingolipid metabolism pathway (hsa00600). They are defined in `R/common.R`.

| Axis | Genes |
|---|---|
| S1P | ACER1, ACER2, ASAH1, ASAH2, SPHK1, SPHK2 |
| C1P | CERK |
| SM | SGMS1, SGMS2 |
| GSL | UGCG, UGT8, B4GALT4, B4GALT5 |

## Scripts

| Script | Description | Figure |
|---|---|---|
| `01_axis_pressure.R` | Regulatory pressure per axis and metacell | – |
| `02_dominant_axis.R` | Threshold diagnostics, dominant axis assignment, UMAP and proportion bar plots | 6b, 6c |
| `03_metabolic_diversity.R` | Per-metacell metabolic diversity and population diversity | 6d, 6e |

### Regulatory pressure

For an axis gene set *G*, target gene *g* and metacell *s*:

P<sub>gs</sub> = Σ<sub>r</sub> w<sub>rg</sub> · A<sub>rs</sub> / Σ<sub>r</sub> |w<sub>rg</sub>|,  with w<sub>rg</sub> = tfmode<sub>rg</sub> × likelihood<sub>rg</sub>

Here *A* is the VIPER activity and *w* comes from the regulon. The axis score is the mean of P<sub>gs</sub> over the genes in *G*, stored as `Pressure_s1p`, `Pressure_c1p`, `Pressure_sm` and `Pressure_gsl`.

### Dominant axis

- range4 = max − min of the four pressures
- margin4 = top-1 − top-2 pressure

Each metacell is assigned to the axis with the highest pressure. It is labeled `Ambiguous` when range4 < 1.793752 or margin4 < 0.1855446. These cut-offs are quantiles of the pressure distributions written by the threshold diagnostics (`threshold_diagnostics/thresholds_q10_q25.txt`).

### Metabolic diversity

- **Per metacell:** softmax (β = 1) of the four pressures, then Shannon entropy.
- **Per condition:** each axis is binned into K = 6 quantile bins, and the joint Shannon entropy of the 4D bin occupancy is computed.

Outputs are written to `RESULTS_DIR/04_sphingolipid_axes`.
