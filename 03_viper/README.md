# VIPER

| Script | Description |
|---|---|
| `01_network_to_regulon.R` | Reformats the ARACNe3 network to 3 columns (regulator, target, MI) and builds the regulon with `aracne2regulon` |
| `02_run_viper.R` | Runs VIPER on the metacell CPM matrix and adds the `VIPER` assay to the metacell object |
| `03_metacell_embeddings.R` | SCTransform on metacell counts (`SCT_mc`) → PCA (30) → UMAP (`umap_sct_mc`); scaled VIPER activity → PCA (30) → UMAP (`umap_viper_mc`) |
| `04_attach_chromvar.R` | Metacell chromVAR assay (`chromvar_mc`): mean of per-nucleus chromVAR z-scores over each metacell's member nuclei; also stores JASPAR motif names |
| `05_plot_metacell_umap.R` | **Figure 6a**: metacell UMAP colored by condition |

## VIPER parameters

| Parameter | Value |
|---|---|
| `method` | `rank` |
| `minsize` | 25 |
| `nes` | TRUE |
| `pleiotropy` | FALSE |

Outputs are written to `RESULTS_DIR/03_viper`. The object used by the next steps is `all_metacell_hdWGCNA_with_VIPER_sct_umap_viper_umap_chvar.rds`.
