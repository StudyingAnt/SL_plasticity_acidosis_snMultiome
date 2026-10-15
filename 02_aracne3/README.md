# ARACNe3 network inference

## Regulator list

The regulators are human transcription factors from the Human TFs database v1.01 (Lambert *et al.*, 2018). Download `TF_names_v_1.01.txt` from https://humantfs.ccbr.utoronto.ca/download.php and set `HUMAN_TF_FILE` in `config/config.R`.

`01_prepare_regulators.R` keeps only the TFs present in the metacell CPM matrix and writes `regulators_in_matrix.txt`.

## ARACNe3

ARACNe3 (https://github.com/califano-lab/ARACNe3) is run on the metacell CPM matrix by `02_run_aracne3.sh`. It needs `config/config.sh`.

| Parameter | Value |
|---|---|
| Bootstrap subnetworks (`-x`) | 30 |
| `--alpha` | 0.05 |
| `--threads` | 8 |
| `--seed` | 1234 |
| `--runid` | `mcf7_metacell` |

```bash
Rscript 02_aracne3/01_prepare_regulators.R
bash 02_aracne3/02_run_aracne3.sh
```

Main output: `RESULTS_DIR/02_aracne3/consolidated-net_mcf7_metacell.tsv`
