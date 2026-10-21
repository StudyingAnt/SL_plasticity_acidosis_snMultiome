#!/usr/bin/env bash
# =============================================================================
# Submit the pipeline as a chain of SLURM jobs. Each job starts only after the
# job(s) it depends on finished successfully.
#
# Usage (from anywhere):
#   bash slurm/run_all.sh [start]
#
#   start = preprocessing  00-2 steps 01-07, then 01-05 (default)
#           chromvar       00-2 step 07, then 01-05 (06 output must exist)
#           metacells      01-05 only (06 and 07 outputs must exist)
#
# Steps that need more memory than the default partition's nodes provide
# (00-2 steps 05-07) use PARTITION_HIMEM if it is set, e.g.
#   PARTITION_HIMEM=long-96core-shared bash slurm/run_all.sh
# SNMULT_CONFIG_R / SNMULT_CONFIG_SH / JOB_TAG are passed on to submit.sh.
# =============================================================================
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "${REPO}"

START="${1:-preprocessing}"
case "${START}" in
  preprocessing|chromvar|metacells) ;;
  *) echo "Unknown start '${START}'. Use: preprocessing | chromvar | metacells" >&2; exit 1 ;;
esac

sub() { bash slurm/submit.sh "$@" | tail -1; }
himem() {
  if [ -n "${PARTITION_HIMEM:-}" ]; then PARTITION="${PARTITION_HIMEM}" sub "$@"; else sub "$@"; fi
}

summary=""
add() { summary="${summary}  $1=$2"$'\n'; }

# ---- 00-2: Seurat / Signac preprocessing ----
p6=""; p7=""
if [ "${START}" = "preprocessing" ]; then
  p1=$(sub 00-2_seurat_preprocessing/01_make_objects.R);                    add 00-2_01_make_objects "${p1}"
  p2=$(AFTER=${p1} sub 00-2_seurat_preprocessing/02_qc_objects.R);          add 00-2_02_qc_objects "${p2}"
  p3=$(AFTER=${p2} sub 00-2_seurat_preprocessing/03_call_peaks.R);          add 00-2_03_call_peaks "${p3}"
  p4=$(AFTER=${p3} sub 00-2_seurat_preprocessing/04_merge_objects.R);       add 00-2_04_merge_objects "${p4}"
  p5=$(AFTER=${p4} himem 00-2_seurat_preprocessing/05_transformations.R);   add 00-2_05_transformations "${p5}"
  p6=$(AFTER=${p5} himem 00-2_seurat_preprocessing/06_clustering_acdc.R);   add 00-2_06_clustering_acdc "${p6}"
  p7=$(AFTER=${p6} himem 00-2_seurat_preprocessing/07_run_chromvar.R);      add 00-2_07_run_chromvar "${p7}"
elif [ "${START}" = "chromvar" ]; then
  p7=$(himem 00-2_seurat_preprocessing/07_run_chromvar.R);                  add 00-2_07_run_chromvar "${p7}"
fi

# ---- 01-05: metacells -> Figure 6 ----
# metacells need 06 only; the metacell chromVAR step needs 07
m=$( if [ -n "${p6}" ]; then AFTER=${p6} sub 01_metacells/01_build_metacells.R; else sub 01_metacells/01_build_metacells.R; fi )
add 01_build_metacells "${m}"
r=$(AFTER=${m} sub 02_aracne3/01_prepare_regulators.R);                     add 02_01_prepare_regulators "${r}"
a=$(AFTER=${r} sub 02_aracne3/02_run_aracne3.sh);                           add 02_02_run_aracne3 "${a}"
n=$(AFTER=${a} sub 03_viper/01_network_to_regulon.R);                       add 03_01_network_to_regulon "${n}"
v=$(AFTER=${n} sub 03_viper/02_run_viper.R);                                add 03_02_run_viper "${v}"
e=$(AFTER=${v} sub 03_viper/03_metacell_embeddings.R);                      add 03_03_metacell_embeddings "${e}"
c=$(AFTER=${e}${p7:+:${p7}} sub 03_viper/04_attach_chromvar.R);             add 03_04_attach_chromvar "${c}"
u=$(AFTER=${e} sub 03_viper/05_plot_metacell_umap.R);                       add 03_05_plot_metacell_umap "${u}"
pr=$(AFTER=${c} sub 04_sphingolipid_axes/01_axis_pressure.R);               add 04_01_axis_pressure "${pr}"
d=$(AFTER=${pr} sub 04_sphingolipid_axes/02_dominant_axis.R);               add 04_02_dominant_axis "${d}"
dv=$(AFTER=${pr} sub 04_sphingolipid_axes/03_metabolic_diversity.R);        add 04_03_metabolic_diversity "${dv}"
ms=$(AFTER=${d} sub 05_axis_regulators/01_msviper_topbottom.R);             add 05_01_msviper_topbottom "${ms}"
l=$(AFTER=${ms} sub 05_axis_regulators/02_regulon_linkage_motif_support.R); add 05_02_regulon_linkage_motif_support "${l}"
g=$(AFTER=${l} sub 05_axis_regulators/03_plot_regulator_dotplots.R);        add 05_03_plot_regulator_dotplots "${g}"

echo "Submitted (start = ${START}):"
printf '%s' "${summary}"
