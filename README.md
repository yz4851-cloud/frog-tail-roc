# Frog tail ROC (v2)

Goal
Find the Regenerative Organizing Cell (ROC) in frog tail scRNA-seq and list its markers.

Data
Place `cleaned_processed_frogtail.h5ad` in `data/`. Optionally add `Supplementary_Table3.csv` with a column `gene`.

Run with Python only
1) `python run_all.py`
2) `python build_report.py`

Outputs
- Figures: `results/figures/` (`umap_clusters.png`, `roc_consensus.png`, `umap_leiden_labels.png`, `umap_skin_clusters.png`, `roc_violin.png`, `umap_prelim.png`)
- Tables: `results/` (`clustering_metrics.csv`, `roc_markers_logreg_wilcoxon.csv`, `overlap_with_supp_table3.csv`)
- Report: `REPORT.md`

Repro
Python 3.10 64-bit is recommended. Create a venv and install packages as in the project instructions, then run the two commands above.

Code Availability
https://github.com/yz4851-cloud/frog-tail-roc
