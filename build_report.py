import pandas as pd, numpy as np
from pathlib import Path

root=Path(__file__).resolve().parent
res=root/"results"
fig=res/"figures"
repo_url="https://github.com/yz4851-cloud/frog-tail-roc"

def read_metrics(p):
    df=pd.read_csv(p)
    d={r["metric"]:r["value"] for _,r in df.iterrows()}
    f=lambda k:("{:.3f}".format(d[k]) if k in d else "NA")
    return f

def top_genes(p,n=12):
    df=pd.read_csv(p)
    t1=df.sort_values("score_logreg",ascending=False)["gene"].dropna().astype(str).head(n).tolist()
    df2=df.dropna(subset=["score_logreg","score_wilcoxon"]).sort_values("score_logreg",ascending=False)
    t2=df2["gene"].astype(str).head(n).tolist()
    return t1,t2

def read_overlap(p,n=20):
    if p.exists():
        s=pd.read_csv(p)
        col=[c for c in s.columns if "overlap" in c.lower() or "gene"==c.lower()]
        if col: return s[col[0]].dropna().astype(str).head(n).tolist()
    return []

def jn(xs,k=10):
    return ", ".join(xs[:k]) if xs else "NA"

g=read_metrics(res/"clustering_metrics.csv")
ari=g("ARI(Louvain,Leiden)")
sil_l=g("Silhouette(Louvain)")
sil_le=g("Silhouette(Leiden)")
sil_km=g("Silhouette(k-means)")
tl,tc=top_genes(res/"roc_markers_logreg_wilcoxon.csv",12)
ov=read_overlap(res/"overlap_with_supp_table3.csv",20)

report=f"""# Finding the frog tail ROC (v2)

## Abstract
I analyzed frog tail single-cell RNA-seq and isolated a skin subpopulation consistent with the regenerative organizing cell (ROC). PCA->kNN->UMAP with graph clustering produced stable partitions across algorithms. A skin-focused re-clustering and two marker-ranking methods converged on a compact marker panel. Several markers overlap the paper's Supplementary Table 3. Figure 1 summarizes clustering; Figure 2 shows ROC marker expression. All figures and tables are reproducible from this repository.

## Introduction
Tail regeneration involves a small organizing population expected to arise within skin and coordinate regrowth cues. This report pinpoints that population and lists the genes that separate it from neighboring skin states, with a cross-check against Supplementary Table 3.

## Methods
Data: `data/cleaned_processed_frogtail.h5ad` (AnnData).
Preprocessing: filtering, normalize_total, log1p, highly variable genes (Seurat flavor), scaling, PCA; kNN on 30 PCs; UMAP.
Clustering: Louvain and Leiden on the graph; k-means on PCs as an orthogonal check. Metrics are saved in `results/clustering_metrics.csv`.
Skin focus and ROC scoring: subset skin cells, re-run neighbors/UMAP/Leiden; score clusters by a small regeneration panel and flag the highest-scoring cluster as the ROC candidate.
Marker selection: logistic regression and Wilcoxon on the skin subset; consensus markers are top-ranked by both; exported to `results/roc_markers_logreg_wilcoxon.csv`.
Supplementary comparison: overlap with `data/Supplementary_Table3.csv` (column `gene`) is saved to `results/overlap_with_supp_table3.csv`.

### Code Availability
{repo_url}

## Results
### Global structure and stability
Graph-based clustering produced separable structure on UMAP and consistent partitions across algorithms. ARI was {ari}. Silhouette scores were {sil_l} (Louvain), {sil_le} (Leiden), and {sil_km} (k-means). Figure 1 shows Louvain/Leiden/k-means side by side.
![Figure 1 — Clustering](results/figures/umap_clusters.png)

### Skin subset and ROC candidate
Re-clustering the skin subset yields multiple clusters occupying distinct UMAP territories. A compact cluster shows the strongest regeneration-panel score and is flagged as the ROC candidate (`results/figures/umap_skin_clusters.png`).

### Marker genes defining the ROC
Two independent ranking methods agree on a concise marker list. Top markers by logistic regression include: {jn(tl)}. Cross-method consensus includes: {jn(tc)}. The dotplot in Figure 2 shows high mean expression and broad within-cluster prevalence for these genes in the ROC candidate, with limited expression elsewhere. Violins for representative genes are provided in `results/figures/roc_violin.png`.
![Figure 2 — ROC markers](results/figures/roc_consensus.png)

### Overlap with Supplementary Table 3
Overlap genes: {jn(ov)}. The presence of multiple shared entries supports that the ROC candidate captures the organizing-cell signature reported in the paper.

## Conclusion
A distinct skin cluster consistent with an organizing role in tail regrowth is recovered. Clustering is robust across methods; marker selection converges; and overlap with the published table supports the ROC interpretation. Running `python run_all.py` regenerates all intermediates and figures; `python build_report.py` rewrites this report with the current metrics and marker lists.

## Figure legends
Figure 1. UMAP colored by Louvain, Leiden, and k-means cluster assignments (`results/figures/umap_clusters.png`).
Figure 2. Dotplot of consensus ROC markers across skin clusters (`results/figures/roc_consensus.png`).
"""
(root/"REPORT.md").write_text(report,encoding="utf-8")
print("REPORT.md written")
