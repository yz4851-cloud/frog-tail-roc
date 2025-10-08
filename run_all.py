import os, sys, numpy as np, pandas as pd, scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score

def find_data():
    p1 = Path("data/cleaned_processed_frogtail.h5ad")
    p2 = Path(r"D:\cleaned_processed_frogtail.h5ad")
    if p1.exists(): return p1
    if p2.exists(): return p2
    raise FileNotFoundError("cleaned_processed_frogtail.h5ad not found in data/ or D:\")

def ensure_dirs():
    Path("results/figures").mkdir(parents=True, exist_ok=True)

def setup(adata):
    adata.obs.head(50).to_csv("results/obs_preview.tsv", sep="\t")
    adata.var.head(50).to_csv("results/var_preview.tsv", sep="\t")

def preprocess(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    if 'mt' not in adata.var.columns:
        adata.var['mt'] = adata.var_names.str.upper().str.startswith(('MT-','MT.'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat')
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    keys = [c for c in ['tissue','batch','condition','time','stage','cell_type'] if c in adata.obs.columns]
    sc.pl.umap(adata, color=keys[:1] if keys else None, show=False)
    plt.savefig("results/figures/umap_prelim.png", bbox_inches="tight")
    plt.close()
    adata.write_h5ad("results/preprocessed.h5ad")
    return adata

def _labels_to_int(a):
    if pd.api.types.is_categorical_dtype(a):
        return a.cat.codes.to_numpy()
    return pd.factorize(a)[0]

def clustering_and_metrics():
    adata = sc.read_h5ad("results/preprocessed.h5ad")
    sc.tl.louvain(adata, resolution=0.6, key_added='louvain_r06')
    sc.tl.leiden(adata, resolution=0.6, key_added='leiden_r06', flavor='igraph', n_iterations=2, directed=False)
    X = adata.obsm['X_pca'][:, :30]
    k = len(adata.obs['leiden_r06'].unique())
    km = KMeans(n_clusters=max(k,2), n_init=25, random_state=0)
    adata.obs['kmeans_k_leiden'] = km.fit_predict(X).astype(str)
    sc.pl.umap(adata, color=['louvain_r06','leiden_r06','kmeans_k_leiden'], ncols=3, show=False)
    plt.savefig("results/figures/umap_clusters.png", bbox_inches="tight")
    plt.close()
    ari_ll = adjusted_rand_score(adata.obs['louvain_r06'], adata.obs['leiden_r06'])
    y_louv = _labels_to_int(adata.obs['louvain_r06'])
    y_leid = _labels_to_int(adata.obs['leiden_r06'])
    y_km = _labels_to_int(adata.obs['kmeans_k_leiden'])
    sil_louv = silhouette_score(X, y_louv)
    sil_leid = silhouette_score(X, y_leid)
    sil_km = silhouette_score(X, y_km)
    pd.DataFrame({
        'metric':['ARI(Louvain,Leiden)','Silhouette(Louvain)','Silhouette(Leiden)','Silhouette(k-means)'],
        'value':[ari_ll,sil_louv,sil_leid,sil_km]
    }).to_csv("results/clustering_metrics.csv", index=False)
    sc.pl.umap(adata, color=['leiden_r06'], legend_loc='on data', show=False)
    plt.savefig("results/figures/umap_leiden_labels.png", bbox_inches="tight")
    plt.close()
    adata.write_h5ad("results/preprocessed.h5ad")

def top_markers(adata_obj, key, group, n=100):
    names = pd.DataFrame(adata_obj.uns[key]['names'])[group][:n].tolist()
    scores = pd.DataFrame(adata_obj.uns[key]['scores'])[group][:n].tolist()
    return pd.DataFrame({'gene': names, 'score': scores})

def markers_and_overlap():
    adata = sc.read_h5ad("results/preprocessed.h5ad")
    if 'tissue' in adata.obs.columns and 'Skin' in set(adata.obs['tissue']):
        skin = adata[adata.obs['tissue'] == 'Skin'].copy()
    else:
        skin = adata.copy()
    n_pcs = min(30, skin.obsm['X_pca'].shape[1]) if 'X_pca' in skin.obsm else 30
    if 'X_pca' not in skin.obsm:
        sc.pp.scale(skin, max_value=10)
        sc.tl.pca(skin, svd_solver='arpack')
    sc.pp.neighbors(skin, n_neighbors=15, n_pcs=n_pcs)
    sc.tl.umap(skin)
    sc.tl.leiden(skin, resolution=0.8, key_added='leiden_skin', flavor='igraph', n_iterations=2, directed=False)
    sc.pl.umap(skin, color=['leiden_skin'], legend_loc='on data', show=False)
    plt.savefig("results/figures/umap_skin_clusters.png", bbox_inches="tight")
    plt.close()
    sc.tl.rank_genes_groups(skin, 'leiden_skin', method='logreg', key_added='rank_logreg')
    sc.tl.rank_genes_groups(skin, 'leiden_skin', method='wilcoxon', key_added='rank_wilcoxon')
    panel = ['fgf8','wnt5b','wnt3a','sox2','sox10','ctgfa','mki67']
    genes_exist = [g for g in panel if g in skin.var_names]
    if genes_exist:
        sc.pp.scale(skin, max_value=10, copy=False)
        s = np.asarray(skin[:, genes_exist].X.mean(axis=1)).ravel()
        skin.obs['roc_score'] = s
        roc_cluster = skin.obs.groupby('leiden_skin')['roc_score'].mean().idxmax()
    else:
        roc_cluster = skin.obs['leiden_skin'].value_counts().idxmax()
    logreg_df = top_markers(skin, 'rank_logreg', roc_cluster, n=100)
    wilc_df = top_markers(skin, 'rank_wilcoxon', roc_cluster, n=100)
    merged = logreg_df.merge(wilc_df, on='gene', how='outer', suffixes=('_logreg','_wilcoxon'))
    merged.to_csv("results/roc_markers_logreg_wilcoxon.csv", index=False)
    consensus = merged.dropna(subset=['score_logreg','score_wilcoxon']).sort_values('score_logreg', ascending=False).head(12)['gene'].tolist()
    if len(consensus) > 0:
        sc.pl.dotplot(skin, consensus, groupby='leiden_skin', show=False)
        plt.savefig("results/figures/roc_consensus.png", bbox_inches="tight")
        plt.close()
        sc.pl.violin(skin, consensus[:6], groupby='leiden_skin', stripplot=False, multi_panel=True, show=False)
        plt.savefig("results/figures/roc_violin.png", bbox_inches="tight")
        plt.close()
    sup = Path("data/Supplementary_Table3.csv")
    if sup.exists():
        sup3 = pd.read_csv(sup)
    else:
        sup3 = pd.DataFrame(columns=['gene'])
    roc_set = set(merged['gene'].dropna().tolist())
    sup_set = set(sup3['gene'].astype(str).str.strip())
    overlap = sorted(roc_set.intersection(sup_set))
    pd.Series(overlap, name='overlap_gene').to_csv("results/overlap_with_supp_table3.csv", index=False)
    skin.write_h5ad("results/skin_subset.h5ad")

if __name__ == "__main__":
    ensure_dirs()
    data_path = find_data()
    adata = sc.read_h5ad(data_path)
    setup(adata)
    adata = preprocess(adata)
    clustering_and_metrics()
    markers_and_overlap()
    print("done")
