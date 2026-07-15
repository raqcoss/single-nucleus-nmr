import os
import scanpy as sc
import harmonypy as hm 
import pandas as pd
import numpy as np
import scipy.sparse as sp
from sklearn.metrics import silhouette_score, silhouette_samples
import matplotlib.pyplot as plt
import warnings
import traceback

def compute_silhouettes(adata, labels_key="cluster", embedding_key="X_pca"):
    """
    Computes global and sample-wise silhouette scores for a given clustering.

    Parameters
    ----------
    adata : AnnData
        Annotated data object.
    labels_key : str, default="cluster"
        The column in `adata.obs` containing the cluster labels.
    embedding_key : str, default="X_pca"
        The key in `adata.obsm` representation to use for calculating distances.

    Returns
    -------
    avg : float
        Average silhouette score across all samples.
    sample_vals : numpy.ndarray
        Silhouette score for each individual cell.
    """
    labels = adata.obs[labels_key].astype(str)
    if labels.nunique() < 2:
        raise ValueError("Need at least 2 clusters to compute silhouette scores")
    if embedding_key not in adata.obsm:
        raise KeyError(f"Embedding '{embedding_key}' not found in adata.obsm")
    embedding = adata.obsm[embedding_key]
    avg = silhouette_score(embedding, labels)
    sample_vals = silhouette_samples(embedding, labels)
    return avg, sample_vals


def plot_silhouette(
    adata, 
    sample_sil_values, 
    labels_key="cluster", 
    avg_sil_score=None, 
    title=None, 
    save_path=None,
    show_plot=True
):
    """
    Generates and optionally displays/saves a silhouette plot for cluster evaluations.

    Parameters
    ----------
    adata : AnnData
        Annotated data object.
    sample_sil_values : numpy.ndarray
        Array containing pre-computed silhouette values for each cell.
    labels_key : str, default="cluster"
        The column in `adata.obs` containing cluster labels.
    avg_sil_score : float, optional
        The pre-calculated average silhouette score to display as a red vertical dashed line.
    title : str, optional
        Title of the generated plot.
    save_path : str, optional
        File path where the generated plot image will be saved.
    show_plot : bool, default=True
        Whether to display the generated plot dynamically in the active session.
    """
    labels = adata.obs[labels_key].astype(str)
    unique_labels = sorted(labels.unique(), key=lambda x: str(x))
    y_lower = 10
    fig, ax1 = plt.subplots(1, 1, figsize=(8, 4))
    for i, label in enumerate(unique_labels):
        label_sil_vals = np.sort(sample_sil_values[labels == label])
        size_cluster = label_sil_vals.shape[0]
        if size_cluster == 0:
            continue
        y_upper = y_lower + size_cluster
        color = plt.cm.nipy_spectral(float(i) / max(len(unique_labels), 1))
        ax1.fill_betweenx(np.arange(y_lower, y_upper), 0, label_sil_vals, facecolor=color, alpha=0.7)
        ax1.text(-0.05, y_lower + 0.5 * size_cluster, str(label), fontsize=8)
        y_lower = y_upper + 10
    title = title or f"Silhouette plot ({labels_key})"
    ax1.set_title(title)
    ax1.set_xlabel("Silhouette coefficient values")
    ax1.set_ylabel("Cluster")
    if avg_sil_score is not None:
        ax1.axvline(x=avg_sil_score, color="red", linestyle="--", label=f"Avg = {avg_sil_score:.3f}")
    ax1.set_xlim([-1, 1])
    ax1.set_xticks(np.arange(-1.0, 1.1, 0.2))
    ax1.legend(loc="lower right")
    plt.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=200)
    
    if show_plot:
        plt.show()
    else:
        plt.close(fig)
    return fig


def process_adata(
    adata, 
    hvgs_n=3000, 
    n_pcs=15, 
    k=25, 
    resolution=0.1, 
    batch_key="sample", 
    cluster_key="supercluster", 
    parent_cluster_key=None,
    seed=42
):
    """
    Executes a standard single-cell processing pipeline.

    Normalization, log transformation, highly variable gene selection, PCA reduction, 
    Harmony batch correction, neighbor construction, Leiden clustering, and UMAP embeddings are 
    applied in sequence.

    Parameters
    ----------
    adata : AnnData
        The input annotated single-cell data matrix to modify in-place.
    hvgs_n : int, default=3000
        Number of highly variable genes to identify.
    n_pcs : int, default=15
        Number of Principal Components to compute.
    k : int, default=25
        The neighborhood size parameter for nearest neighbors construction.
    resolution : float, default=0.1
        The resolution parameter used for Leiden clustering.
    batch_key : str, default="sample"
        The variable in `adata.obs` defining batch labels for Harmony.
    cluster_key : str, default="supercluster"
        The key where final clustering labels will be stored in `adata.obs`.
    parent_cluster_key : str, optional
        If provided, prefixes the output cluster IDs with parent labels.
    """
    # Store raw counts in a layer for later use in DE analysis
    adata.layers["counts"] = adata.X.copy()
    
    # Normalize & log transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Defining embeddings
    sc.pp.highly_variable_genes(adata, n_top_genes=hvgs_n, subset=False, flavor="seurat_v3", batch_key=batch_key)

    # Standard normalization and log transformation
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, svd_solver='arpack', n_comps=n_pcs, use_highly_variable=True)
    
    harmony_out = hm.run_harmony(
        adata.obsm['X_pca'], 
        adata.obs, 
        batch_key, 
        max_iter_harmony=20, 
        max_iter_kmeans=3, 
        theta=1.5, 
        lamb=1, 
        sigma=0.1, 
        verbose=True, 
        random_state=42 
    )
    adata.obsm['X_pca_harmony'] = harmony_out.Z_corr 
    sc.pp.neighbors(adata, n_neighbors=k, use_rep="X_pca_harmony", random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key, random_state=seed)
    
    if parent_cluster_key is not None:
        adata.obs[cluster_key] = adata.obs[parent_cluster_key].astype(str) + "/" + adata.obs[cluster_key].astype(str)
    sc.tl.umap(adata, random_state=seed)


def run_recursive_hpo_subclustering_with_layers(
    adata,
    batch_key="replicate",
    max_depth=3,
    sil_threshold=0.02,
    min_cells=50,
    show_plots=True,
    seed=42
):
    """
    Performs recursive unsupervised hyperparameter-optimized (HPO) subclustering.

    This function dynamically sweeps Leiden resolutions and evaluates local split quality 
    at each hierarchical node using silhouette scores. Valid splits are accepted and 
    recursively evaluated down to the configured max_depth threshold.

    Parameters
    ----------
    adata : AnnData
        Annotated data object. Must contain 'counts' in `adata.layers`.
    batch_key : str, default="replicate"
        Metadata column defining batches used to perform local Harmony integrations.
    max_depth : int, default=3
        The maximum recursion depth allowed.
    sil_threshold : float, default=0.02
        The minimum average silhouette score delta required to accept a sub-split.
    min_cells : int, default=50
        The minimum size threshold. Nodes smaller than this are not split further.
    show_plots : bool, default=True
        Controls whether intermediate silhouette and UMAP diagnostic figures are displayed 
        via matplotlib's standard popups. Saving of images to disk is preserved in either state.
    """
    os.makedirs('figures', exist_ok=True)
    region_name = globals().get('region', 'analysis')
    tracking_df = pd.DataFrame(
        index=adata.obs_names,
        columns=["supercluster", "cluster", "subcluster"]
    )
    print("--- Computing Superclusters ---")
    sc.tl.leiden(
        adata,
        resolution=0.1,
        key_added="temp_supercluster"
    )
    tracking_df["supercluster"] = adata.obs["temp_supercluster"].astype(str)
    tracking_df["cluster"] = tracking_df["supercluster"]
    tracking_df["subcluster"] = tracking_df["supercluster"]

    def sanitize_node_name(name):
        return name.replace('/', '_').replace(' ', '_')

    def evaluate_and_split_node(
        cell_names,
        parent_path,
        current_depth
    ):
        if current_depth >= max_depth:
            return
        if len(cell_names) < min_cells:
            print(f"[{parent_path}] ⏭ Too few cells ({len(cell_names)})")
            return
        print(f"[{parent_path}] Sweeping HPO Depth {current_depth + 1}")
        try:
            adata_sub = adata[cell_names].copy()
            if adata_sub.n_obs < min_cells:
                return
            if batch_key not in adata_sub.obs:
                raise ValueError(f"{batch_key} not found in obs")
            batch = adata_sub.obs[batch_key].copy()
            if pd.api.types.is_numeric_dtype(batch):
                batch = batch.astype(str)
            batch = batch.astype("category").cat.remove_unused_categories()
            adata_sub.obs[batch_key] = batch
            if "counts" not in adata_sub.layers:
                raise ValueError("counts layer missing")
            raw = adata_sub.layers["counts"]
            if sp.issparse(raw):
                gene_var = np.asarray(raw.power(2).mean(axis=0) - np.square(raw.mean(axis=0))).ravel()
            else:
                gene_var = raw.var(axis=0)
            keep = gene_var > 1e-4
            if keep.sum() < 20:
                print(f"[{parent_path}] ⏭ No variable genes")
                return
            adata_sub = adata_sub[:, keep].copy()
            n_genes = min(3000, adata_sub.n_vars, max(50, adata_sub.n_obs // 2))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sc.pp.highly_variable_genes(
                    adata_sub,
                    flavor="seurat_v3",
                    layer="counts",
                    n_top_genes=n_genes,
                    subset=True
                )
            if adata_sub.n_vars < 15:
                print(f"[{parent_path}] ⏭ Insufficient HVGs")
                return
            n_comps = min(50, adata_sub.n_obs - 1, adata_sub.n_vars - 1)
            sc.pp.pca(adata_sub, n_comps=n_comps)
            neighbors_rep = "X_pca"
            batch = adata_sub.obs[batch_key]
            batch_sizes = batch.value_counts()
            valid_batches = batch_sizes[batch_sizes >= 2].index
            if len(valid_batches) >= 2:
                keep = batch.isin(valid_batches).values
                adata_sub = adata_sub[keep].copy()
                adata_sub.obs[batch_key] = adata_sub.obs[batch_key].astype(str).astype("category")
                try:
                    harmony_out = hm.run_harmony(
                        adata_sub.obsm['X_pca'], 
                        adata_sub.obs, 
                        batch_key, 
                        max_iter_harmony=20, 
                        max_iter_kmeans=3, 
                        theta=1.5, 
                        lamb=1, 
                        sigma=0.1, 
                        verbose=True,
                        random_state=seed
                    )
                    adata_sub.obsm['X_pca_harmony'] = harmony_out.Z_corr 
                    neighbors_rep = "X_pca_harmony"
                    print(f"[{parent_path}] Harmony applied")
                except Exception as e:
                    print(f"[{parent_path}] Harmony skipped ({e})")
            else:
                print(f"[{parent_path}] Single batch → no Harmony")
            sc.pp.neighbors(
                adata_sub,
                use_rep=neighbors_rep,
                metric="cosine",
                n_neighbors=max(15, int(np.log(adata_sub.n_obs)*10))
            )
            best_sil = -1
            best_res = None
            best_labels = None
            best_key = None
            for r in [0.1, 0.3, 0.5]:
                key = f"tmp_{r}"
                sc.tl.leiden(adata_sub, resolution=r, key_added=key)
                labels = adata_sub.obs[key].astype(str)
                n_clusters = labels.nunique()
                if n_clusters < 2 or n_clusters >= adata_sub.n_obs:
                    continue
                try:
                    sil = silhouette_score(adata_sub.obsm[neighbors_rep], labels)
                    if sil > best_sil:
                        best_sil = sil
                        best_res = r
                        best_labels = labels.copy()
                        best_key = key
                except Exception:
                    continue
        except Exception:
            traceback.print_exc()
            return
        if best_res is None or best_sil < sil_threshold:
            print(f"[{parent_path}] Stable node")
            return
        print(f"[{parent_path}] Split accepted → res={best_res} sil={best_sil:.3f}")
        node_name = sanitize_node_name(parent_path)
        label_col = f"{node_name}_labels"
        sil_col = f"{node_name}_silhouette"
        adata_sub.obs[label_col] = best_labels
        adata_sub.obs[sil_col] = silhouette_samples(adata_sub.obsm[neighbors_rep], best_labels)
        if "X_umap" not in adata_sub.obsm:
            sc.tl.umap(adata_sub, random_state=seed)
            
        plot_silhouette(
            adata_sub,
            adata_sub.obs[sil_col].values,
            labels_key=label_col,
            avg_sil_score=best_sil,
            title=f"{parent_path} depth {current_depth} silhouette",
            save_path=f"figures/{region_name}_{node_name}_depth{current_depth}_silhouette.png",
            show_plot=show_plots
        )
        
        # Determine whether to block cell execution output
        show_umap = None if show_plots else False
        sc.pl.umap(
            adata_sub,
            color=[label_col, sil_col],
            cmap="coolwarm",
            title=[f"{parent_path} labels", f"{parent_path} silhouette"],
            legend_loc='on data',
            save=f"_{region_name}_{node_name}_depth{current_depth}.png",
            show=show_umap
        )
        target = "cluster" if current_depth == 1 else "subcluster"
        for cid in best_labels.unique():
            child = best_labels.index[best_labels == cid]
            child_path = f"{parent_path}/{cid}"
            tracking_df.loc[child, target] = child_path
            if target == "cluster":
                tracking_df.loc[child, "subcluster"] = child_path
            evaluate_and_split_node(child, child_path, current_depth + 1)
            
    for super_id in tracking_df["supercluster"].unique():
        print("" + "=" * 60)
        print(f"Processing Supercluster {super_id}")
        print("=" * 60)
        cells = tracking_df.index[tracking_df["supercluster"] == super_id]
        evaluate_and_split_node(cells, str(super_id), 1)
        
    adata.obs["supercluster"] = tracking_df["supercluster"].astype("category")
    adata.obs["cluster"] = tracking_df["cluster"].astype("category")
    adata.obs["subcluster"] = tracking_df["subcluster"].astype("category")
    adata.obs.drop(columns=["temp_supercluster"], errors="ignore", inplace=True)
    print("🎉 Recursive Custering with HPO completed")


def run_sc_annotation_pipeline(
    adata, 
    batch_key="sample",
    replicate_key="replicate",
    hvgs_n=3000, 
    n_pcs=20, 
    supercluster_k=25, 
    supercluster_resolution=0.1,
    max_depth=3,
    sil_threshold=0.05,
    min_cells=50,
    show_plots=False, 
    random_state=42
):
    """
    An integrated orchestration function to run the complete clustering and subclustering pipeline.

    This function acts as a wrapper. It first executes standard global preprocessing 
    and superclustering via `process_adata`, and then initiates the recursive 
    hyperparameter-optimized (HPO) subclustering down into subclusters.

    Parameters
    ----------
    adata : AnnData
        The input annotated single-cell data matrix to process in-place.
    batch_key : str, default="sample"
        The categorical key in `adata.obs` representing the batch used for global Harmony integration.
    replicate_key : str, default="replicate"
        The categorical key in `adata.obs` representing replicates for local subclustering integration.
    hvgs_n : int, default=3000
        Number of highly variable genes to use for the global integration steps.
    n_pcs : int, default=20
        Number of global Principal Components to retain.
    supercluster_k : int, default=25
        The neighborhood size for global superclustering.
    supercluster_resolution : float, default=0.1
        The clustering resolution for global superclustering.
    max_depth : int, default=3
        The maximum permitted recursion depth for subclustering.
    sil_threshold : float, default=0.05
        The threshold change in silhouette score to validate a local node split.
    min_cells : int, default=50
        Smallest size group of cells considered for a split evaluation.
    show_plots : bool, default=False
        Whether to show generated matplotlib plotting layouts. Saves plot figures regardless.

    Returns
    -------
    adata : AnnData
        The fully-annotated object modified in place with added observations: 
        `supercluster`, `cluster`, and `subcluster`.
    """
    print("STEP 1: Processing Global Embeddings & Superclusters")
    process_adata(
        adata, 
        hvgs_n=hvgs_n, 
        n_pcs=n_pcs, 
        k=supercluster_k, 
        resolution=supercluster_resolution, 
        batch_key=batch_key, 
        cluster_key="supercluster", 
        seed=random_state
    )
    
    print("\nSTEP 2: Running Recursive HPO Subclustering")
    run_recursive_hpo_subclustering_with_layers(
        adata, 
        batch_key=replicate_key, 
        max_depth=max_depth, 
        sil_threshold=sil_threshold,
        min_cells=min_cells,
        show_plots=show_plots,
        seed=random_state
    )
    
    return adata