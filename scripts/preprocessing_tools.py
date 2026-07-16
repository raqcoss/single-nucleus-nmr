import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scanpy as sc
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
import seaborn as sns
import os

def subsample_dataset(dataset_path: str, random_seed=71, fraction=0.05, min_per_cluster=34,
                      filter_by=None, add_info_to_stratify_by=None, class_col_name="cluster_id"):
    """
    Subsample an AnnData dataset in backed mode (don't modify backed object).
    Stratifies sampling by grouping cells and ensuring each group is represented.

    Parameters
    ----------
    dataset_path : str
        Path to the AnnData dataset in h5ad format (backed mode).
    random_seed : int, optional
        Random seed for reproducibility. Default is 71.
    fraction : float, optional
        Fraction of cells to sample from each stratification group. Default is 0.05 (5%).
    min_per_cluster : int, optional
        Minimum number of cells to sample from each stratification group. Default is 34.
    filter_by : dict, optional
        Dictionary of column names and regex patterns to filter cells before sampling. Default is None.
    add_info_to_stratify_by : tuple, optional
        Tuple of (column name, mapping dict) to add additional info for stratification. Default is None.
    class_col_name : str, optional
        Column name in obs that contains cluster labels (or the mapped stratification column after cluster_map).
        Default is "cluster_id".

    Returns
    -------
    AnnData
        Subsampled AnnData object.
    """
    dirname = os.path.dirname(dataset_path)
    basename = os.path.basename(dataset_path)
    subset_path = os.path.join(dirname, f"{basename.split(r'.h5')[0]}_subset{random_seed}.h5ad")

    if os.path.exists(subset_path):
        print("Loading existing subset...")
        human_data = sc.read(subset_path)
        print(f"Sampled dataset: {human_data.shape[0]} × {human_data.shape[1]}")
        return human_data

    print("Creating new subset...")
    adata_backed = sc.read_h5ad(dataset_path, backed="r")
    n_obs = adata_backed.n_obs

    # Build combined boolean mask (safe for multiple filters)
    if filter_by:
        mask = np.ones(n_obs, dtype=bool)
        for key, value in filter_by.items():
            if key not in adata_backed.obs.columns:
                raise NameError(f"Column '{key}' not found in obs")
            mask &= adata_backed.obs[key].str.contains(value, regex=True, na=False).values
    else:
        mask = np.ones(n_obs, dtype=bool)

    valid_positions = np.where(mask)[0]  # integer positions in the backed AnnData
    n_cells = len(valid_positions)
    target_size = int(n_cells * fraction)
    print(f"Original: {n_obs} cells, filtered: {n_cells} cells, target {target_size}")

    # Get cluster labels by positional indexing to avoid label/position ambiguity
    original_clusters = adata_backed.obs[class_col_name].iloc[valid_positions].reset_index(drop=True)

    # Map stratify values but fall back to original cluster if mapping misses
    if add_info_to_stratify_by:
        new_col, mapping_dict = add_info_to_stratify_by
        stratify_clusters = original_clusters.map(mapping_dict).fillna('Other').astype(str)
        print(f"Stratifying by mapped column with {stratify_clusters.nunique()} unique values")
    else:
        stratify_clusters = original_clusters.astype(str)

    np.random.seed(random_seed)
    sampled_indices = []
    summary = []

    # Use np.where on the values to guarantee integer positional arrays
    unique_groups = stratify_clusters.unique()
    for clust in unique_groups:
        pos_in_filtered = np.where(stratify_clusters.values == clust)[0]  # integer positions 0..len(valid_positions)-1
        if pos_in_filtered.size == 0:
            continue
        actual_indices = valid_positions[pos_in_filtered]  # convert to dataset row indices
        n_total = len(actual_indices)

        n_sample = int(n_total * fraction) + min_per_cluster
        n_sample = min(n_sample, n_total)

        chosen = np.random.choice(actual_indices, size=n_sample, replace=False)
        chosen_indexes = [int(x) for x in chosen]
        sampled_indices.extend(chosen_indexes)
        summary.append([clust, n_total, n_sample])

    sampled_indices = sorted(sampled_indices)
    if len(sampled_indices) == 0:
        raise ValueError("No cells sampled (check filters / fraction).")

    # Load only sampled cells into memory
    human_data = adata_backed[sampled_indices, :].to_memory()
    print(f"Final sampled dataset: {human_data.shape[0]} × {human_data.shape[1]}")

    df_summary = pd.DataFrame(summary, columns=["cluster", "total_cells", "sampled_cells"])
    df_summary["fraction_sampled"] = df_summary["sampled_cells"] / df_summary["total_cells"]
    print(df_summary.sort_values("total_cells").head(5))
    print("...")
    print(df_summary.sort_values("total_cells").tail(5))

    human_data.write_h5ad(subset_path)
    print(f"Sampled dataset saved to {subset_path}")
    return human_data

def describe_obs_by_cluster(adata, cluster_key, columns):
    """Computes descriptive statistics (mean, min, max, std) for numerical metadata grouped by clusters.

    This utility aggregates cell-level metrics (e.g., quality control metrics like total counts, 
    number of detected genes, or doublet scores) across defined cluster identities. It returns a 
    pivoted, multi-index summary table representing the distribution of those values within 
    each individual cluster.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix containing cell-level metadata in `adata.obs`.
    cluster_key : str
        The column name in `adata.obs` representing cell cluster assignments (e.g., "leiden", "subcluster").
    columns : list of str
        The continuous/numerical columns in `adata.obs` to summarize (e.g., `['total_counts', 'n_genes_by_counts']`).

    Returns
    -------
    stats : pandas.DataFrame
        A DataFrame with the cluster IDs as the row index. The columns are a pandas MultiIndex 
        where the first level corresponds to the summarized variables (`columns`) and the second 
        level corresponds to the computed statistics: `'mean'`, `'min'`, `'max'`, and `'std'`.

    Raises
    ------
    KeyError
        If `cluster_key` or any of the specified `columns` are missing from `adata.obs`.
    TypeError
        If the targeted summary columns contain non-numeric data types.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc3k()
    >>> sc.pp.calculate_qc_metrics(adata, inplace=True)
    >>> sc.tl.leiden(adata)
    >>> summary_df = describe_obs_by_cluster(
    ...     adata, 
    ...     cluster_key="leiden", 
    ...     columns=["total_counts", "n_genes_by_counts"]
    ... )
    >>> # Access the mean of total_counts for cluster '0'
    >>> summary_df.loc['0', ('total_counts', 'mean')]
    """
    df = adata.obs[[cluster_key] + columns].copy()
    df[cluster_key] = df[cluster_key].astype(str)

    grouped = df.groupby(cluster_key)   # pandas groupby

    stats = grouped.agg(["mean", "min", "max", "std"])
    return stats

def check_pc_correlation(adata, n_pcs=None, pc_col_name="X_pca"):
    """Checks statistical correlation between Principal Components (PCs) and cell-level QC metrics.

    This analysis helps identify potential technical or biological confounders (such as library size, 
    mitochondrial content, ribosomal content, or species) that might be driving major axes of 
    variation in the dataset. Strong correlations suggest that downstream regression or batch 
    correction may be necessary.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing computed dimensionality reduction coordinates in `adata.obsm` 
        and cell metadata in `adata.obs`.
    n_pcs : int, optional
        The number of top PCs to check. If None, checks all available PCs inside `adata.obsm[pc_col_name]`.
    pc_col_name : str, default 'X_pca'
        The key in `adata.obsm` representing the PCA coordinate array to evaluate.

    Returns
    -------
    corr : pandas.DataFrame
        A long-format DataFrame with correlation statistics containing the following columns:
        
        * `'PC'`: String identifier of the evaluated Principal Component (e.g., `"PC01"`).
        * `'covariate'`: Name of the QC metric or metadata column tested.
        * `'r'`: Spearman rank correlation coefficient (values range from -1.0 to 1.0).
        * `'p'`: Raw p-value calculated from the correlation test.
        * `'q'`: False Discovery Rate (FDR) adjusted p-value using the Benjamini-Hochberg method.

    Raises
    ------
    KeyError
        If `pc_col_name` is missing from `adata.obsm` or if expected QC columns/`species` 
        are not present in `adata.obs`.

    Notes
    -----
    - This function calculates **Spearman rank correlation** because it is non-parametric and robust 
      to outliers in both PC scores and technical metrics (such as read depth).
    - Categorical variables like `'species'` are automatically mapped to numerical codes to 
      facilitate correlation calculations.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc3k()
    >>> sc.pp.calculate_qc_metrics(adata, inplace=True)
    >>> sc.tl.pca(adata)
    >>> adata.obs['species'] = 'human'  # Mock covariate
    >>> corr_df = check_pc_correlation(adata, n_pcs=5)
    >>> # Filter for highly significant confounding variables (FDR < 0.05)
    >>> strong_confounders = corr_df[corr_df['q'] < 0.05]
    """
    from scipy.stats import spearmanr, pearsonr
    from statsmodels.stats.multitest import multipletests
    import matplotlib.pyplot as plt

    if n_pcs is None:
        n_pcs_check = adata.obsm[pc_col_name].shape[1]
    else: n_pcs_check = min(n_pcs, adata.obsm[pc_col_name].shape[1])
    pcs = pd.DataFrame(
        adata.obsm[pc_col_name][:, :n_pcs_check],
        index=adata.obs_names,
        columns=["PC"+str(i+1).zfill(2) for i in range(n_pcs_check)]
    )

    qc_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt", 'pct_counts_ribo']
    qc = adata.obs[qc_cols].astype(float)
    qc['species'] = adata.obs['species'].astype('category').cat.codes
    qc_cols.append('species')

    rows = []
    for pc in pcs.columns:
        for cov in qc_cols:
            # Spearman es más robusto; cambia a pearsonr si prefieres lineal
            r, p = spearmanr(pcs[pc].values, qc[cov].values, nan_policy="omit")
            rows.append({"PC": pc, "covariate": cov, "r": r, "p": p})

    corr = pd.DataFrame(rows)
    corr["q"] = multipletests(corr["p"], method="fdr_bh")[1]
    return corr