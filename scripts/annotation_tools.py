import gseapy as gp
from gseapy.plot import gseaplot, gseaplot2
import pandas as pd
import numpy as np
def annotate_with_celltypist(adata, model_name: str, subcluster_col: str = "leiden_4.0", min_prop: float = 0.2, plot_umap: bool = True):
    """
    Annotate single-cell data using a CellTypist model and plot the results on a UMAP.

    This function attempts to run CellTypist annotation on the input AnnData object. 
    If the specified model is not found locally, it will attempt to download it 
    automatically. If successful, it adds the raw and majority-voting predictions 
    to the AnnData's `.obs` and generates a UMAP plot.

    Parameters
    ----------
    vdata : AnnData
        The AnnData object containing the single-cell dataset to annotate.
    model_name : str
        The name of the CellTypist model (e.g., 'Immune_All_LowResolution') 
        without the '.pkl' extension.

    Returns
    -------
    AnnData
        The annotated AnnData object with celltypist predictions added to `.obs`.
    """
    import celltypist
    from celltypist import models
    import scanpy as sc
    model_file = f"{model_name}.pkl"
    
    try:
        # Attempt annotation with the local model
        predictions = celltypist.annotate(
            adata, 
            model=model_file, 
            majority_voting=True, 
            over_clustering=subcluster_col, 
            min_prop=min_prop
        )
    except Exception: 
        print(f"Model '{model_file}' not found. Trying to download it...")
        try: 
            # Attempt to download the model
            models.download_models(model=model_file)
            
            # Retry annotation after downloading
            predictions = celltypist.annotate(
                adata, 
                model=model_file, 
                majority_voting=True, 
                over_clustering='leiden_4.0', 
                min_prop=0.2
            )
        except Exception as e: 
            print(f"Model '{model_name}' is not available to download. Skipping. Error: {e}")
            return adata

    # Map predicted labels to the AnnData metadata
    adata.obs[f'celltypist_{model_name}'] = adata.obs_names.map(predictions.predicted_labels.predicted_labels)
    adata.obs[f'celltypist_{model_name}_mv'] = adata.obs_names.map(predictions.predicted_labels.majority_voting)
    
    # Plot the results
    if plot_umap:
        sc.pl.umap(adata, color=[f'celltypist_{model_name}', f'celltypist_{model_name}_mv'])
    
    return adata


def run_gsea_on_group(
    clus_markers,
    gene_sets="Allen_Brain_Atlas_10x_scRNA_2021",
    min_perm=1000,
    seed=6,
):
    """
    Performs Gene Set Enrichment Analysis (GSEA) Prerank on a set of cluster markers.

    This function calculates a custom ranking metric combining the sign of log-fold 
    change and the -log10 of the adjusted p-value. It cleans the ranked list 
    (removing duplicate genes, NaNs, and infinities), runs `gseapy.prerank`, and 
    safely attempts to plot the top enriched pathway using either `gseaplot` or 
    `gseaplot2` as a fallback.

    Parameters
    ----------
    clus_markers : pandas.DataFrame
        DataFrame of differential expression markers. Must contain columns: 
        'names' (gene symbols), 'pvals_adj' (adjusted p-values), and 
        'logfoldchanges' (log fold-changes).
    gene_sets : str or list, default="Allen_Brain_Atlas_10x_scRNA_2021"
        The library name or list of gene sets to query via GSEAPy/Enrichr.
    min_perm : int, default=1000
        Number of permutations to run for GSEA.
    seed : int, default=6
        Random seed for reproducibility during permutations.

    Returns
    -------
    pre_res : gseapy.Prerank
        The raw fitted GSEAPy Prerank object containing complete results.
    out_df : pandas.DataFrame
        A cleaned, sorted summary DataFrame of terms, False Discovery Rates (fdr), 
        Enrichment Scores (es), and Normalized Enrichment Scores (nes). Empty if 
        no results are found.
    term_to_graph : str or None
        The name of the top-ranked enriched term (by FDR) selected for plotting. 
        None if no results are found.
    gsea_plot : matplotlib.figure.Figure or None
        The generated enrichment plot object. None if no results are found.

    Notes
    -----
    The custom ranking metric is defined as:
    $$Rank = -\\log_{10}(pval\\_adj) \\times logfoldchanges$$
    To avoid math errors from perfectly zero p-values (common in tools like Scanpy), 
    adjusted p-values are clipped at a lower bound of $10^{-300}$ ($1e-300$) 
    prior to log transformation.
    """
    # Defensive copy
    clus_markers = clus_markers.copy()

    # Avoid log(0); Scanpy can produce pvals_adj == 0
    eps = 1e-300
    clus_markers["Rank"] = (
        -np.log10(clus_markers["pvals_adj"].clip(lower=eps))
        * clus_markers["logfoldchanges"]
    )

    # Remove invalid entries
    clus_markers = (
        clus_markers
        .replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["names", "Rank"])
        .drop_duplicates(subset="names")
        .sort_values("Rank", ascending=False)
        .reset_index(drop=True)
    )

    ranking = clus_markers[["names", "Rank"]]

    pre_res = gp.prerank(
        rnk=ranking,
        gene_sets=gene_sets,
        seed=seed,
        permutation_num=min_perm,
        verbose=False,
    )

    if not pre_res.results:
        return pre_res, pd.DataFrame(), None, None

    out_df = (
        pd.DataFrame.from_dict(pre_res.results, orient="index")
        .reset_index(names="Term")
        .loc[:, ["Term", "fdr", "es", "nes"]]
        .sort_values("fdr", na_position="last")
        .reset_index(drop=True)
    )

    # Safely pick top term
    term_to_graph = out_df.loc[out_df["fdr"].notna(), "Term"].iloc[0]
    res = pre_res.results[term_to_graph]
    
    assert term_to_graph in pre_res.results
    assert isinstance(pre_res.ranking, (pd.Series, pd.DataFrame, np.ndarray))


    try: 
        gsea_plot = gseaplot(
            term=term_to_graph,
            hits=res["hits"],
            nes=res["nes"],
            pval=res["pval"],
            fdr=res["fdr"],
            RES=res["RES"],
            rank_metric=pre_res.ranking)
    except:
        gsea_plot = gseaplot2(
            terms=[term_to_graph],
            hits=[res["hits"]],
            RESs=[res["RES"]],
            rank_metric=pre_res.ranking,
        )

    return pre_res, out_df, term_to_graph, gsea_plot


import matplotlib.pyplot as plt

from math import ceil

def plot_gsea_grid(
    df_all,
    group_col="group",
    gene_sets="Allen_Brain_Atlas_10x_scRNA_2021",
    ncols=3,
    save_path=None,
    **gsea_kwargs
):
    """
    Runs GSEA for each unique group in a DataFrame and plots the top enriched terms in a grid.

    This function loops over all unique groups (e.g., clusters) present in `df_all`,
    calculates enrichment using `run_gsea_on_group`, extracts the single top-ranked pathway
    by FDR for each group, and renders their enrichment profiles side-by-side in a 
    scannable grid.

    Parameters
    ----------
    df_all : pandas.DataFrame
        The combined differential expression data containing all groups/clusters. Must
        contain the column specified in `group_col`.
    group_col : str, default="group"
        The name of the column in `df_all` containing the group/cluster assignments.
    gene_sets : str or list, default="Allen_Brain_Atlas_10x_scRNA_2021"
        The gene set database to use for GSEA enrichment.
    ncols : int, default=3
        Number of columns in the plotting grid.
    save_path : str, optional
        Filepath to save the final grid image (e.g., "gsea_grid.png").
    **gsea_kwargs : dict
        Additional keyword arguments (like `min_perm` or `seed`) passed directly 
        to `run_gsea_on_group`.

    Returns
    -------
    all_results : dict
        A dictionary mapping each group/cluster to its respective GSEA results structure:
        `{ group_id: {"pre_res": Prerank_object, "term": top_term_string} }`
    fig : matplotlib.figure.Figure
        The generated figure object containing the grid plot.
    """
    # 1. Run GSEA sequentially on each unique group
    all_results = {}
    unique_groups = df_all[group_col].unique()

    for group in unique_groups:
        group_df = df_all[df_all[group_col] == group]
        pre_res, out_df, term_to_graph, _ = run_gsea_on_group(
            group_df,
            gene_sets=gene_sets,
            **gsea_kwargs
        )

        if term_to_graph is not None:
            all_results[group] = {
                "pre_res": pre_res,
                "term": term_to_graph,
            }

    # If absolutely no groups had enrichments, stop early to avoid empty subplots
    n = len(all_results)
    if n == 0:
        print("Warning: No enriched pathways found for any group.")
        return all_results, None

    # 2. Calculate dynamic grid dimensions
    nrows = ceil(n / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(6 * ncols, 4 * nrows),
    )
    
    # Ensure axes is a flat array even if it's 1x1 or 1xN
    axes = np.array(axes).flatten()

    # 3. Plot enrichment metrics into individual grid panels
    for ax, (group, d) in zip(axes, all_results.items()):
        pre_res = d["pre_res"]
        term = d["term"]
        res = pre_res.results[term]

        # Plot Running Enrichment Score (RES)
        ax.plot(res["RES"], linewidth=2)

        # Plot Hit indicators at the bottom 15% of the y-axis
        ymin, ymax = ax.get_ylim()
        ax.vlines(
            res["hits"],
            ymin,
            ymin + (ymax - ymin) * 0.15,
            linewidth=0.4,
        )

        # Reference baseline
        ax.axhline(0, linestyle="--", linewidth=0.8)

        # Title containing enrichment statistics
        ax.set_title(
            f"Cluster {group}\n"
            f"{term[:50]}\n"
            f"NES={res['nes']:.2f}, FDR={res['fdr']:.2e}",
            fontsize=9,
        )
        ax.set_xlabel("Rank")
        ax.set_ylabel("ES")

    # 4. Hide empty panels (if n is not a perfect multiple of ncols)
    for ax in axes[n:]:
        ax.axis("off")

    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        
    plt.show()

    return all_results, fig

def parse_markers(adata, marker_table: pd.DataFrame, top_hvg: int = 5000):
    """
    Filter and split positive and negative marker genes from a metadata table 
    against the (optionally HVG-filtered) gene names in an AnnData object.

    This function parses comma-separated marker gene strings in a DataFrame, 
    validates them against the genes present in `adata` (optionally restricting 
    the search space to the top N highly variable genes), and maps them to their 
    respective cell types/classes.

    Parameters
    ----------
    adata : AnnData
        The AnnData object containing the single-cell dataset. If `top_hvg` is 
        provided, `adata.var` must contain a `'highly_variable_rank'` column.
    marker_table : pd.DataFrame
        A pandas DataFrame that contains at least the following columns:
        - 'class': Name of the cell type or cluster.
        - 'pos_markers': Comma-separated string of positive marker genes (or NaN).
        - 'neg_markers': Comma-separated string of negative marker genes (or NaN).
    top_hvg : int or None, default 5000
        The rank cutoff for highly variable genes to restrict the search. Only 
        genes with a `'highly_variable_rank'` less than or equal to this value 
        will be considered. If set to `None`, all genes in `adata.var_names` 
        are used.

    Returns
    -------
    pos_markers : dict
        A dictionary mapping each 'class' to a list of matching positive marker genes.
    neg_markers : dict
        A dictionary mapping each 'class' to a list of matching negative marker genes.
    not_found_genes_overall : list of str
        A list of all marker genes defined in the table that were either missing 
        from `adata` or filtered out because they exceeded the `top_hvg` threshold.
    """
    pos_markers = {}
    not_found_genes_overall = []

    if top_hvg is not None:
        var_names = set(adata.var_names[adata.var["highly_variable_rank"] <= top_hvg])  
    else:
        var_names = set(adata.var_names)

    for i, lst in enumerate(marker_table.pos_markers):
        found_genes = []
        not_found_genes = []
        if lst is not np.nan:
            for gene in lst.split(','):
                gene = gene.strip()
                if gene in var_names:
                    found_genes.append(gene)
                else:
                    not_found_genes.append(gene)
            
            pos_markers[marker_table['class'][i]] = found_genes
            not_found_genes_overall.extend(not_found_genes)
    print(f"Positive markers found in adata.var:")
    print(pos_markers)

    neg_markers = {}
    for i, lst in enumerate(marker_table.neg_markers):
        found_genes = []
        not_found_genes = []
        if lst is np.nan:
            neg_markers[marker_table['class'][i]] = []
            continue
        else:
            for gene in lst.split(','):
                gene = gene.strip()
                if gene in var_names:
                    found_genes.append(gene)
                else:
                    not_found_genes.append(gene)
            
            neg_markers[marker_table['class'][i]] = found_genes
            not_found_genes_overall.extend(not_found_genes)
    print(f"Negative markers found in adata.var:")
    print(neg_markers)
    print(f"Not found genes: {len(not_found_genes_overall)}")
    print(not_found_genes_overall)
    
    return pos_markers, neg_markers, not_found_genes_overall



# 1. PREPARE THE NET SCORING FUNCTION
def apply_dual_marker_scoring(adata, pos_dict, neg_dict):
    """
    Calculates a Net Score: Positive Marker Signal - Negative Marker Signal.
    """
    all_cell_types = list(pos_dict.keys())
    
    for cell_type in all_cell_types:
        # Get markers, ensuring they exist in the dataset
        pos_list = [g for g in pos_dict[cell_type] if g in adata.var_names]
        neg_list = [g for g in neg_dict.get(cell_type, []) if g in adata.var_names]
        
        # Calculate Positive Score
        if pos_list:
            sc.tl.score_genes(adata, gene_list=pos_list, score_name=f"{cell_type}_pos_score")
        else:
            adata.obs[f"{cell_type}_pos_score"] = 0
            
        # Calculate Negative Score (Penalty)
        if neg_list:
            sc.tl.score_genes(adata, gene_list=neg_list, score_name=f"{cell_type}_neg_score")
        else:
            adata.obs[f"{cell_type}_neg_score"] = 0
            
        # Compute Net Score
        # Logic: Net = Positive - Negative. 
        # We also apply a 'Relativity' clip: if negative signal is stronger than positive, score is 0.
        adata.obs[f"{cell_type}_net_score"] = adata.obs[f"{cell_type}_pos_score"] - adata.obs[f"{cell_type}_neg_score"]
        adata.obs.loc[adata.obs[f"{cell_type}_net_score"] < 0, f"{cell_type}_net_score"] = 0

    # 2. ASSIGN INITIAL LABELS
    net_score_cols = [f"{ct}_net_score" for ct in all_cell_types]
    
    # Find the best fitting cell type for each cell
    adata.obs['marker_label'] = adata.obs[net_score_cols].idxmax(axis=1).str.replace("_net_score", "")
    adata.obs['max_net_score'] = adata.obs[net_score_cols].max(axis=1)
    
    # 3. DEFINE HIGH-CONFIDENCE ANCHORS
    # A cell is "High Confidence" if it has a strong net score (e.g., > 0.1)
    # This excludes cells that have high negative marker expression or weak positive signals.
    adata.obs['reject_prediction'] = adata.obs['max_net_score'] < 0.5
    
    return adata

# Execute the scoring
adata = apply_dual_marker_scoring(adata, celltype_pos_markers, celltype_neg_markers)

