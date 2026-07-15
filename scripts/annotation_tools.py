import gseapy as gp
from gseapy.plot import gseaplot, gseaplot2
import pandas as pd
import numpy as np

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