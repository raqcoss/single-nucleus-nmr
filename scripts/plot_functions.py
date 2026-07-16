import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_qc_violin_with_thres(adata, qc_cuts=None, region_col = "region"):
    """Plots violin plots with jittered data points and red threshold lines for single-cell QC metrics.

    This function generates a single row of four violin plots (number of genes, total counts, 
    mitochondrial percentage, and ribosomal percentage) for a given AnnData object. It overlays 
    individual cells as jittered points and draws horizontal dashed lines indicating the 
    filtering thresholds.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing cell-level metadata in `adata.obs`.
    qc_cuts : dict, optional
        Dictionary specifying the filtering thresholds. If None, defaults to:
        {
            "n_genes_by_counts": 3,
            "min_total_counts": 200,
            "max_total_counts": 22000,
            "pct_counts_mt": 3,
            "pct_counts_ribo": 7
        }
    region_col : str, default 'region'
        The column name in `adata.obs` representing the tissue or region. Used
        exclusively to extract the label for the print statement and the overall figure title.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created matplotlib figure object.
    axes : numpy.ndarray of matplotlib.axes.Axes
        An array containing the four subplots' axes.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc3k()
    >>> adata.obs['region'] = 'Peripheral Blood'
    >>> # Calculate QC metrics if not already present
    >>> sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    >>> fig, axes = plot_qc_violin_with_thres(adata)
    """
    if qc_cuts is None:
        qc_cuts = {"n_genes_by_counts":3, 
        "min_total_counts":200, 
        "max_total_counts":22000,
        "pct_counts_mt": 3, 
        'pct_counts_ribo': 7}
    tissue = adata.obs[region_col].iloc[0]
    print(f"Plotting {tissue}")
    plot_keys = ["n_genes_by_counts", "total_counts", "pct_counts_mt", 'pct_counts_ribo']
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    for i, key in enumerate(plot_keys):
        # Create violin plot with seaborn
        data_to_plot = adata.obs[[key]].copy()
        sns.violinplot(data=data_to_plot, y=key, ax=axes[i], color='skyblue')
        
        # Add jitter points on top
        sns.stripplot(data=data_to_plot, y=key, ax=axes[i], color='black', alpha=0.4, size=1, jitter=1 )
        
        # Add threshold lines
        if key == "total_counts":
            axes[i].axhline(y=qc_cuts['min_total_counts'], color='red', linestyle='--', linewidth=2, alpha=1, label=f'min: {qc_cuts["min_total_counts"]}')
            axes[i].axhline(y=qc_cuts['max_total_counts'], color='red', linestyle='--', linewidth=2, alpha=0.7, label=f'max: {qc_cuts["max_total_counts"]}')
            axes[i].legend()
        elif key in qc_cuts:
            threshold = qc_cuts[key]
            axes[i].axhline(y=threshold, color='red', linestyle='--', linewidth=2, alpha=0.8, label=f'threshold: {threshold}')
            axes[i].legend()
        
        axes[i].set_title(f'{key}', fontsize=12)
    
    fig.suptitle(f'{tissue}', fontsize=14, fontweight='bold')
    return fig, axes

import matplotlib.pyplot as plt
import scanpy as sc
from adjustText import adjust_text

def plot_highly_variable_genes(
    adata, 
    save_suffix: str, 
    n_top_genes: int = 20, 
    figsize: tuple = (6, 4)
):
    """
    Plot a scatter plot of mean vs. normalized variance, highlighting and labeling 
    the top highly variable genes.

    Parameters
    ----------
    adata : AnnData
        The AnnData object containing gene metadata and `'highly_variable_rank'` in `.var`.
    save_suffix : str
        Suffix for the saved file name (e.g., 'cortex', 'region_1'). 
        The file will be saved as 'scatter_{save_suffix}_hvg.png'.
    n_top_genes : int, default 20
        The number of top highly variable genes to label on the plot.
    figsize : tuple, default (6, 4)
        The width and height of the figure in inches.
    """
    import os

    # Identify the top highly variable genes based on rank
    tophvg = adata.var.sort_values('highly_variable_rank').head(n_top_genes).index.to_list()
    
    with plt.rc_context({"figure.figsize": figsize}):
        x = "means"
        y = "variances_norm"
        color = "is_highly_variable"

        # Ensure the color column exists as a string categorical for discrete coloring
        adata.var["is_highly_variable"] = (
            adata.var["highly_variable"]
            .astype(bool)
            .astype(str)
        )

        # Generate the base scatter plot
        ax = sc.pl.scatter(
            adata,
            x=x,
            y=y,
            color=color,
            show=False,
            size=15
        )

        # Log-transform x axis
        ax.set_xscale("log")

        # Move plot title from Axes to Legend
        ax.set_title("")
        legend = ax.get_legend()
        if legend is not None:
            legend.set_title("Highly Variable")

        # Place the text labels for the top genes
        texts = []
        for gene in tophvg:
            x_loc = adata.var.loc[gene, x]
            y_loc = adata.var.loc[gene, y]

            texts.append(
                ax.text(
                    x_loc,
                    y_loc,
                    gene,
                    color="k",
                    fontsize=10
                )
            )

        # Reposition labels to prevent overlaps
        adjust_text(
            texts,
            expand=(1.2, 1.2),
            arrowprops=dict(color="gray", lw=1),
            ax=ax,
        )
        
        # Save figure
        os.makedirs("figures", exist_ok=True)
        save_path = f"figures/scatter_{save_suffix}_hvg.png"
        plt.savefig(
            save_path,
            dpi=300,
            bbox_inches="tight"
        )
        print(f"Plot saved successfully to {save_path}")

        plt.show()

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plot_3d_pca_gene(
    adata, 
    gene: str, 
    save_prefix: str, 
    cmap: str = "viridis", 
    figsize: tuple = (8, 7)
):
    """
    Plot cells in a 3D PCA space colored by the expression of a specific gene.

    Parameters
    ----------
    adata : AnnData
        The AnnData object containing 'X_pca' in `.obsm` and the target gene in `.var_names`.
    gene : str
        The name of the gene to color the cells by (e.g., 'ADARB2').
    save_prefix : str
        Prefix for the saved file name (e.g., 'cortex' or 'region_1'). 
        The file will be saved as '{save_prefix}_PCA_3D_{gene}.png'.
    cmap : str, default 'viridis'
        The colormap used to represent gene expression intensity.
    figsize : tuple, default (8, 7)
        The width and height of the figure in inches.

    Raises
    ------
    ValueError
        If 'X_pca' is missing from `.obsm` or the specified gene is not found in `.var_names`.
    """
    # 1. Verify PCA coordinates exist
    if "X_pca" not in adata.obsm:
        raise ValueError("PCA coordinates 'X_pca' not found in adata.obsm. Please run sc.tl.pca first.")
        
    # Extract first 3 PCs
    pcs = adata.obsm["X_pca"][:, :3]

    # 2. Extract and format gene expression values safely
    if gene in adata.var_names:
        expr = adata[:, gene].X
        # Convert sparse matrix (csr_matrix) to dense numpy array if needed
        if hasattr(expr, "toarray"):
            expr = expr.toarray()
        expr = np.ravel(expr)
    else:
        raise ValueError(f"Gene '{gene}' not found in adata.var_names")

    # 3. Create the 3D plot
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')

    scatter = ax.scatter(
        pcs[:, 0],
        pcs[:, 1],
        pcs[:, 2],
        c=expr,
        cmap=cmap,
        s=3
    )

    # 4. Axis labels and titles
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")
    ax.set_title(f"3D PCA colored by {gene} expression ({save_prefix})")

    # 5. Colorbar
    cbar = plt.colorbar(scatter, ax=ax, pad=0.1)
    cbar.set_label(f"{gene} expression")

    plt.tight_layout()

    # 6. Save and show the figure
    save_path = f"{save_prefix}_PCA_3D_{gene}.png"
    plt.savefig(
        save_path,
        dpi=300,
        bbox_inches="tight"
    )
    print(f"Plot saved successfully to {save_path}")
    
    plt.show()


from matplotlib.lines import Line2D
def umap_alpha(adata, color_col='leiden', alpha_col ='total_counts', size = 5, normalize_a= True, embedding = "X_umap"):
    """Plots a UMAP scatter visualization where cell transparency (alpha) is driven by a continuous metadata value.

    This function displays an embedding space (e.g., UMAP) using cluster colors assigned to cells
    while varying the transparency (`alpha`) of individual points based on a continuous variable in 
    `adata.obs` (such as total counts, QC metrics, or marker expression values).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing the embedding in `adata.obsm` and categorical/continuous 
        variables in `adata.obs`.
    color_col : str, default 'leiden'
        The categorical metadata column in `adata.obs` used to group and color cells.
    alpha_col : str, default 'total_counts'
        The continuous metadata column in `adata.obs` used to determine individual cell alpha values.
    size : float, default 5
        Marker size of the individual scatter points.
    normalize_a : bool, default True
        If True, normalizes `alpha_col` values linearly to fall between a visible minimum (0.05) and 
        maximum (1.0) range to prevent high values from dominating or zero-values from completely vanishing.
    embedding : str, default 'X_umap'
        The coordinate multidimensional array key stored in `adata.obsm` (e.g., 'X_umap', 'X_tsne').

    Raises
    ------
    KeyError
        If `adata.uns` does not contain the color palette key `{color_col}_colors` generated by 
        previous scanpy plotting steps (e.g. `leiden_colors` for `leiden`).

    Notes
    -----
    - This function requires that cluster color mapping has already been configured in `adata` by 
      running a standard plotting command first (e.g., `sc.pl.umap(adata, color=color_col)`), which 
      generates the color array in `adata.uns["{color_col}_colors"]`.
    - It assumes that cluster labels can be sorted numerically after conversion to string/integer 
      for standard legend ordering.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc3k()
    >>> sc.pp.neighbors(adata)
    >>> sc.tl.umap(adata)
    >>> sc.tl.leiden(adata)
    >>> # Plot once to generate 'leiden_colors' in adata.uns
    >>> sc.pl.umap(adata, color='leiden', show=False)
    >>> # Plot with transparency driven by UMI count
    >>> umap_alpha(adata, color_col='leiden', alpha_col='total_counts')
    """
    emb = adata.obsm[embedding]
    clusters = adata.obs[color_col]      # or any cluster column
    values = adata.obs[alpha_col]  # column that drives alpha
        # normalize alpha into [0.05, 1]
    alpha = values.values.astype(float)

    if normalize_a:
        alpha = (alpha - alpha.min()) / (alpha.max() - alpha.min() + 1e-9)
        alpha = 0.05 + 0.95 * alpha

    # assign colors per cluster (Scanpy already stores them after sc.pl.umap)
    palette = adata.uns["leiden_colors"]
    # Convert cluster categories to integers for sorting
    cluster_labels = list(clusters.unique())
    cluster_labels_int = sorted(cluster_labels, key=lambda x: int(x))

        # Reorder colors according to numeric order
    color_map = dict(zip(cluster_labels, palette))
    ordered_colors = [color_map[k] for k in cluster_labels_int]

    fig, ax = plt.subplots(figsize=(7, 6))

    # plot scatter in the same numeric order
    for clust in cluster_labels_int:
        color = color_map[clust]
        idx = clusters == clust
        ax.scatter(
            emb[idx, 0],
            emb[idx, 1],
            s=size,
            c=color,
            alpha=alpha[idx],
            linewidth=0,
        )

    # build proxy legend
    handles = [
        Line2D([0], [0], marker="o", color=color, markersize=6, linestyle="",
            markerfacecolor=color, markeredgewidth=0)
        for color in ordered_colors
    ]

    ax.legend(
        handles,
        cluster_labels_int,
        title="Clusters",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False
    )

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"UMAP with per-cell alpha (driven by {alpha_col})")

    plt.tight_layout()
    plt.show()


def boxplot_obs_by_cluster(adata, cluster_key, columns):
    """Generates sequential boxplots showing the distribution of numerical metadata columns grouped by cluster.

    This utility iterates through a list of continuous cell-level metadata columns (e.g., quality control metrics 
    or classification scores) and creates an individual boxplot figure for each. The clusters on the x-axis 
    are dynamically sorted in ascending numerical order.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix containing cell-level metadata in `adata.obs`.
    cluster_key : str
        The column name in `adata.obs` representing cell cluster assignments (e.g., "leiden", "subcluster").
        Values in this column will be cast to strings and sorted numerically.
    columns : list of str
        The numerical columns in `adata.obs` whose distributions will be visualized (e.g., `['total_counts', 'pct_counts_mt']`).

    Returns
    -------
    None
        This function directly displays the generated matplotlib plots inline using `plt.show()` 
        and does not return any objects.

    Raises
    ------
    KeyError
        If `cluster_key` or any of the specified `columns` are missing from `adata.obs`.
    ValueError
        If cluster names cannot be cleanly parsed as integers (e.g., alphanumeric strings like "A1", "B1") 
        during the numerical sorting process.

    Notes
    -----
    - Missing or `NaN` values within the specified `columns` are automatically dropped before plotting 
      on a per-cluster basis.
    - The boxes are filled with a solid default color (`patch_artist=True`), and a light horizontal grid 
      is added to improve readability.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc3k()
    >>> sc.pp.calculate_qc_metrics(adata, inplace=True)
    >>> sc.tl.leiden(adata)
    >>> # Generates two distinct figures sequentially
    >>> boxplot_obs_by_cluster(
    ...     adata, 
    ...     cluster_key="leiden", 
    ...     columns=["total_counts", "n_genes_by_counts"]
    ... )
    """
    df = adata.obs[[cluster_key] + columns].copy()
    df[cluster_key] = df[cluster_key].astype(str)

    # sort cluster labels numerically
    cluster_order = sorted(df[cluster_key].unique(), key=lambda x: int(x))

    for col in columns:
        plt.figure(figsize=(10, 5))
        plt.title(f"Distribution of {col} by cluster")

        data = [df.loc[df[cluster_key] == cl, col].dropna().values
                for cl in cluster_order]

        plt.boxplot(
            data,
            labels=cluster_order,
            notch=False,
            patch_artist=True
        )

        plt.xlabel("Cluster")
        plt.ylabel(col)
        plt.grid(axis="y", linestyle="--", alpha=0.4)
        plt.tight_layout()
        plt.show()


def plot_pc_qc_correlation(corr, save_path=None):
    """Plots a heatmap of Spearman correlation coefficients between PCs and QC covariates.

    Parameters
    ----------
    corr : pandas.DataFrame
        Long-format correlation results returned by `check_pc_correlation` containing 
        `'PC'`, `'covariate'`, and `'r'` columns.
    save_path : str, optional
        File path where the generated figure should be saved (e.g., `'plots/pc_qc_corr.png'`). 
        If None, the plot is not written to disk.

    Returns
    -------
    None
        Displays the generated heatmap inline using `plt.show()`.

    Examples
    --------
    >>> corr_df = check_pc_correlation(adata)
    >>> plot_pc_qc_correlation(corr_df, save_path="qc_correlation_heatmap.pdf")
    """
    mat = corr.pivot(index="covariate", columns="PC", values="r")
    plt.figure(figsize=(10,3))
    im = plt.imshow(mat.values, aspect="auto", interpolation="nearest")
    plt.yticks(range(mat.shape[0]), mat.index)
    plt.xticks(range(mat.shape[1]), mat.columns, rotation=90)
    plt.colorbar(im, label="Spearman r")
    plt.title("Correlation PC vs QC")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    plt.show()

def plot_composition(adata, region, celltype_col, save_path=None):
    """
    Plots a single stacked bar representing the composition of cell types in an AnnData object.

    This function calculates cell type proportions, filters out unobserved categories 
    (crucial for subsets), maps them to their original color palette from `adata.uns`, 
    and displays a stacked bar chart annotated with percentages and a legend showing raw counts.

    Parameters
    ----------
    adata : AnnData
        The annotated data object containing single-cell data.
    region : str
        The name of the region or sample being plotted (used for the plot title).
    celltype_col : str
        The column in `adata.obs` containing cell type classifications.
    save_path : str, optional
        The file path where the plot should be saved. If None, the plot is not saved.
    """
    # 1. Map colors using the stable category-to-color mapping first
    color_map = dict(
        zip(
            adata.obs[celltype_col].cat.categories,
            adata.uns[f"{celltype_col}_colors"],
        )
    )

    # 2. Count cells and explicitly drop zero-count categories for subsets
    counts = (
        adata.obs[celltype_col]
        .value_counts()
        .loc[lambda x: x > 0]  # <-- THE FOOLPROOF FIX FOR SUBSETS
        .sort_values(ascending=False)
    )

    cluster_order = counts.index.tolist()
    proportions = counts / counts.sum()

    # 3. Pull colors cleanly matching only the observed cell types
    colors = [color_map.get(c, "lightgray") for c in cluster_order]

    # 4. Plot one stacked bar
    fig, ax = plt.subplots(figsize=(3, 10))
    bottom = 0

    for celltype, prop, raw, color in zip(
        cluster_order, proportions, counts, colors
    ):
        ax.bar(
            x=0, height=prop, bottom=bottom, color=color, width=0.8, label=celltype
        )

        # annotate percentage
        if prop >= 0.008:
            ax.text(
                0,
                bottom + prop / 2,
                f"{prop*100:.1f}%",
                ha="center",
                va="center",
                fontsize=6,
            )
        bottom += prop

    # 5. Formatting
    ax.set_xticks([0])
    ax.set_xticklabels(["All cells"])
    ax.set_ylabel("Proportion of cells")
    ax.set_title(f"{celltype_col.capitalize()} composition ({region})")

    # 6. Legend with counts (order perfectly matching the observed subset)
    handles, labels = ax.get_legend_handles_labels()
    new_labels = [f"{label}: {counts[label]}" for label in labels]

    ax.legend(handles, new_labels, bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    plt.show()

def plot_marker_genes_umap(
    adata, 
    filtered_markers: dict, 
    group_by: str = 'leiden', 
    ncols: int = 6, 
    color_map: str = 'YlOrBr'
):
    """
    Plot UMAPs of a grouping variable (e.g., leiden clusters) alongside specified marker genes.

    This function automatically formats the plot titles to show which cluster/group 
    each marker gene belongs to (e.g., 'Cluster_0: GeneA') and plots them in a grid.

    Parameters
    ----------
    vdata : AnnData
        The AnnData object containing the single-cell dataset and UMAP coordinates.
    filtered_markers : dict
        A dictionary where keys represent groups/clusters and values are lists of 
        marker genes associated with that group (e.g., {'T-cells': ['CD3D', 'CD3E']}).
    group_by : str, default 'leiden'
        The categorical observation key in `vdata.obs` to plot alongside the genes.
    ncols : int, default 6
        Number of columns in the multi-panel UMAP grid.
    color_map : str, default 'YlOrBr'
        The colormap to use for continuous gene expression values.
    """
    import scanpy as sc
    # Initialize lists with the main categorical grouping (e.g., 'leiden')
    chosen_genes = [group_by]
    titles = [group_by]
    
    # Dynamically build the list of genes to plot and their corresponding titles
    for group_key, gene_list in filtered_markers.items():
        chosen_genes.extend(gene_list)
        titles.extend([f"{group_key}: {gene}" for gene in gene_list])
    
    # Generate the multi-panel UMAP plot
    sc.pl.umap(
        adata, 
        color=chosen_genes, 
        title=titles, 
        ncols=ncols, 
        color_map=color_map
    )

import os
import matplotlib.pyplot as plt
import seaborn as sns

def plot_score_distribution_by_group(
    adata,
    group_name: str = "marker_label",
    score_column: str = "max_net_score",
    threshold: float = 0.5,
    save_prefix: str = "region",
    save_dir: str = "figures",
    palette: str = "Spectral",
    figsize: tuple = (14, 7),
):
    """
    Generate a hybrid Violin-Strip plot showing the distribution of scores across
    different cell groups, including sample size annotations and a custom threshold line.

    Parameters
    ----------
    adata : AnnData
        The AnnData object containing metadata in `.obs`.
    group_name : str, default 'marker_label'
        The categorical observation column in `adata.obs` to group cells by.
    score_column : str, default 'max_net_score'
        The numerical column in `adata.obs` representing scores to plot on the y-axis.
    threshold : float, default 0.5
        The numerical cut-off threshold line to display horizontally on the plot.
    save_prefix : str, default 'region'
        Name prefix (such as region/sample name) used in the plot title and output filename.
    save_dir : str, default 'figures'
        Directory where the generated plot will be saved.
    palette : str, default 'Spectral'
        Seaborn palette used for coloring the violin plots.
    figsize : tuple, default (14, 7)
        The width and height of the generated figure.
    """
    # 1. Calculate cell counts and determine plotting order (sorted by group size)
    group_counts = adata.obs.groupby(group_name)[score_column].count().sort_values(ascending=False)
    order = group_counts.index
    n_cells_per_type = group_counts.values

    # Set up figure
    fig, ax = plt.subplots(figsize=figsize)

    # 2. Draw Violin Plot
    sns.violinplot(
        data=adata.obs, 
        x=group_name, 
        y=score_column, 
        order=order, 
        palette=palette,
        inner='quartile',
        cut=0,
        alpha=0.4,
        ax=ax
    )

    # 3. Draw Jittered Individual Cells (Strip Plot)
    sns.stripplot(
        data=adata.obs, 
        x=group_name, 
        y=score_column, 
        order=order, 
        color='black', 
        size=1.5,
        alpha=0.3,
        jitter=True,
        ax=ax
    )

    # 4. Add dynamic cell counts above each violin plot
    for i, (cell_type, count) in enumerate(zip(order, n_cells_per_type)):
        max_val = adata.obs[adata.obs[group_name] == cell_type][score_column].max()
        # Fallback to 0 if max is null or empty
        max_val = max_val if not pd.isna(max_val) else 0.0
        
        ax.text(
            i, 
            max_val + 0.05, 
            f'{count}', 
            ha='center', 
            va='bottom', 
            fontsize=10, 
            fontweight='bold'
        )

    # 5. Add Horizontal Cutoff Threshold
    ax.axhline(
        threshold, 
        ls='--', 
        color='red', 
        alpha=0.8, 
        label=f'{score_column.replace("_", " ").title()} Threshold ({threshold:.2f})', 
        linewidth=2
    )

    # 6. Apply Styling and Labels
    ax.set_title(f'{save_prefix} Max Net Scores Distribution with Individual Cells', fontsize=16)
    ax.set_ylabel(score_column.replace('_', ' ').title(), fontsize=12)
    ax.set_xlabel('Assigned Cell Type', fontsize=12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.legend(bbox_to_anchor=(1.00, 1), loc='upper left')

    plt.tight_layout()

    # 7. Safe Save and Export
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, f'{save_prefix}_score_distribution_jitter.png')
        plt.savefig(save_path, dpi=300)
        print(f"Plot successfully saved to: {save_path}")

    plt.show()
    import os
import matplotlib.pyplot as plt
import pandas as pd

def plot_annotation_confidence_proportions(
    adata,
    group_col: str = "marker_label",
    confidence_col: str = "is_confident",
    threshold: float = 0.5,
    save_prefix: str = "region",
    save_dir: str = "figures",
    figsize: tuple = (13, 7),
):
    """
    Generate a stacked bar plot showing the proportion of confident annotations 
    per cell type, with total cell counts and success percentages annotated above each bar.

    Parameters
    ----------
    adata : AnnData
        The AnnData object containing metadata in `.obs`.
    group_col : str, default 'marker_label'
        The column in `adata.obs` to group on the x-axis (e.g., cell types).
    confidence_col : str, default 'is_confident'
        The boolean or categorical column indicating whether the annotation is confident/acceptable.
    threshold : float, default 0.5
        The cell score threshold value to display in the plot's title.
    save_prefix : str, default 'region'
        Name prefix (such as region/sample name) used in the plot title and output filename.
    save_dir : str, default 'figures'
        Directory where the generated plot will be saved. Set to None to skip saving.
    figsize : tuple, default (13, 7)
        The width and height of the generated figure.
    """
    # 1. Prepare data: Create a cross-tabulation of groups vs. confidence
    stack_data = pd.crosstab(adata.obs[group_col], adata.obs[confidence_col])

    # 2. Rename columns dynamically to reflect confidence levels
    # Assumes a boolean structure: False (Index 0) -> Not Acceptable, True (Index 1) -> Acceptable
    stack_data.columns = ['Not Acceptable', 'Acceptable']

    # 3. Calculate totals for sorting
    stack_data['Total'] = stack_data.sum(axis=1)
    stack_data = stack_data.sort_values('Total', ascending=False)

    # 4. Remove 'Total' to isolate bars for plotting
    plot_df = stack_data.drop(columns=['Total'])

    # 5. Plot vertical stacked bars
    fig, ax = plt.subplots(figsize=figsize)
    plot_df.plot(
        kind='bar', 
        stacked=True, 
        color=['#E69F00', '#0072B2'], 
        edgecolor='white',
        ax=ax
    )

    # 6. Apply text labels above the bars & decorative hatches to the 'Not Acceptable' slice
    n_categories = len(stack_data)
    for i, (name, row) in enumerate(stack_data.iterrows()):
        total = row['Total']
        anchor = row['Acceptable']
        percentage = (anchor / total) * 100 if total > 0 else 0.0
        
        # Position label text comfortably above the bar peak
        ax.text(
            i, 
            total + (total * 0.02), 
            f'{int(anchor)}/{int(total)}\n({percentage:.1f}%)', 
            ha='center', 
            va='bottom', 
            fontsize=10, 
            fontweight='regular'
        )
        
        # Apply texture pattern safely to the base (Not Acceptable) bars
        if i < len(ax.patches):
            ax.patches[i].set_hatch('//')

    # Styling configurations
    ax.set_title(
        f'Proportion of Confident Annotations per Celltype in {save_prefix}\n'
        f'(min cell type score threshold = {threshold:.2f})', 
        fontsize=16, 
        pad=25
    )
    ax.set_ylabel('Number of cells', fontsize=12)
    ax.set_xlabel('Cell type', fontsize=12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.legend(title='Confidence Level', frameon=False)

    # Allocate a 15% head room above the tallest bar for annotations
    ax.set_ylim(0, stack_data['Total'].max() * 1.15)

    plt.tight_layout()

    # 7. Safe Save and Export
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, f'{save_prefix}_conf_prop_per_celltype.png')
        plt.savefig(save_path, dpi=300)
        print(f"Plot successfully saved to: {save_path}")

    plt.show()