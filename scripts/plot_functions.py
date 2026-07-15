import pandas as pd
import matplotlib.pyplot as plt
# import scanpy as sc

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