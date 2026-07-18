# Naked Mole-Rat Brain snRNA-Seq Analysis Pipeline

This repository contains the bioinformatics pipeline and custom scripts developed for the processing, quality control, clustering, and cell-type annotation of single-nucleus RNA-sequencing (snRNA-seq) data derived from the naked mole-rat (*Heterocephalus glaber*) brain. This project is part of a colaboration between researchers from the National Autonomous University of Mexico (UNAM) and the Monterrey Institute of Technology and Higher Education (ITESM). 

---

## Table of Contents
- [Overview](#-overview)
- [Repository Structure](#-repository-structure)
- [Pipeline Workflow](#-pipeline-workflow)
  - [A. Preprocessing (From reads to high-quality count matrix)](#a-preprocessing)
  - [B. Pipeline](#b-pipeline)
  - [B1. Integration Benchmark and Clustering Optimization](#b1-optimization)
  - [B2. Annotation Methods](#b2-annotation)
    - [Cluster-Aware Dual Marker Score (CADuMS) Annotation](#cluster-aware-dual-marker-score-cadums-annotation)
    - [TrainCell](#traincell)
    - [Cross Species Label Transfer (XSLT)](#cross-species-label-transfer-xslt)
  - [C. Manual Refinement](#c-manual-refinement)
  - [Z. Archive](#z-archive)
- [Installation & Environment Setup](#-installation--environment-setup)
- [Usage Instructions](#-usage-instructions)
- [Contributing](#-contributing)
- [License](#-license)

---

## Overview

The naked mole-rat exhibits extraordinary physiological traits, including extreme longevity, cancer resistance, and hypoxia tolerance. Specifically, they presume to avoid neurodegeneration, which drives the motivation behind this study to discover the cellular mechanisms within the cerebral cortex, hippocampus, and midbrain regions using single-nucleus transcriptomics. 


This repository standardizes the workflow from raw `.fastq` reads to stable clusters and their cell-type annotations. It also addresses the unique challenges of analyzing the naked mole-rat transcriptome, such as the limited availability of a genome reference, ambient RNA contamination in nuclear preparations, and identifying novel cell types in this non-model organism.

---

## Repository Structure

```text
```directory

├── A_Preprocessing
│   ├── 01_sequencing_metrics.ipynb
│   ├── 02_cellbender_qc.ipynb
│   ├── 03_preprocessing.ipynb
│   └── orthology/
├── B1_Optimization
│   ├── 01_batch_correction_bm.ipynb
│   ├── 02_recursive_subclustering.ipynb
│   └── figures/
├── B2_Annotation
│   ├── CADuMS/
│   ├── CellTrain/
│   ├── celltypist/
│   ├── gseapy/
│   ├── scanvi/
│   └── XSLT/
├── B_Pipeline
│   ├── 01_classic_processing.ipynb
│   ├── 02_optimal_processing.ipynb
│   └── 03_join_clusters_and_annotations.ipynb
├── C_Refinement
│   ├── 01_cortex_cadmus_refinement.ipynb
│   ├── 02_hippocampus_cadmus_refinement.ipynb
│   ├── 03_saving_midbrain.ipynb
│   ├── 04_midbrain_annotation.ipynb
│   ├── 05_all_region_refinement.ipynb
│   ├── 06_check_refined_annotations.ipynb
├── data/
├── envs/                         # Conda environment specifications
├── LICENSE
├── pyproject.toml
├── README.md
├── results/
├── scripts/                      # Quick acess to custom functions and scripts
└── Z_Archive      
    └── seurat/                   # Alternative analysis in seurat

```

---

## Pipeline Workflow

### A. Preprocessing (From reads to high-quality count matrix)
Read alignment to the *Heterocephalus glaber* reference genome and generation of raw/filtered gene-expression matrices.

- Quality assessment of raw sequence reads to evaluate base qualities, GC content, adapter contamination, and duplication levels: `run_fasqc.sh`
  
- Alignment of reads to reference genome using `cellranger` from 10X Genomics:  `run_cellranger.sh`, `01_sequencing_metrics.ipynb`
  
- Elimination of background RNA contamination typical in single-nucleus preparations using a deep generative model, `CellBender`: `run_cellbender.sh`, `02_cellbender_qc.ipynb`

- Filter low quality genes and cells: `preprocessing.ipynb`

### B. Pipeline
- Classic Pipeline (Normalize, HVG, Scale, PCA, KNN, Clustering & UMAP): `01_classic_processing.ipynb`
  
- Improved Pipeline after choosing the most adequate integration technique and applying recursive subclustering (both specified in B1): `02_optimal_processing.ipynb`
  
- Join defined subclusters and annotations (specified in B2, selected method: CADuMS Annotation): `03_join_clusters_and_annotations.ipynb`

### B1. Integration Benchmark and Clustering Optimization

Three integration methods (to correct the effect given by different replicates) were evaluated: Harmony, BBKNN and scVI. Although unintegrated dataset achieved acceptable replicate mixing and low cell-type mixing (both desirable characteristics) the Harmony integration slightly increased the replicate mixing without affecting the biological conservation given by the cell-type mixing. This procedure is shown in `01_batch_correction_bm.ipynb`.

Then, a recursive clustering pipeline was developed to derive 3-level grouppings, first "superclusters", then "clusters" and finally "subclusters". An additional hyper-parameter optimization was implemented within each clustering step to adjust the *resolution* parameter in the Leiden algorithm to achive the larger silhouette score. This techinique is detailed in `02_recursive_subclustering.ipynb`

### B2. Annotation Methods

Several annotaton methods were used to define cell-types, including celltypist, GSEA and scANVI. Additionally some custom annotation methods were developed to better address the challenges derived from unmatched species references on extrapolating the expected gene expression for every cell-type label on the naked-mole rat brain dataset. 

The following methods were developed within this study. Please cite them accordingly (Article publication still pending).

#### Cluster-Aware Dual Marker Score (CADuMS) Annotation

A maker-based method that uses gene scores per panel of possitively and negatively associated markers to a cell-type diven a set of expected celltype classes, and defines the max-scored celltype per cell. Then unconfidently matched cells are marked if the `max total score` is lower than the smallest 20th percentile. After that, a consensus within each subcluster is applied considering the following:

- Homogeneity: Most cells in the subcluster belong to mostly one cell-type.
- Dual Identity: The second majority cell-type is at least 90% of the fist majority cell-type.
- Group Confidence: The number of unconfident cells within a subcluster is less than 30%.

This way class labels are merged (with a ` | ` sign) along dual identity groups and marked with an asterisk (*) if they are not confident. Above that, if the group can not comply with an homogeneus composition it is called "Ambiguous" and if it also fails on group confidence or has dual identity it remains as "Unassigned".

##### Method Requirements

- Biological knowledge of the sample to define expected cell-types.
- Table of possitively and negatively associated gene markers per expected cell-types
- Manual Parameter setting

* **Notebook:** `B2_Annotation/CADuMS/CADuMS.ipynb`

#### CellTrain
This method is also a marker-based method which automatically defines the cell identities per each data point. Then a logistic regression model is trained and tested with this data and later used to predict the same original dataset. Although validation infers that the model is learning the structure of the data, it is not fully recommended as it introduces many sources of bias on the second prediction of the same dataset. This method is displayed here as precedent to hopefully be used on other future naked mole rat datasets, after having validated the true naked-mole rat brain cell-types. An important note: This method is derived from the *celltypist* method described in Domínguez Conde *et al.*, 2022.

#### Cross Species Label Transfer (XSLT)
Finally an integration of the annotated dataset for human adult and developing brain from the Linnarson's Lab with the naked-mole rat dataset was developed to perform a label transfer on similar cells defined to be in the same cluster as the other species cells. Unfortunately many cells of the naked-mole rat remained unannotaded with this method, but it serve as precedent of comparing the two species and can be further used to compare gene expression on similar cell-types across species.

### C. Manual Refinement

Manual refinement was performed using as input the labels obtained from the CADuMS annotation tool, then per each group, the gene expression was evaluated per subluster on each class to see a congruent marker gene expression. Most expressed genes within the group were also investigated to validate their identity and resolve dual identities or unconfident annotations. 
This refinement was first performed on each region independently and then all regions were analysed together to define the final annotation. In the case of midbrain, only few cells could be recovered so a deep custom cleaning is performed in `03_saving_midbrain.ipynb` and a manual annotation in `04_midbrain_annotation.ipynb`

### Z. Archive

This data was also analysed with Seurat package in R. You can find the old scripts in this folder.

---

## Installation & Environment Setup

### Prerequisites

Ensure you have [Conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/) installed.

```bash
# Clone the repository
git clone [https://github.com/raqcoss/single-nucleus-nmr.git](https://github.com/raqcoss/single-nucleus-nmr.git)
cd single-nucleus-nmr

```
## References

1. Domínguez Conde, C., Xu, C., Jarvis, L., Rainbow, D., Wells, S., Gomes, T., Howlett, S., Suchanek, O., Polanski, K., King, H., Mamanova, L., Huang, N., Szabo, P., Richardson, L., Bolt, L., Fasouli, E., Mahbubani, K., Prete, M., Tuck, L., … Teichmann, S. (2022). Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science (New York, N.Y.), 376(6594), eabl5197. https://doi.org/10.1126/science.abl5197
2. Xu, C., Prete, M., Webb, S., Jardine, L., Stewart, B. J., Hoo, R., He, P., Meyer, K. B., & Teichmann, S. A. (2023). Automatic cell-type harmonization and integration across Human Cell Atlas datasets. Cell, 186(26), 5876-5891.e20. https://doi.org/10.1016/j.cell.2023.11.026
3. Braun, E., Danan-Gotthold, M., Borm, L. E., Lee, K. W., Vinsland, E., Lönnerberg, P., Hu, L., Li, X., He, X., Andrusivová, Ž., Lundeberg, J., Barker, R. A., Arenas, E., Sundström, E., & Linnarsson, S. (2023). Comprehensive cell atlas of the first-trimester developing human brain. Science, 382(6667), eadf1226. https://doi.org/10.1126/science.adf1226
4. Siletti, K., Hodge, R., Mossi Albiach, A., Lee, K. W., Ding, S.-L., Hu, L., Lönnerberg, P., Bakken, T., Casper, T., Clark, M., Dee, N., Gloe, J., Hirschstein, D., Shapovalova, N. V., Keene, C. D., Nyhus, J., Tung, H., Yanny, A. M., Arenas, E., … Linnarsson, S. (2023). Transcriptomic diversity of cell types across the adult human brain. Science, 382(6667), eadd7046. https://doi.org/10.1126/science.add7046


## Authors & Citation

* **Raquel Cossío Ramírez** - *Pipeline development & alternative annotation methods development*  - [GitHub](https://github.com/raqcoss)

If you use this pipeline or code in your research, please cite:

> [Insert your publication citation or preprint info here]

If you have issues or doubts please [send me an email](mailto:raquel.cossior@gmail.com)
