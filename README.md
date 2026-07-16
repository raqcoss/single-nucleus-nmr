# Naked Mole-Rat Brain snRNA-Seq Analysis Pipeline

This repository contains the bioinformatics pipeline and custom scripts developed for the processing, quality control, clustering, and cell-type annotation of single-nucleus RNA-sequencing (snRNA-seq) data derived from the naked mole-rat (*Heterocephalus glaber*) brain. This project is a colaboration between researchers from the National Autonomous University of Mexico (UNAM) and the Monterrey Institute of Technology and Higher Education (ITESM). 

---

## Table of Contents
- [Overview](#-overview)
- [Repository Structure](#-repository-structure)
- [Pipeline Workflow](#-pipeline-workflow)
  - [A. Preprocessing (From reads to high-quality count matrix)](#A_Preprocessing)
  - [B. Pipeline](#B_Pipeline)
  - [B1. Integration Benchmark and Clustering Optimization](#B_Optimization)
  - [B2. Annotation Methods](#B2_Annotation)
  - [C. Manual Refinement](#C_Refinement)
  - [Z. Archive](#Z_Archive)
- [Installation & Environment Setup](#-installation--environment-setup)
- [Usage Instructions](#-usage-instructions)
- [Contributing](#-contributing)
- [License](#-license)

---

## Overview

The naked mole-rat exhibits extraordinary physiological traits, including extreme longevity, cancer resistance, and hypoxia tolerance. Specifically, they presume to avoid neurodegeneration, which drives the motivation behind this study to discover the cellular mechanisms within the cerebral cortex, hippocampus, and midbrain regions using single-nucleus transcriptomics. 

This repository standardizes the workflow from raw `.fastq` reads to stable clusters and their cell-type annotations while addressing the unique challenges associated with analyzing the naked mole-rat transcriptome, such as the limited availability of a genome reference, ambient RNA contamination in nuclear preparations, and identifying novel cell types in this non-model organism.

---

## Repository Structure

```text
```directory

├── A_Preprocessing
│   ├── 01_sequencing_metrics.ipynb
│   ├── 02_cellbender_qc.ipynb
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

- Quality assessment of raw sequence reads to evaluate base qualities, GC content, adapter contamination, and duplication levels: `run_fasqc.sh`
- Alignment of reads to reference genome using `cellranger` from 10X Genomics:  `run_cellranger.sh`, `01_sequencing_metrics.ipynb`
- Denoising ambient RNA using `CellBender`: `run_cellbender.sh`using
- Filter low quality genes and cells: `main_preprocessing.ipynb' (Note: this notebook also performs the preprocessing and merge of human and NMR datasets to perform comparisons between species)

### B. Pipeline

Read alignment to the *Heterocephalus glaber* reference genome and generation of raw/filtered gene-expression matrices.

* **Tools:** `CellRanger` (v7.0.0 or specified version)
* **Script:** `scripts/02_cellranger_count.sh`
* *Note: Custom modifications made to the GTF file (e.g., handling of unannotated untranslated regions or lncRNAs) should be detailed here.*

### B1. Integration Benchmark and Clustering Optimization

Elimination of background/ambient RNA contamination typical in single-nucleus preparations using a deep generative model.

* **Tools:** `CellBender` v0.3.0 (`remove-background`)
* **Script:** `scripts/03_cellbender.sh`
* *Note: [Specify expected-cells and total-droplets parameters used here]*

### B2. Annotation Methods

Filtering of low-quality nuclei based on data-driven thresholds (e.g., minimum/maximum UMI counts, unique gene counts, and mitochondrial gene percentage bounds tailored to nuclear sequencing).

* **Tools:** `Scanpy` (Python)
* **Notebook:** `notebooks/04_downstream_qc.Rmd`

### C. Manual Refinement

Hyperparameter Optimization (HPO) for clustering resolution combined with an iterative, recursive clustering framework to separate major lineages (neurons, glia, vascular cells) down to distinct sub-populations.

* **Methodology:** [Detail your specific HPO or recursive split metrics here, e.g., Silhouette score, ROGUE, or resolution sweep]
* **Notebook:** `notebooks/05_recursive_clustering.Rmd`

### Z. Archive

Multi-pronged approach for cell-type classification using integrated automated tools alongside manual verification with canonical marker genes.

* **Tools Utilized:** [e.g., SingleR, Azimuth, CellTypist, or ScType]
* **Notebook:** `notebooks/06_cell_type_annotation.Rmd`

---

## Installation & Environment Setup

### Prerequisites

Ensure you have [Conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/) installed.

```bash
# Clone the repository
git clone [https://github.com/yourusername/nmr-brain-snrna-seq.git](https://github.com/yourusername/nmr-brain-snrna-seq.git)
cd nmr-brain-snrna-seq

# Create the conda environment
conda env create -f environment.yml

# Activate environment
conda activate nmr-snrna-env

```

---

## Usage Instructions

1. **Configure Paths:** Edit `config/config.yaml` to match your local cluster environment and dataset directories.
2. **Run Read QC:**
```bash
bash scripts/01_fastqc.sh

```


3. **Execute Mapping & Quantification:**
```bash
bash scripts/02_cellranger_count.sh

```


4. **Remove Background Contamination:**
```bash
bash scripts/03_cellbender.sh

```


5. **Interactive Analysis:** Open and run the step-by-step notebooks located within the `notebooks/` directory sequentially.

---

## Authors & Citation

* **Raquel Cossío Ramírez** - *Pipeline development* - [GitHub](https://github.com/raqcoss)

If you use this pipeline or code in your research, please cite:

> [Insert your publication citation or preprint info here]

If you have issues or doubts please [send me an email](mailto:raquel.cossior@gmail.com)
