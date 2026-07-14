# Naked Mole-Rat Brain snRNA-Seq Analysis Pipeline

This repository contains the bioinformatics pipeline and custom scripts developed for the processing, quality control, clustering, and cell-type annotation of single-nucleus RNA-sequencing (snRNA-seq) data derived from the naked mole-rat (*Heterocephalus glaber*) brain. This project is a colaboration between researchers from the National Autonomous University of Mexico (UNAM) and the Monterrey Institute of Technology and Higher Education (ITESM). 

---

## Table of Contents
- [Overview](#-overview)
- [Repository Structure](#-repository-structure)
- [Pipeline Workflow](#-pipeline-workflow)
  - [1. Raw Data Quality Control (FastQC)](#1-raw-data-quality-control-fastqc)
  - [2. Alignment and Quantification (CellRanger)](#2-alignment-and-quantification-cellranger)
  - [3. Ambient RNA Removal (CellBender)](#3-ambient-rna-removal-cellbender)
  - [4. Downstream Quality Control & Filtering](#4-downstream-quality-control--filtering)
  - [5. HPO & Recursive Clustering](#5-hpo--recursive-clustering)
  - [6. Cell-Type Annotation](#6-cell-type-annotation)
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


```

```text
Markdown template generated successfully.

```directory
├── README.md
├── LICENSE
├── environment.yml               # Conda environment specifications
├── config/
│   └── config.yaml               # Pipeline hyperparameters and paths
├── data/
│   ├── reference/                # Naked mole-rat genome assembly and GTF
│   └── raw/                      # Placeholder or symlinks for raw FASTQ data
├── scripts/
│   ├── 01_fastqc.sh              # Quality control on raw reads
│   ├── 02_cellranger_count.sh    # Reference building and cellranger count
│   └── 03_cellbender.sh          # Ambient RNA background removal
├── notebooks/
│   ├── 04_downstream_qc.Rmd      # Seurat/Scanpy based cellular QC filtering
│   ├── 05_recursive_clustering.Rmd# HPO and multi-level clustering logic
│   └── 06_cell_type_annotation.Rmd# Annotation using automatic tools and markers
└── results/                  
    ├── figures/                  # Generated plots (UMAPs, dotplots, etc.)
    └── objects/                  # Saved processed objects (.rds or .h5ad)

```

---

## Pipeline Workflow

### 1. Raw Data Quality Control (FastQC)

Initial quality assessment of raw sequence reads to evaluate base qualities, GC content, adapter contamination, and duplication levels.

* **Tools:** `FastQC` v0.11.9, `MultiQC` v1.11
* **Script:** `scripts/01_fastqc.sh`
* *Note: [Add your specific sequencing parameters or notes here]*

### 2. Alignment and Quantification (CellRanger)

Read alignment to the *Heterocephalus glaber* reference genome and generation of raw/filtered gene-expression matrices.

* **Tools:** `CellRanger` (v7.0.0 or specified version)
* **Script:** `scripts/02_cellranger_count.sh`
* *Note: Custom modifications made to the GTF file (e.g., handling of unannotated untranslated regions or lncRNAs) should be detailed here.*

### 3. Ambient RNA Removal (CellBender)

Elimination of background/ambient RNA contamination typical in single-nucleus preparations using a deep generative model.

* **Tools:** `CellBender` v0.3.0 (`remove-background`)
* **Script:** `scripts/03_cellbender.sh`
* *Note: [Specify expected-cells and total-droplets parameters used here]*

### 4. Downstream Quality Control & Filtering

Filtering of low-quality nuclei based on data-driven thresholds (e.g., minimum/maximum UMI counts, unique gene counts, and mitochondrial gene percentage bounds tailored to nuclear sequencing).

* **Tools:** `Scanpy` (Python)
* **Notebook:** `notebooks/04_downstream_qc.Rmd`

### 5. HPO & Recursive Clustering

Hyperparameter Optimization (HPO) for clustering resolution combined with an iterative, recursive clustering framework to separate major lineages (neurons, glia, vascular cells) down to distinct sub-populations.

* **Methodology:** [Detail your specific HPO or recursive split metrics here, e.g., Silhouette score, ROGUE, or resolution sweep]
* **Notebook:** `notebooks/05_recursive_clustering.Rmd`

### 6. Cell-Type Annotation

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
