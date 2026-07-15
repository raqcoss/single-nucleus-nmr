

library(devtools)
library(BiocManager)
# Install SoupX from CRAN
install.packages('SoupX')
# To use latelest development version of SoupX
devtools::install_github("constantAmateur/SoupX",ref='devel')

library(Seurat)
library(SoupX)
library(SeuratDisk)
Convert("/home/raquelcr/scanpy/axolotl/Axolotl_midbrain_s2_raw.h5ad", dest = "h5seurat", overwrite = TRUE)
sc = load10X('/home/raquelcr/scanpy/axolotl/Axolotl_midbrain_s2_raw.h5ad')
sc = autoEstCont(sc)
out = adjustCounts(sc)