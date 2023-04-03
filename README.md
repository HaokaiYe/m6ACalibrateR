# m6ACalibrateR

### Introduction
N6-methyladenosine (m6A) is the most prevalent and functionally significant mRNA modification in eukaryotes. However, discrepancies in m6A maps between studies have prompted concerns regarding the reliability of their biological validity, primarily attributed to non-specific antibody enrichment during immunoprecipitation (IP), which leads to false positives. *m6ACalibrateR*, a novel machine learning-based computational method, is designed to address this discrepancies by calibrating transcriptome-wide m6A maps. By integrating *genomic features*, we identify and eliminate non-specific antibody enrichment-induced false positives in *MeRIP-seq*, generating a high-accuracy m6A epitranscriptome map.

### Installation
To install m6ACalibrateR from Github, please use the following command in R console.
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("HaokaiYe/m6ACalibrateR")
```

### Usage
Example m6A coordinates can be found in inst/extdata:
```
x <- readRDS(system.file("extdata", "peaks.rds", package = "m6ACalibrateR"))
```
It is recommended to use the function **encGeo** to generate the encoding. Different encodings can be selected by the parameter *type*:
```
library(m6ACalibrateR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

# Calibrate m6A maps with default parameters
calibrated_m6A <- m6ACalibrate(x, txdb, bsgenome, "ensemble")
calibrated_m6A
```
