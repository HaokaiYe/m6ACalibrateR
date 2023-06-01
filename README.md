# m6ACalibrateR

### Introduction
N6-methyladenosine (m6A) is the most prevalent and functionally significant mRNA modification in eukaryotes. However, discrepancies in m6A maps between studies have prompted concerns regarding the reliability of their biological validity, primarily attributed to non-specific antibody enrichment during immunoprecipitation (IP), which leads to false positives. *m6ACalibrateR*, a novel machine learning-based computational method, is designed to address this discrepancies by calibrating transcriptome-wide m6A maps. By integrating *genomic features*, we identify and eliminate non-specific antibody enrichment-induced false positives in *MeRIP-seq*, generating a high-accuracy m6A epitranscriptome map.

### Installation
To install m6ACalibrateR from Github, please use the following command in R console.
``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("HaokaiYe/m6ACalibrateR")
```

### Usage
First, load the package into R.
``` r
library(m6ACalibrateR)
```

Example 1: Using a gene annotation GFF/GTF file for transcript annotation
``` r
# Load example data
x <- readRDS(system.file("extdata", "peaks.rds", package = "m6ACalibrateR"))

# Specify gtf file path
gtf_file_path <- system.file("extdata", "annotation.gtf", package = "m6ACalibrateR")

# Calibrate m6A maps
calibrated_m6A <- m6ACalibrate(x, gff = gtf_file_path, genome = "hg38", anti_type = "ensemble")
calibrated_m6A
```


Example 2: Using a TxDb object for transcript annotation and a BSgenome object for the reference genome
``` r
# Load example data
x <- readRDS(system.file("extdata", "peaks.rds", package = "m6ACalibrateR"))

# Load TxDb and BSgenome
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

# Calibrate m6A maps
calibrated_m6A <- m6ACalibrate(x, txdb = txdb, genome = bsgenome, anti_type = "ensemble")
calibrated_m6A
```


Example 3: Changing antibody type and false positive threshold
``` r
# Load example data
x <- readRDS(system.file("extdata", "peaks.rds", package = "m6ACalibrateR"))

# Load TxDb and BSgenome
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

# Calibrate m6A maps
calibrated_m6A <- m6ACalibrate(x, txdb = txdb, genome = bsgenome, anti_type = "Abcam", FP_threshold = 0.4)
calibrated_m6A
```