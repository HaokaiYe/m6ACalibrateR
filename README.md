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

Example 1: Using a GFF/GTF file for transcript annotation
``` r
# Load example m6A peak data
peaks <- readRDS(system.file("extdata", "peaks.rds", package = "m6ACalibrateR"))

# Path to the gene annotation file
gtf_path <- system.file("extdata", "annotation.gtf", package = "m6ACalibrateR")

# Perform calibration using an ensemble model
calibrated_peaks <- m6ACalibrate(peaks,
                                 gff = gtf_path,
                                 genome = "hg38",
                                 model_type = "ensemble")
print(calibrated_peaks)
```


Example 2: Using a TxDb object for transcript annotation and a BSgenome object for the reference genome
``` r
# Reload example m6A peak data (assuming 'peaks' variable already loaded)

# Load necessary libraries for TxDb and BSgenome
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

# Define TxDb and BSgenome objects
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

# Calibrate m6A maps using a specific false positive threshold
calibrated_peaks_ensemble <- m6ACalibrate(peaks,
                                          txdb = txdb,
                                          genome = bsgenome,
                                          model_type = "ensemble",
                                          FP_threshold = 0.6)
print(calibrated_peaks_ensemble)
```


Example 3: Customizing antibody type and filter mode
``` r
# Reload example m6A peak data (assuming 'peaks' variable already loaded)

# Calibrate using a specific antibody (e.g., "Abcam") and a filter mode to exclude flanking regions around false positives
calibrated_peaks_abcam <- m6ACalibrate(peaks,
                                       txdb = txdb,
                                       genome = bsgenome,
                                       model_type = "Abcam",
                                       filter_mode = "flankExclusion",
                                       flank_width = 100)
print(calibrated_peaks_abcam)
```
