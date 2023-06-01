#' Calibrate putative m6A maps obtained from peak-calling.
#'
#' @description This function takes in putative m6A maps as input and generates calibrated and high-confidence m6A maps as output.
#'
#' @details \code{m6ACalibrate} extracts DRACH motifs from m6A peaks within exons and utilizes a selected model to identify false positives within these motifs.
#' It removes all m6A peaks overlapping with false positives, while retaining only those containing DRACH motifs.
#' The function outputs a high-precision m6A map after this calibration process.
#'
#' @param x A \code{\link{GRanges}} object for the genomic location of m6A sites.
#' @param txdb A \code{\link{TxDb}} object for the transcript annotation.
#'
#' The \code{TxDb} can be obtained from either Bioconductor or from the GFF/GTF files using the function \code{\link{makeTxDbFromGFF}}.
#'
#' @param genome A \code{character} or a \code{\link{BSgenome}} for the reference genome.
#'
#' The character should be the UCSC genome name which is acceptable by \code{\link{getBSgenome}} or/and \code{\link{makeTxDbFromUCSC}}; example: \code{"hg38"}.
#'
#' @param gff optional, a \code{character} which specifies the directory toward a gene annotation GFF/GTF file, it is applied when the \code{TxDb} object is not available; default \code{= NULL}.
#'
#' @param anti_type A \code{character} specifying the type of m6A-specific antibody, can be one of \code{c("Abcam", "NEB", "SYSY", "ensemble")}; default \code{= "ensemble"}.
#' @param FP_threshold A \code{numeric} value specifying the probability cutoff in predicting false positives at each m6A site.
#'
#' It should be a numeric value between 0 and 1, where smaller values indicate stricter filtering; default \code{= 0.5}.
#'
#' @importFrom randomForest randomForest
#' @importFrom magrittr "%>%"
#' @import GenomicFeatures
#' @import GenomicRanges
#' @import BSgenome
#' @import IRanges
#' @import AnnotationDbi
#'
#' @return A \code{GRanges} object containing the calibrated m6A sites.
#'
#' @examples
#'
#' library(m6ACalibrateR)
#'
#' # Example 1: Using a gene annotation GFF/GTF file for transcript annotation
#' # Load example data
#' x <- readRDS(system.file("extdata", "peaks.rds", package = "m6ACalibrateR"))
#'
#' # Specify gtf file path
#' gtf_file_path <- system.file("extdata", "annotation.gtf", package = "m6ACalibrateR")
#'
#' # Calibrate m6A maps
#' calibrated_m6A <- m6ACalibrate(x, gff = gtf_file_path, genome = "hg38", anti_type = "ensemble")
#' calibrated_m6A
#'
#'
#' # Example 2: Using a TxDb object for transcript annotation and a BSgenome object for the reference genome
#' # Load example data
#' x <- readRDS(system.file("extdata", "peaks.rds", package = "m6ACalibrateR"))
#'
#' # Load TxDb and BSgenome
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#'
#' # Calibrate m6A maps
#' calibrated_m6A <- m6ACalibrate(x, txdb = txdb, genome = bsgenome, anti_type = "ensemble")
#' calibrated_m6A
#'
#'
#' # Example 3: Changing antibody type and false positive threshold
#' # Load example data
#' x <- readRDS(system.file("extdata", "peaks.rds", package = "m6ACalibrateR"))
#'
#' # Load TxDb and BSgenome
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#'
#' # Calibrate m6A maps
#' calibrated_m6A <- m6ACalibrate(x, txdb = txdb, genome = bsgenome, anti_type = "Abcam", FP_threshold = 0.4)
#' calibrated_m6A
#'
#
#'
#' @export

m6ACalibrate <- function(x,
                         txdb = NULL,
                         genome = NULL,
                         gff = NULL,
                         anti_type = c("Abcam", "NEB", "SYSY", "ensemble"),
                         FP_threshold = 0.5) {

  # Check input validity
  if(!is(x, "GRanges")) {
    stop("Argument 'x' must be a GRanges object.", call. = FALSE)
  }

  if (missing(anti_type)) {
    anti_type <- "ensemble"
  } else if (length(anti_type) > 1) {
    stop("Only one 'anti_type' value allowed.")
  } else if(!is.character(anti_type) | !anti_type %in% c("Abcam", "NEB", "SYSY", "ensemble")) {
    stop("Argument 'anti_type' must be a character and must be one of 'Abcam', 'NEB', 'SYSY', or 'ensemble'.", call. = FALSE)
  }

  if (!is.numeric(FP_threshold) || FP_threshold < 0 || FP_threshold > 1) {
    stop("Argument 'FP_threshold' must be a numeric value between 0 and 1.", call. = FALSE)
  }

  # Prepare transcript annotation
  if (is.null(gff) & is.null(txdb) & is.null(genome)){
    stop("Require one of the argument in txdb, gff, and genome for transcript annotation.")
  }

  if(!is.null(txdb) & !is(txdb, "TxDb")) {
    stop("Argument 'txdb' must be a TxDb object.", call. = FALSE)
  }

  if (is.null(txdb) & is.null(gff) & is.character(genome)) {
    txdb <- makeTxDbFromUCSC(genome)
  }

  if (is.null(txdb) & !is.null(gff)) {
    txdb <- makeTxDbFromGFF(gff)
  }

  # Prepare reference genome
  if (is.character(genome)){
    genome <- getBSgenome(genome)
  } else if (is.null(genome)) {
    if(!is.null(txdb)){
      organism <- organism(txdb)
      genome <- getBSgenome(organism)
    } else {
      stop("Neither genome nor txdb provided. Please provide one of these or both.")
    }
  } else if(!is(genome, "BSgenome")) {
    stop("Argument 'genome' must be a character or a BSgenome object.", call. = FALSE)
  }


  message("## Loading random forest models ...")

  rf_Abcam <- readRDS(system.file("data", "Abcam.rds", package = "m6ACalibrateR"))
  rf_NEB <- readRDS(system.file("data", "NEB.rds", package = "m6ACalibrateR"))
  rf_SYSY <- readRDS(system.file("data", "SYSY.rds", package = "m6ACalibrateR"))

  # Extract all DRACH motif on exons
  exbtx <- exonsBy(txdb, by = "tx")
  motif_all <- sort(sampleSequence("DRACH", exbtx, genome) - 2)

  # Find DRACH motifs overlapping with input m6A sites
  x_motif <- subsetByOverlaps(motif_all, x)

  message("## Extracting genomic features ...")

  gfeatures <- feature_extraction(x_motif, txdb)

  message("## Predicting the confidences of m6A sites ...")

  if(anti_type == "Abcam") {
    rf_pred <- predict(rf_Abcam, newdata = gfeatures)
    filter_index <- which(rf_pred == 0)
  } else if(anti_type == "NEB") {
    rf_pred <- predict(rf_NEB, newdata = gfeatures)
    filter_index <- which(rf_pred == 0)
  } else if(anti_type == "SYSY") {
    rf_pred <- predict(rf_SYSY, newdata = gfeatures)
    filter_index <- which(rf_pred == 0)
  } else {
    rf_Abcam_pred <- predict(rf_Abcam, newdata = gfeatures, type = "prob")
    rf_NEB_pred <- predict(rf_NEB, newdata = gfeatures, type = "prob")
    rf_SYSY_pred <- predict(rf_SYSY, newdata = gfeatures, type = "prob")
    predp0 = data.frame(abcam = rf_Abcam_pred[,1], neb = rf_NEB_pred[,1], sysy = rf_SYSY_pred[,1])
    filter_index <- which(rowMeans(predp0) >= FP_threshold)
  }

  message("## Filtering out low-confidence m6A sites ...")

  motif_x <- subsetByOverlaps(x, motif_all)
  filter_x_motif <- x_motif[filter_index]
  filter_x_index <- unique(queryHits(findOverlaps(motif_x, filter_x_motif)))

  # Check if any m6A sites were filtered out
  if(length(filter_x_index) > 0) {
    x_retain <- motif_x[-filter_x_index]
  } else {
    x_retain <- motif_x
  }

  # Count the number of filtered m6A sites
  num_filtered <- length(x) - length(x_retain)
  if(num_filtered > 0) {
    message(paste0("## ", num_filtered, " m6A sites were filtered out due to low-confidence predictions."))
  } else {
    message("## All m6A sites were high-confidence maps.")
  }

  return(x_retain)
}
