#' Calibrate putative m6A maps obtained from peak-calling.
#'
#' @description \code{m6ACalibrate} is specifically designed to refine putative m6A maps identified from MeRIP-Seq data into calibrated, high-confidence m6A maps.
#'
#' @details \code{m6ACalibrate} operates as an advanced site calibration tool, specifically tailored for refining and scoring putative m6A peaks obtained from MeRIP-Seq data analysis.
#'
#' Designed to seamlessly integrate into existing m6A profiling pipelines, this function acts as a critical subsequent step following the initial peak calling process.
#'
#' By accurately extracting and calibrating motifs from m6A peaks, and utilizing a selected model to meticulously identify and filter out false positives, \code{m6ACalibrate} significantly enhances the precision of m6A site identification.
#'
#' The ability to choose from different filtering modes (\code{"directFilter"}, \code{"strictFilter"}, and \code{"flankExclusion"}) allows for customized calibration that can adapt to various analysis requirements, thus enabling a more tailored approach to site validation.
#'
#'
#' @param x A \code{\link{GRanges}} or \code{\link{GRangesList}} object for the genomic location of m6A sites.
#'
#' @param txdb A \code{\link{TxDb}} object for the transcript annotation.
#' The \code{TxDb} can be obtained from either Bioconductor or from the GFF/GTF files using the function \code{\link{makeTxDbFromGFF}}.
#'
#' @param genome A \code{character} or a \code{\link{BSgenome}} for the reference genome.
#' The character should be the UCSC genome name which is acceptable by \code{\link{getBSgenome}} or/and \code{\link{makeTxDbFromUCSC}}; example: \code{"hg38"}.
#'
#' @param gff Optional, a \code{character} specifying the path to a gene annotation GFF/GTF file.
#' It is applied when the \code{TxDb} object is not available; default \code{= NULL}.
#'
#' @param model_type A \code{character} specifying the model trained on data from three commonly used m6A-specific antibodies. Can be one of \code{c("Abcam", "NEB", "SYSY", "ensemble")} with the following implications:
#' \itemize{
#'   \item \code{"Abcam"}: Abcam (ab190886).
#'   \item \code{"NEB"}: NEB (E1610S).
#'   \item \code{"SYSY"}: Synaptic Systems (202003).
#'   \item \code{"ensemble"}: The mean of the predictions from the Abcam, NEB, and SYSY models.
#'}; default \code{= "ensemble"}.
#'
#' @param FP_threshold A \code{numeric} between 0 and 1 specifying the cutoff for predicting false positives.
#' Applicable only with \code{model_type = "ensemble"}; defaults to training set proportions for other models.
#' Smaller values indicate stricter filtering; default \code{= 0.5}.
#'
#' @param filter_mode A \code{character} specifying the method used to filter false positives.
#' Can be one of \code{c("directFilter", "strictFilter", "flankExclusion")} with the following implications:
#' \itemize{
#'   \item \code{"directFilter"}: Directly removes all m6A peaks that overlap with predicted false positives.
#'   \item \code{"strictFilter"}: Removes all m6A peaks overlapping with false positives and additionally excludes peaks that do not contain any of the specified motifs (e.g., DRACH motifs). The specified motif is determined by the \code{motif_sequence} parameter.
#'   \item \code{"flankExclusion"}: Excludes a specified flank width around each false positive without directly removing the peaks that overlap with false positives. The flank width is determined by the \code{flank_width} parameter.
#'}; default \code{= "strictFilter"}.
#'
#' @param flank_width A \code{numeric} specifying the width of the flanking region to be excluded around false positives.
#' Only applicable when \code{filter_mode = "flankExclusion"}; default \code{= 50}.
#'
#' @param motif_sequence A \code{character} specifying the motif sequence used in predicting false positives; default \code{= "DRACH"}.
#'
#' @param save_motif_prob A \code{logical} indicating whether to include predicted motif probabilities in the output; default \code{= FALSE}.
#'
#'
#' @importFrom randomForest randomForest
#' @importFrom magrittr "%>%"
#' @import GenomicFeatures
#' @import GenomicRanges
#' @import BSgenome
#' @import IRanges
#' @import AnnotationDbi
#'
#'
#' @return Depending on the input and parameters, this function returns:
#' \itemize{
#'   \item If \code{save_motif_prob = FALSE}, a \code{\link{GRanges}} or \code{\link{GRangesList}} object containing the calibrated m6A locations, matching the format of the input 'x'.
#'   \item If \code{save_motif_prob = TRUE}, a list object containing two elements:
#'     \describe{
#'       \item{motif}{A \code{\link{GRanges}} object containing the predicted motifs. Each motif includes its false positive probabilities (\code{pred_FP_prob}) and annotations indicating from which region of 'x' it was derived (\code{sourceHits}).}
#'       \item{calibrated_m6A}{A \code{\link{GRanges}} or \code{\link{GRangesList}} object containing the calibrated m6A sites, matching the format of the input 'x'.}
#'     }
#' }
#'
#'
#' @examples
#'
#' library(m6ACalibrateR)
#'
#' # Example 1: Using a GFF/GTF file for transcript annotation
#' # Load example m6A peak data
#' peaks <- readRDS(system.file("extdata", "peaks.rds", package = "m6ACalibrateR"))
#'
#' # Path to the gene annotation file
#' gtf_path <- system.file("extdata", "annotation.gtf", package = "m6ACalibrateR")
#'
#' # Perform calibration using an ensemble model
#' calibrated_peaks <- m6ACalibrate(peaks, gff = gtf_path, genome = "hg38", model_type = "ensemble")
#' print(calibrated_peaks)
#'
#'
#' # Example 2: Using a TxDb object for transcript annotation and a BSgenome object for the reference genome
#' # Reload example m6A peak data (assuming 'peaks' variable already loaded)
#'
#' # Load necessary libraries for TxDb and BSgenome
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Define TxDb and BSgenome objects
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#'
#' # Calibrate m6A maps using a specific false positive threshold
#' calibrated_peaks_ensemble <- m6ACalibrate(peaks, txdb = txdb, genome = bsgenome, model_type = "ensemble", FP_threshold = 0.6)
#' print(calibrated_peaks_ensemble)
#'
#'
#' # Example 3: Customizing antibody type and filter mode
#' # Reload example m6A peak data (assuming 'peaks' variable already loaded)
#'
#' # Calibrate using a specific antibody (e.g., "Abcam") and a filter mode to exclude flanking regions around false positives
#' calibrated_peaks_abcam <- m6ACalibrate(peaks, txdb = txdb, genome = bsgenome, model_type = "Abcam", filter_mode = "flankExclusion", flank_width = 100)
#' print(calibrated_peaks_abcam)
#'
#'
#' @export

m6ACalibrate <- function(x,
                         txdb = NULL,
                         genome = NULL,
                         gff = NULL,
                         model_type = c("ensemble", "Abcam", "NEB", "SYSY"),
                         FP_threshold = 0.5,
                         filter_mode = c("strictFilter", "directFilter", "flankExclusion"),
                         flank_width = 50,
                         motif_sequence = "DRACH",
                         save_motif_prob = FALSE) {

  # Validate inputs
  model_type <- match.arg(model_type)
  filter_mode <- match.arg(filter_mode)
  stopifnot(flank_width > 0)
  stopifnot(is(x, "GRangesList") || is(x, "GRanges"))
  stopifnot(is.numeric(FP_threshold) && FP_threshold >= 0 && FP_threshold <= 1)

  # Prepare TxDb and genome
  txdb <- prepareTxDb(txdb, gff, genome)
  genome <- prepareGenome(genome, txdb)

  # Process motifs and genomic features
  message("## Extracting motifs and genomic features ...")
  motif_within_x <- sampleSequence(motif_sequence, x + 2, genome) - 2
  geo_features <- extractGeoFeature(motif_within_x, txdb)

  # Predict the confidences of m6A sites and filter
  motif_within_x <- predictMotifConfidence(motif_within_x, geo_features, model_type, FP_threshold)
  x_calibrated <- filterFP(x, motif_within_x, filter_mode, flank_width)

  # Return result
  if (save_motif_prob) {
    list(motif = motif_within_x, calibrated_m6A = x_calibrated)
  } else {
    x_calibrated
  }
}
