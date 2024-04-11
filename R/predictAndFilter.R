
predictMotifConfidence <- function(motif, features, model_type, FP_threshold) {
  # For non-ensemble cases, only load the specified model
  if (model_type != "ensemble") {
    rf_model <- readRDS(system.file("data", paste0(model_type, ".rds"), package = "m6ACalibrateR"))
    message("## Predicting m6A site confidences using the ", model_type, " model...")

    # Predicting class and probability
    motif$pred_class <- predict(rf_model, newdata = features)
    motif$pred_FP_prob <- predict(rf_model,
                                           newdata = features,
                                           type = "prob")[, 1]
  } else {
    # For the ensemble case, load all models and average the prediction probabilities
    model_types <- c("Abcam", "NEB", "SYSY")
    rf_models <- lapply(model_types, function(type)
      readRDS(system.file("data", paste0(type, ".rds"), package = "m6ACalibrateR"))
    )
    message("## Predicting m6A site confidences using the ensemble method...")

    # Averaging prediction probabilities across all models
    pred_FP_probs <- sapply(rf_models, function(model)
      predict(model, newdata = features, type = "prob")[, 1]
    )

    # Determining class and averaging probabilities
    motif$pred_class <- 1 - as.integer(rowMeans(pred_FP_probs) >= FP_threshold)
    motif$pred_FP_prob <- rowMeans(pred_FP_probs)
  }
  return(motif)
}

filterFP <- function(x, motif, filter_mode, flank_width) {
  message("## Filtering out low-confidence m6A sites ...")

  # Identify motifs predicted as false positives
  FP_motif <- motif[motif$pred_class == 0]

  # Filter modes
  if (filter_mode == "directFilter") {
    return(subsetByOverlaps(x, FP_motif, invert = TRUE))

  } else if (filter_mode == "strictFilter") {
    x_overlap_motif <- subsetByOverlaps(x, motif)
    return(subsetByOverlaps(x_overlap_motif, FP_motif, invert = TRUE))

  } else if (filter_mode == "flankExclusion") {
    FP_flank <- adjustFlank(FP_motif, flank_width)

    if (is(x, "GRangesList")) {
      x_unlist <- unlist(x)
      disjoin_region <- disjoin(c(x_unlist, FP_flank))
      x_calibrated <- subsetByOverlaps(disjoin_region, FP_flank, invert = TRUE)

      ov <- findOverlaps(x_calibrated, x)
      listId <- factor(names(x)[subjectHits(ov)], levels = names(x))
      x_calibrated <- split(x_calibrated, listId)

      return(x_calibrated[elementNROWS(x_calibrated) > 0])

    } else {
      disjoin_region <- disjoin(c(x, FP_flank))
      return(subsetByOverlaps(disjoin_region, FP_flank, invert = TRUE))
    }
  }
}

adjustFlank <- function(motif, width) {
  # Adjust motif start positions to not go below 1 after flank subtraction
  flanked <- motif + width
  start(flanked) <- pmax(start(flanked), 1)
  return(flanked)
}

