
extractGeoFeature <- function(x, txdb){

  ################################################################################
  ####                                exons                                   ####
  ################################################################################
  message("#### Extracting properties of exons ...")
  exons <- exons(txdb)
  log2_length_exons <- log2(extractRegionLength(x, exons) + 1)
  relativePOS_exons <- extractRegionRelativePosition(x, exons)
  rm(exons)

  ################################################################################
  ####                              exonicCDS                                 ####
  ################################################################################
  message("#### Extracting properties of exonic CDS ...")
  exonicCDS <- cdsBy(txdb, by = "tx")
  log2_length_exonicCDS <- log2(extractRegionLength(x, exonicCDS) + 1)
  relativePOS_exonicCDS <- extractRegionRelativePosition(x, exonicCDS)
  rm(exonicCDS)

  ################################################################################
  ####                            exonicGenes                                 ####
  ################################################################################
  message("#### Extracting properties of exonic genes ...")
  exonicGenes <- exonsBy(txdb, by = "gene")
  log2_length_exonicGenes <- log2(extractRegionLength(x, exonicGenes) + 1)
  relativePOS_exonicGenes <- extractRegionRelativePosition(x, exonicGenes)

  ################################################################################
  ####                          exonicThreePrimeUTR                           ####
  ################################################################################
  message("#### Extracting properties of exonic 3'UTR ...")
  exonicThreePrimeUTR <- threeUTRsByTranscript(txdb)
  log2_length_exonicThreePrimeUTR <- log2(extractRegionLength(x, exonicThreePrimeUTR) + 1)
  rm(exonicThreePrimeUTR)

  ################################################################################
  ####                          exonicTranscripts                             ####
  ################################################################################
  message("#### Extracting properties of exonic transcripts ...")
  exonicTranscripts <- exonsBy(txdb, by = "tx")
  log2_length_exonicTranscripts <- log2(extractRegionLength(x, exonicTranscripts) + 1)
  relativePOS_exonicTranscripts <- extractRegionRelativePosition(x, exonicTranscripts)
  rm(exonicTranscripts)

  ################################################################################
  ####                               fullCDS                                  ####
  ################################################################################
  message("#### Extracting properties of full CDS ...")
  fullCDS <- unlist(range(cdsBy(txdb, by = "tx")))
  log2_length_fullCDS <- log2(extractRegionLength(x, fullCDS) + 1)
  relativePOS_fullCDS <- extractRegionRelativePosition(x, fullCDS)
  rm(fullCDS)

  ################################################################################
  ####                              fullGenes                                 ####
  ################################################################################
  message("#### Extracting properties on full genes ...")
  fullGenes <- unlist(genes(txdb, single.strand.genes.only = FALSE))
  log2_length_fullGenes <- log2(extractRegionLength(x, fullGenes) + 1)
  relativePOS_fullGenes <- extractRegionRelativePosition(x, fullGenes)
  rm(fullGenes)

  ################################################################################
  ####                           fullThreePrimeUTR                            ####
  ################################################################################
  message("#### Extracting properties of full 3'UTR ...")
  fullThreePrimeUTR <- unlist(range(threeUTRsByTranscript(txdb)))
  log2_length_fullThreePrimeUTR <- log2(extractRegionLength(x, fullThreePrimeUTR) + 1)
  rm(fullThreePrimeUTR)

  ################################################################################
  ####                            fullTranscripts                             ####
  ################################################################################
  message("#### Extracting properties of full transcripts ...")
  fullTranscripts <- transcripts(txdb)
  log2_length_fullTranscripts <- log2(extractRegionLength(x, fullTranscripts) + 1)
  relativePOS_fullTranscripts <- extractRegionRelativePosition(x, fullTranscripts)
  rm(fullTranscripts)

  ################################################################################
  ####                            GeneExonNumber                              ####
  ################################################################################
  message("#### Extracting properties meta-genomic statistics ...")
  log2_GeneExonNumber <- log2(extractRegionProperty(x,
                                                    property = elementNROWS(exonicGenes),
                                                    region = range(exonicGenes),
                                                    nomapValue = "0") + 1)
  rm(exonicGenes)

  ################################################################################
  ####                            TxIsoformNumber                             ####
  ################################################################################
  txbg <- transcriptsBy(txdb, by = "gene")
  log2_TxIsoformNumber <- log2(extractRegionProperty(x,
                                                     property = elementNROWS(txbg),
                                                     region = range(txbg),
                                                     nomapValue = "0") + 1)
  rm(txbg)

  ################################################################################
  ####                            MetaTxTopology                              ####
  ################################################################################
  MetaTxTopology <- suppressWarnings(topologyOnTranscripts(x, txdb))



  X <- as.data.frame(cbind(log2_length_exons, relativePOS_exons,
             log2_length_exonicTranscripts, relativePOS_exonicTranscripts,
             log2_length_fullGenes, relativePOS_fullGenes,
             log2_length_exonicCDS, relativePOS_exonicCDS,
             log2_length_exonicGenes, relativePOS_exonicGenes,
             log2_length_fullTranscripts, relativePOS_fullTranscripts,
             log2_length_fullCDS, relativePOS_fullCDS,
             log2_GeneExonNumber,
             log2_TxIsoformNumber,
             log2_length_exonicThreePrimeUTR,
             log2_length_fullThreePrimeUTR,
             MetaTxTopology))

  X[is.na(X)] <- 0
  return(X)
}




################################################################################
####                        extractRegionLength                             ####
################################################################################
extractRegionLength <- function(x,
                                region = NULL,
                                ambiguityMethod = c("mean", "sum", "min", "max"),
                                maxgap = -1L,
                                minoverlap = 0L,
                                type = c("any", "start", "end", "within", "equal"),
                                nomapValue = c("NA", "0", "nearest"),
                                ignore.strand = FALSE) {
  stopifnot(is(x, "GRanges"))
  ambiguityMethod <- match.arg(ambiguityMethod)
  nomapValue <- match.arg(nomapValue)
  type <- match.arg(type)

  if (is.null(region)) {
    length_property <- width(x)
  } else if (is(region, "GRanges")) {
    length_property <- width(region)
  } else if (is(region, "GRangesList")) {
    region <- region[elementNROWS(region) != 0]
    length_property <- sum(width(region))
  } else{
    stop("`region` should be either `GRanges` or `GRangesList`")
  }

  if (is.null(region)) {
    return(length_property)
  } else{
    region_property <- extractRegionProperty(
      x = x,
      property = length_property,
      region = region,
      ambiguityMethod = ambiguityMethod,
      maxgap = maxgap,
      minoverlap = minoverlap,
      type = type,
      nomapValue = nomapValue,
      ignore.strand = ignore.strand
    )
    return(region_property)
  }
}


################################################################################
####                     extractRegionRelativePosition                      ####
################################################################################
extractRegionRelativePosition <- function(x,
                                          region = NULL,
                                          ambiguityMethod = c("mean", "sum", "min", "max"),
                                          nomapValue = c("NA", "0"),
                                          ignore.strand = FALSE) {
  stopifnot(is(x, "GRanges"))
  ambiguityMethod <- match.arg(ambiguityMethod)
  nomapValue <- match.arg(nomapValue)
  nomapValue <- eval(parse(text = nomapValue))

  if (is.null(region)) {
    rrp_property <- rep(nomapValue, length(x))
  } else if (is(region, "GRanges") | is(region, "GRangesList")) {
    if (is(region, "GRanges")) {
      region_grl <- split(region, seq_along(region))
    } else{
      region <- region[elementNROWS(region) != 0]
      more_strand_region <-
        which(elementNROWS(runValue(strand(region))) == 1)
      region <- grl_resolve_multi_strand(region)
      names(region) <- seq_along(region)
      region_grl <- region
    }
    rrp_property <- rep(nomapValue, length(x))
    map2tx <-
      mapToTranscripts(x, region_grl, ignore.strand = ignore.strand)
    relpos <-
      start(map2tx) / sum(width(region_grl))[map2tx$transcriptsHits]
    weighted_relpos <-
      tapply(relpos, map2tx$xHits, eval(parse(text = ambiguityMethod[1])))
    rm(relpos)
    rrp_property[as.numeric(names(weighted_relpos))] <-
      weighted_relpos
    rm(map2tx, weighted_relpos)
  } else{
    stop("`region` should be either `GRanges` or `GRangesList`")
  }
  return(rrp_property)
}


################################################################################
####                        extractRegionProperty                           ####
################################################################################
extractRegionProperty <- function(x,
         region,
         property,
         ambiguityMethod = c("auto", "mean", "sum", "min", "max"),
         maxgap = 0L,
         minoverlap = 0L,
         type = c("any", "start", "end", "within", "equal"),
         nomapValue = c("NA", "0", "FALSE", "nearest"),
         ignore.strand = FALSE) {
  stopifnot(is(x, "GRanges"))
  stopifnot(is(region, "GRanges") | is(region, "GRangesList"))
  stopifnot(length(property) == length(region))
  ambiguityMethod <- match.arg(ambiguityMethod)
  type <- match.arg(type)
  nomapValue <- match.arg(nomapValue)
  return_property <- rep(NA, length(x))

  fol <-
    findOverlaps(x,
                 region,
                 maxgap = maxgap,
                 minoverlap = minoverlap,
                 type = type)
  property_mapped <- property[subjectHits(fol)]

  if (!is.null(property_mapped)) {
    if (ambiguityMethod == "auto") {
      if (is.logical(property)) {
        weighted_properties <- tapply(property_mapped, queryHits(fol), any)
      } else{
        weighted_properties <- tapply(property_mapped, queryHits(fol), mean)
      }
    } else{
      weighted_properties <-
        tapply(property_mapped, queryHits(fol), eval(parse(text = ambiguityMethod)))
    }
    return_property[as.numeric(names(weighted_properties))] <-
      weighted_properties
  }

  if (anyNA(return_property)) {
    if (nomapValue == "nearest") {
      if (!is.null(property_mapped)) {
        if (is(region, "GRanges")) {
          nomap_indx <- which(is.na(return_property))
          return_property[nomap_indx] <-
            property[nearest(x[nomap_indx], region)]
          return_property[is.na(return_property)] <- 0
        } else{
          region_ranges <- range_unique_grl(region)
          nomap_indx <- which(is.na(return_property))
          return_property[nomap_indx] <-
            property[nearest(x[nomap_indx], unlist(region_ranges))]
          rm(region_ranges, nomap_indx)
          return_property[is.na(return_property)] <- 0
        }
      }
    } else if (nomapValue == "0") {
      return_property[is.na(return_property)] <- 0
    } else if (nomapValue == "FALSE") {
      return_property[is.na(return_property)] <- FALSE
    }
  }
  return(return_property)
}

################################################################################
####                        grl_resolve_multi_strand                        ####
################################################################################
grl_resolve_multi_strand <- function(grl) {
  indx_multi <- elementNROWS(runValue(strand(grl))) > 1
  if (any(indx_multi)) {
    gr_ambiguous <- unlist(grl[indx_multi])
    names(gr_ambiguous) <-
      paste0(names(gr_ambiguous), strand(gr_ambiguous))
    gr_resolved <- split(gr_ambiguous, names(gr_ambiguous))
    rm(gr_ambiguous)
    return(c(grl[!indx_multi], gr_resolved))
  } else{
    return(grl)
  }
}

################################################################################
####                         topologyOnTranscripts                          ####
################################################################################
topologyOnTranscripts <- function(x,
                                  txdb,
                                  region_weights = c(1/3,1/3,1/3),
                                  ambiguityMethod = c("mean", "sum", "min", "max"),
                                  ignore.strand=FALSE){

  u5bytx <- (fiveUTRsByTranscript(txdb))
  topology <- extractRegionRelativePosition(x,
                                            u5bytx,
                                            ambiguityMethod=ambiguityMethod,
                                            nomapValue="NA",
                                            ignore.strand=ignore.strand)*region_weights[1]
  rm(u5bytx)
  cdsbytx <- (cdsBy(txdb, by = "tx"))
  cdsrps <- extractRegionRelativePosition(x,
                                          cdsbytx,
                                          ambiguityMethod=ambiguityMethod,
                                          nomapValue="NA",
                                          ignore.strand=ignore.strand)
  rm(cdsbytx)
  indx <- !is.na(cdsrps)
  topology[indx] <- cdsrps[indx]*region_weights[2] + region_weights[1]
  rm(indx,cdsrps)

  u3bytx <- (threeUTRsByTranscript(txdb))
  u3rps <- extractRegionRelativePosition(x,
                                         u3bytx,
                                         ambiguityMethod=ambiguityMethod,
                                         nomapValue="NA",
                                         ignore.strand=ignore.strand)
  rm(u3bytx)
  indx <- !is.na(u3rps)
  topology[indx] <- u3rps[indx]*region_weights[3] + region_weights[2] + region_weights[1]
  rm(indx,u3rps)

  return(topology)
}


################################################################################
####                           sampleSequence                               ####
################################################################################

sampleSequence <- function(motif, region, sequence, fixed = FALSE){
  require(BSgenome)
  require(GenomicFeatures)
  stopifnot(is(region, "GRangesList")|is(region, "GRanges"))

  isGrl <- is(region, "GRangesList")
  if(isGrl) {
    region <- unlist(region)
    region$listId <- sub("\\..*$", "", names(region))
  }
  #region <- reduce(region)

  region_dnass <- getSeq(x=sequence,
                         names=seqnames(region),
                         start=start(region),
                         end=end(region),
                         strand=strand(region),
                         as.character=FALSE)

  indx <- paste0("reg_", seq_along(region))
  regions_GRL <- split(region, indx)
  regions_GRL <- regions_GRL[indx]
  rm(indx)
  vmp <- vmatchPattern(motif, region_dnass, fixed = fixed)
  rm(region_dnass)
  vmp_gr <- GRanges(seqnames = rep(names(regions_GRL), elementNROWS(vmp)), ranges = unlist(vmp))
  rm(vmp)
  motif_on_regions <- mapFromTranscripts(vmp_gr,regions_GRL)
  rm(vmp_gr, regions_GRL)

  transcriptsHits <- motif_on_regions$transcriptsHits
  mcols(motif_on_regions) <- NULL
  if (isGrl) {
    motif_on_regions$sourceHits <- region$listId[transcriptsHits]
    names(motif_on_regions) <- motif_on_regions$sourceHits
  } else {
    motif_on_regions$sourceHits <- transcriptsHits
  }
  rm(transcriptsHits)

  seqlengths(motif_on_regions) <- seqlengths(region)
  return(motif_on_regions)
}




