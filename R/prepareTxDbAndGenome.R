
prepareTxDb <- function(txdb, gff, genome) {
  # Directly validate the txdb object
  if (!is.null(txdb)) {
    if (!is(txdb, "TxDb")) {
      stop("The 'txdb' argument must be a 'TxDb' object.", call. = FALSE)
    }
    return(txdb)
  }

  # Create a TxDb object using a GFF file
  if (!is.null(gff)) {
    return(makeTxDbFromGFF(gff))
  }

  # Create a TxDb object using a UCSC genome name
  if (is.character(genome)) {
    return(makeTxDbFromUCSC(genome))
  }

  # If no valid input is provided, throw an error
  stop("One of 'txdb', 'gff', or 'genome' must be provided.", call. = FALSE)
}

prepareGenome <- function(genome, txdb) {
  # If genome is already a BSgenome object, return it directly
  if (is(genome, "BSgenome")) {
    return(genome)
  }

  # If genome is a character string, attempt to directly load the corresponding BSgenome object
  if (is.character(genome)) {
    tryCatch({
      return(getBSgenome(genome))
    }, error = function(e) {
      stop("Failed to load BSgenome for ", genome, ": ", e$message, call. = FALSE)
    })
  }

  # If genome is NULL and txdb is provided, attempt to infer the organism from txdb and load the BSgenome
  if (is.null(genome) && !is.null(txdb)) {
    genomeInfo <- genome(txdb)[1]
    tryCatch({
      return(getBSgenome(genomeInfo))
    }, error = function(e) {
      stop("Failed to load BSgenome for organism derived from txdb: ", e$message, call. = FALSE)
    })
  }

  # If no genome information is provided, and either txdb is not provided or the genome cannot be inferred from txdb
  stop("Argument 'genome' must be a character name of a BSgenome data package, a BSgenome object, or NULL with a valid 'txdb' provided.", call. = FALSE)
}
