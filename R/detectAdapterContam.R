##' Detect sequencing adapter contamination
##'
##' For each read or pair of read, search for specific
##' Illumina adapter sequences in the read. Flag if at least
##' one read has significant overlap with adapter.
##'
##' @param lreads List of reads as ShortRead objects
##' @param save_dir Save directory of a pipeline run
##' @return Boolean vector indicating vector contamination for each read
##' @export
##' @keywords internal
detectAdapterContam <- function(lreads, save_dir=NULL) {
  ## init
  if (missing(save_dir)) save_dir <- file.path(getConfig("save_dir"))
  setUpDirs(save_dir, overwrite="overwrite")
  loginfo("detectAdapterContam.R/detectAdapterContam: detecting adapter contamination...")

  ## get config parameters
  force_paired_end_adapter <- getConfig.logical("detectAdapterContam.force_paired_end_adapter")
  paired_ends <- getConfig.logical("paired_ends")
  prepend_str <- getConfig("prepend_str")

  ## cutoff parameter, estimated by bootstrapping with estimateCutoffs()
  adapter_cutoff_score <- 13.87229

  ## detect adapter contam in forward reads
  adapter_seqs <- getAdapterSeqs(paired_ends, force_paired_end_adapter, pair_num=1)
  adapter_contaminated_1 <- isAdapter(lreads[[1]], score_cutoff=adapter_cutoff_score,
                                      adapter_seqs=adapter_seqs)
  ids <- as.character(id(lreads[[1]])[adapter_contaminated_1])
  saveWithID(ids, 'adapter_contaminated_1', id=prepend_str,
             save_dir=file.path(save_dir, "results"), compress=FALSE)
  adapter_contaminated <- adapter_contaminated_1

  ## detect adapter contam in reverse reads
  if (paired_ends) {
    adapter_seqs <- getAdapterSeqs(paired_ends, force_paired_end_adapter, pair_num=2)
    adapter_contaminated_2 <- isAdapter(lreads[[2]], score_cutoff=adapter_cutoff_score,
                                        adapter_seqs=adapter_seqs)
    ids <- as.character(id(lreads[[2]])[adapter_contaminated_2])
    saveWithID(ids, 'adapter_contaminated_2', id=prepend_str,
               save_dir=file.path(save_dir, "results"), compress=FALSE)
    adapter_contaminated <- adapter_contaminated_1 & adapter_contaminated_2
  }
  
  loginfo(paste("detectAdapterContam.R/detectAdapterContam: fraction of contaminated reads=", mean(adapter_contaminated)))
  return(adapter_contaminated)
}

mergeDetectAdapterContam <- function (indirs, outdir, prepend_str, paired_ends) {
  ## set the files to merge
  tfilenames <- c("adapter_contaminated_1")
  if (paired_ends) tfilenames <- c(tfilenames, "adapter_contaminated_2")

  for (tfilename in tfilenames) {
    ## guess filename from tfilename
    infiles <- sapply(indirs, function(indir) getObjectFilename(file.path(indir, "results"), tfilename))

    ## merge files
    loginfo(paste("detectAdapterContam.R/mergeDetectAdapterContam: merging file=", tfilename, "..."))
    allda <- unlist(lapply(infiles, function(z) get(load(z))))
    saveWithID(allda, tfilename, id=prepend_str, save_dir=file.path(outdir, "results"), compress=TRUE)
  }
}

##' Read list of Illumina adapter seqs from package data
##' 
##' @param paired_ends Dow we have paired ends reads?
##' @param force_paired_end_adapter Force paired end adapters for single end reads?
##' @param pair_num 1 for forward read, 2 for reverse read
##' @return The adapter seq as string
##' @keywords internal
getAdapterSeqs <- function(paired_ends, force_paired_end_adapter, pair_num=1) {
  if (paired_ends) {
    if (pair_num==1) {
      adapter_file <- "extdata/adapter_data/Illumina_adapters.curated.paired_end_1.fasta"
    } else if (pair_num==2) {
      adapter_file <- "extdata/adapter_data/Illumina_adapters.curated.paired_end_2.fasta"
    } else {
      stop( "Invalid pair number. Allowed values: 1,2" )
    }
  } else {
    if (force_paired_end_adapter) {
      adapter_file <- "extdata/adapter_data/Illumina_adapters.curated.paired_end_1.fasta"      
    } else {
      adapter_file <- "extdata/adapter_data/Illumina_adapters.curated.single_end.fasta"
    }
  }
  
  adapter_seqs <- readDNAStringSet(getPackageFile(adapter_file))
  adapter_seqs <- as.character(adapter_seqs)[[1]]

  return(adapter_seqs)
}

##' Detect adapter contamination
##'
##' Does a Needleman-Wunsch like small-in-large alignment
##' of the adapter vs each read. Flag read if score exceeds threshold
##' 
##' @param reads Set of reads as ShortRead object
##' @param score_cutoff Alignment score threshold that needs to be exceeded to be flagged as adapter.
##'                     Usually this value is determined empirically by getAdpaterThreshold()
##' @param adapter_seqs One or more adapter sequences
##' @return boolean vector indicating adapter contamination
##' @export
##' @keywords internal
isAdapter <- function(reads,
                      score_cutoff,
                      adapter_seqs){

  if(class(reads) == "ShortReadQ") reads <- sread(reads)
  
  res <- pairwiseAlignment(reads,
                           adapter_seqs,
                           type       = "overlap",
                           gapOpening = -Inf,
                           scoreOnly  = TRUE)

  res <- res > score_cutoff
}


##' Estimate an adopter alignment cutoff scor
##'
##' Empirically estimate a threshold that discriminates
##' random reads from reads with adapter contamination
##' @param read_len The read length
##' @param n Number of samples
##' @keywords internal
getRandomAlignCutoff <- function(read_len, n) {
  adapter_seqs <- readFasta(getPackageFile("extdata/adapter_data/Illumina_adapters.paired_end.fasta"))
  adapter_seqs<- as.character(sread(adapter_seqs))[1]

  ## generate random reads for null alignment hypothesis
  randomStrings <- matrix(sample(DNA_ALPHABET[1:4], read_len*n, replace=TRUE), nrow=n)
  randomStrings <- apply(randomStrings, 1, paste, collapse = "")
  randomStrings <- DNAStringSet(randomStrings)

  randomAlignScores <- pairwiseAlignment(randomStrings,
                                         adapter_seqs,
                                         type       = "overlap",
                                         gapOpening = -Inf,
                                         scoreOnly  = TRUE)
  return(randomAlignScores)
}

estimateCutoffs <- function() {
  readlengths <- c(50, 75, 100, 125)
  n <- 200000 ## sample size
  pval.cutoff <- 1e-4 ## false positive rate
  p <- 12 ## replicates

  cutoffs <- do.call(cbind, lapply(readlengths, function(rl) {
    z <- mclapply(1:p, function(i) getRandomAlignCutoff(rl, n))
    sapply(z, quantile, 1-pval.cutoff)
  }))
  colnames(cutoffs) <- readlengths
  
  print(cutoffs)
  ## cutoffs seem to be independent of read length
  ## cutoff = 13.87229
}
