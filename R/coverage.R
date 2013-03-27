##' Calculate read coverage
##'
##' @title  Calculate read coverage
##' @return Nothing
##' @author Jens Reeder
##' @keywords internal
##' @export
calculateCoverage <- function() {
  safeExecute({
    loginfo("coverage.R/calculateCoverage: starting...")

    ## get parameters
    paired_ends <- getConfig.logical("paired_ends") 
    extendReads <- getConfig.logical("coverage.extendReads")
    fragmentLength <- getConfig.numeric("coverage.fragmentLength")
    maxFragmentLength <- getConfig.numeric("coverage.maxFragmentLength")
    prepend_str <- getConfig("prepend_str")
    chunkdirs <- getChunkDirs()
    
    ## estimate fragmentLength from first chunk, if needed
    if (!paired_ends && extendReads && is.null(fragmentLength)) {
      chunkdir <- chunkdirs[1]
      bamfile <- file.path(chunkdir, "bams", paste(prepend_str, "analyzed.bam", sep="."))
      cov <- computeCoverage(bamfile, extendReads, paired_ends, maxFragmentLength=maxFragmentLength)
      fragmentLength <- attr(cov, "fragmentLength")
      loginfo(paste("coverage.R/calculateCoverage: estimated fragmentLength=", fragmentLength))
    }
        
    ## compute coverage for each chunk
    listIterator.init(chunkdirs)
    processChunks(listIterator.next, function(chunkdir, save_dir) {
      bamfile <- file.path(chunkdir, "bams", paste(prepend_str, "analyzed.bam", sep="."))
      cov <- computeCoverage(bamfile, extendReads, paired_ends, fragmentLength, maxFragmentLength)
      saveCoverage(cov, chunkdir, prepend_str)
    })
    
    ## merge coverage
    save_dir <- getConfig("save_dir")
    mergeCoverage(chunkdirs, save_dir, prepend_str)
  }, memtracer=getConfig.logical("debug.tracemem"))
}

##' Merge coverage files
##'
##' @title Merge coverage files
##' @param indirs A character vector, indicating which directories have to be merged
##' @param outdir A character string indicating the output directory (which must exist)
##' @param prepend_str A character string, containing a prefix going to be appended on all output result files
##' @return Nothing
##' @author Gregoire Pau
##' @keywords internal
##' @export
mergeCoverage <- function(indirs, outdir, prepend_str) {
  ## merge coverages
  cov <- NULL
  for (indir in indirs) {
    filename <- dir(file.path(indir, "results"), pattern="coverage.RData", full.names=TRUE)
    cov0 <- get(load(filename))
    if (is.null(cov)) cov <- cov0
    else cov <- cov + cov0
  }

  ## restore seqlengths
  seqlengths(cov) <- seqlengths(cov0)

  ## save coverage
  saveCoverage(cov, outdir, prepend_str)
}

##' Compute the coverage vector given a bamfile
##'
##' @title Compute the coverage vector given a bamfile
##' @param bamfile A character string indicating the path of bam file
##' @param extendReads A logical, indicating whether reads should be extended
##' @param paired_ends A logical, indicating whether reads are paired
##' @param fragmentLength An integer, indicating the new size of reads when extendReads is TRUE
##' and paired_ends is FALSE. If NULL, read size is estimated using estimate.mean.fraglen from the chipseq package.
##' @param maxFragmentLength An optional integer, specifying the maximal size of fragments.
##' Longer fragments will be disregarded when computing coverage.
##' @return A SimpleRleList object containing the coverage
##' @author Gregoire Pau
##' @keywords internal
##' @export
computeCoverage <- function(bamfile, extendReads=FALSE, paired_ends=FALSE, fragmentLength=NULL,
                            maxFragmentLength=NULL) {
  ## is paired end?
  if (paired_ends) {
    reads <- readRNASeqEnds(bamfile, remove.strandness=TRUE)
    reads <- consolidateByRead(reads)
    if (extendReads) reads <- unlist(range(reads))
  }
  else {
    reads <- readRNASeqEnds(bamfile, remove.strandness=FALSE)
    reads <- unlist(reads)
    if (extendReads) {
      if (is.null(fragmentLength)) {
        fragmentLength <- median(chipseq::estimate.mean.fraglen(reads, method="coverage"), na.rm=TRUE)
      }
      reads <- resize(reads, width=fragmentLength)
    }
  }

  ## keep reads shorter than maxFragmentLength
  if (!is.null(maxFragmentLength)) reads <- reads[width(reads)<maxFragmentLength]

  ## compute coverage
  cov <- coverage(reads)
  seqlengths(cov) <- seqlengths(reads)
 
  ## return coverage
  attr(cov, "fragmentLength") <- fragmentLength
  cov
}

saveCoverage <- function(cov, save_dir, prepend_str) {
  ## save coverage as a RData file
  saveWithID(cov, "coverage", prepend_str, save_dir=file.path(save_dir, "results"), compress=FALSE)

  ## save coverage as a BigWig file
  bwfile <- file.path(save_dir, "results", paste(prepend_str, "coverage", "bw", sep="."))
  rd <- as(cov, "RangedData")
  seqlengths(rd) <- seqlengths(cov)
  rtracklayer::export(rd, con=bwfile, format="bw")
  
  invisible(TRUE)
}
