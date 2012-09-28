##' Calculate read coverage
##'
##' @title  Calculate read coverage
##' @return Nothing
##' @author Jens Reeder
##' @export
##' @keywords internal
calculateCoverage <- function() {
  ## get parameters
  save.dir <- getConfig("save_dir")
  extendReads <- getConfig.logical("coverage.extendReads")
  paired_ends <- getConfig.logical("paired_ends") 
  fragmentLength <- getConfig.numeric("coverage.fragmentLength")
    
  safeExecute({
    bam.file <- file.path(save.dir, 'bams', paste(getConfig('prepend_str'), 'analyzed.bam', sep="."))
    storeCoverage(bam.file, extendReads, paired_ends, fragmentLength)
  }, memtracer=getConfig.logical("debug.tracemem"))
}

##' Save coverage information 
##'
##' Saves the coverage as RData and BigWig
##' 
##' @title Save coverage information
##' @param bam.file Name of bam file.
##' @param extendReads A logical indicating whether reads should be extended. If TRUE, reads are extended using mate-pair information if paired_ends is TRUE and using fragmentLength otherwise. Default is FALSE.
##' @param paired_ends A logical indicating whether reads are paired. Default is FALSE.
##' @param fragmentLength An optional numeric, indicating the fragment length. If NULL, fragment length is estimated from the reads using \code{chipseq::estimate.mean.fraglen}.
##' Default is NULL.
##' @return coverage
##' @author Cory Barr, Jens Reeder
##' @importFrom chipseq estimate.mean.fraglen
##' @export
##' @keywords internal
storeCoverage <-function(bam.file, extendReads=FALSE, paired_ends=FALSE, fragmentLength=NULL){
  ## load reads
  loginfo("coverage.R/storeCoverage: Calculating coverage from the merged BAM files...")
  reads <- readRNASeqEnds(bam.file, remove.strandness=FALSE)
  
  ## is paired end? 
  if (paired_ends) {
    ## reads <- readBamGappedAlignmentPairs(bam.file) ## readBamGappedAlignmentPairs() doesn't work on a mixture of concordant/halfmapping read pairs
    ## reads <- grglist(reads, drop.D.ranges=TRUE)
    reads <- consolidateByRead(reads)
    if (extendReads) reads <- unlist(range(reads))
  }
  else {
    reads <- unlist(reads)
    if (extendReads) {
      if (is.null(fragmentLength)) {
        fragmentLength <- median(chipseq::estimate.mean.fraglen(reads, method = "coverage"), na.rm=TRUE)
        loginfo(paste("coverage.R/storeCoverage: estimated fragmentLength=", fragmentLength))
      }
      reads <- resize(reads, width=fragmentLength)
    }
  }

  ## compute coverage
  cov <- coverage(reads)
  cov <- as(cov, "RangedData")
  seqlengths(cov) <- seqlengths(reads)
  
  ## make BigWig file
  bw_file <- file.path(getConfig('save_dir'), "results",
                       paste(getConfig('prepend_str'), "coverage", "bw", sep="."))
  rtracklayer::export(cov, con = bw_file, format= "bw")

  filename <- saveWithID(cov, "coverage",
                          getConfig('prepend_str'),
                          save_dir=file.path(getConfig('save_dir'), "results"),
                          compress=FALSE)
  loginfo("coverage.R/storeCoverage: ...done")
  return(filename)
}
