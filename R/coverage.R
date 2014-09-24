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
    loginfo("coverage.R/calculateCoverage: done")
  }, memtracer=getConfig.logical("debug.tracemem"))
}

##' Check coverage for sparseness
##'
##' Some Rle related operations become very slow when they are dealing with
##' data that violates their sparseness assumption. This method provides an estimate about
##' whether the data is dense or sparse. More precicely it checks
##' if the fraction of the number of runs over the total length is smaller than a threshold
##' 
##' @title isSparse
##' @param cov A cov object as SimpleRleList
##' @param threshold Fraction of number of runs over total length
##' @return Boolean whether this object is dense or sparse
##' @importMethodsFrom IRanges runLength unlist length
##' @author Jens Reeder
isSparse <- function(cov, threshold=0.1) {
  density <- sum(as.numeric(unlist(lapply(runLength(cov), length)))) /
             sum(as.numeric(lapply(cov, length)))
  return(density < threshold)
}

##' Merge coverage files
##'
##' Merges coverage objects, usually SipleRleLists, in a tree-reduce
##' fashion. The coverage object dynamically switches to a SimpleIntegerList,
##' once the data becomes too dense. 
##' 
##' @title Merge coverage files
##' @param indirs A character vector, indicating which directories have to be merged
##' @param outdir A character string indicating the output directory (which must exist)
##' @param prepend_str A character string, containing a prefix going to be appended on all output result files
##' @return Nothing
##' @author Jens Reeder
##' @keywords internal
##' @export
mergeCoverage <- function(indirs, outdir, prepend_str) {

  .merge <- function(covs){
    l <- length(covs)
    if (l == 1){
      return( covs[[1]]) 
    } else {
      cov.a = .merge(covs[seq(1, l, 2)]) ## take every second, starting at 1
      cov.b = .merge(covs[seq(2, l, 2)]) ## and starting at 2
      cov.combined <- cov.a + cov.b

      if( class(cov.combined) != "SimpleIntegerList"
         && ! isSparse(cov.combined)){
        cov.combined <- as(cov.combined, "SimpleIntegerList")
      }
      return (cov.combined)
    }
  }
  
  covs <- lapply(indirs,
         function(indir){
           filename <- dir(file.path(indir, "results"), pattern="coverage.RData", full.names=TRUE)
           cov.next <- get(load(filename))
           ## SimpleRleLists can be bigger, useful for WGS
           return( as(cov.next, 'SimpleRleList'))
         }
         )
    
  coverage <- .merge(covs)
  
  ## restore seqlengths and fragmentLength
  seqlengths(coverage) <- seqlengths(covs[[1]])
  attr(coverage, "fragmentLength") <- attr(covs[[1]], "fragmentLength")
  saveCoverage(coverage, outdir, prepend_str)
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
##' @importMethodsFrom BiocGenerics table
##' @importFrom rtracklayer export
##' @importMethodsFrom chipseq estimate.mean.fraglen
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
        ## fragmentLength is estimated on each chromosome by estimate.mean.fraglen()
        ## global fragmentLength estimation is computed using a weighted mean on reads per chromosomes
        ## (to minimize inaccurate fragment estimation on chromosomes with few reads)
        fraglen.chr <- estimate.mean.fraglen(reads, method="coverage")
        ctfrag.chr <- table(seqnames(reads))
        fragmentLength <- round(sum(fraglen.chr * ctfrag.chr[names(fraglen.chr)])/sum(ctfrag.chr))
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
  rtracklayer::export(cov, con=bwfile, format="bw")

  ## fragmentLength
  fragmentLength <- attr(cov, "fragmentLength")
  if (is.null(fragmentLength)) fragmentLength <- NA
  
  ## compute coverage stats
  cm <- sapply(1:length(cov), function(i) mean(cov[[i]]))
  chr.length <- as.numeric(seqlengths(cov))
  coverage.mean.genome <- sum(chr.length*cm)/sum(chr.length)
  coverage.mean <- setNames(cm, paste0("coverage.mean.", names(cov)))
  cf1 <- sapply(1:length(cov), function(i) mean(cov[[i]]>=1))
  cf2 <- sapply(1:length(cov), function(i) mean(cov[[i]]>=2))
  cf5 <- sapply(1:length(cov), function(i) mean(cov[[i]]>=5))
  coverage.fraction.atleast1read.genome <- sum(chr.length*cf1)/sum(chr.length)
  coverage.fraction.atleast2read.genome <- sum(chr.length*cf2)/sum(chr.length)
  coverage.fraction.atleast5read.genome <- sum(chr.length*cf5)/sum(chr.length)
  coverage.fraction.atleast1read <- setNames(cf1, paste0("coverage.fraction.atleast1read.", names(cov)))
  coverage.fraction.atleast2read <- setNames(cf2, paste0("coverage.fraction.atleast2read.", names(cov)))
  coverage.fraction.atleast5read <- setNames(cf5, paste0("coverage.fraction.atleast5read.", names(cov)))
  dat <- c(estimated.fragmentLength=fragmentLength,
           coverage.mean.genome=coverage.mean.genome,
           coverage.fraction.atleast1read.genome=coverage.fraction.atleast1read.genome,
           coverage.fraction.atleast2read.genome=coverage.fraction.atleast2read.genome,
           coverage.fraction.atleast5read.genome=coverage.fraction.atleast5read.genome,
           coverage.mean,
           coverage.fraction.atleast1read,
           coverage.fraction.atleast2read,
           coverage.fraction.atleast5read
           )
  
  ## save summary results
  df <- data.frame(name=names(dat), value=dat)
  saveWithID(df, "summary_coverage", id=prepend_str, save_dir=file.path(save_dir, "results"), format="tab")
  
  invisible(TRUE)
}
