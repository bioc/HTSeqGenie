##' Filter reads by quality
##'
##' Filtering reads by quality score.
##' Discards reads that have more than
##' a fraction of X nucleotides with a score below Y.
##'
##' X and Y are controlled by global config variables
##' X: filterQuality.minFrac
##' Y: filterQuality.minQuality
##' 
##' @param lreads A list of ShortReadQ objects
##' @return A list of quality filterered ShortReadQ objects
##' @keywords internal
##' @export
filterQuality <- function(lreads) {
  loginfo("filterQuality.R/filterQuality: filtering reads...")
  
  ## get config parameters
  paired_ends <- getConfig.logical("paired_ends")
  minquality <- getConfig.numeric("filterQuality.minQuality")
  minfrac <- getConfig.numeric("filterQuality.minFrac")
  minlength <- getConfig.numeric("filterQuality.minLength")
  if (is.null(minlength)) minlength <- 1

  ## trim reads based on low quality tail for illumina 1.5 and 1.8
  quality_encoding <- getConfig('quality_encoding')
  if (quality_encoding == 'illumina1.5') {
    lreads <- trimTailsByQuality(lreads, minqual='B')
  } else if (quality_encoding == 'illumina1.8') {
    lreads <- trimTailsByQuality(lreads, minqual='#')
  }
  
  ## remove reads that got trimmed too short or already came in short
  z <- filterByLength(lreads, minlength, paired_ends)
  lreads <- lapply(lreads, function(reads) reads[z])
  
  comb <- isAboveQualityThresh(lreads[[1]], minquality=minquality,
                                           minfrac=minfrac)
  if (paired_ends) {
    comb <- comb & isAboveQualityThresh(lreads[[2]], minquality=minquality,
                                                    minfrac=minfrac)
  }
  lreads <- lapply(lreads, function(reads) reads[comb])
  loginfo(paste("filterQuality.R/filterQuality: fraction of bad reads=", mean(!comb)))
  loginfo("filterQuality.R/filterQuality: done")
  
  return(lreads)
}

##' Check for high quality reads
##'
##' Checks whether reads have more than a fraction
##' of minFrac nucleotides with a score below minquality.
##'
##' @param reads A set of reads as ShortReadQ object
##' @param minquality Minimal quality score
##' @param minfrac Fraction of positions that need to be over
##'                minquality to be considered a good read.
##' @return A boolean vector indicating whether read is considered high quality.
##' @keywords internal
##' @export
##' @importMethodsFrom Biostrings unlist
##' @importMethodsFrom IRanges  unlist as.vector as.factor
isAboveQualityThresh <- function(reads, minquality, minfrac) {
  ## Going through the process of making a single concatted qual
  ## score so can use Bioc qual-score converter
  qualities <- quality(quality(reads))
  bss <- BStringSet(unlist(qualities))
  constructor <- get(class(quality(reads)))
  single_fqq <- constructor(bss)
  raw_quals <- as(single_fqq, "numeric")

  ## filter on quality
  above_thresh <- raw_quals >= minquality
  ends <- cumsum(width(qualities))
  starts <- ends - width(qualities) + 1
  thresh_views <- Views(Rle(above_thresh), start=starts, end=ends)
  means <- viewMeans(thresh_views)
  z <- means > minfrac

  ## return filter
  z
}

##' Filter reads by length
##'
##' Checks whether reads have at least a length of minlength.
##' Useful values are zero the rid of empty reads or 12 to match
##' the gsnap k-mer size.
##'
##' @param lreads A set of reads as ShortReadQ object
##' @param minlength Minimum length
##' @param paired Indicates whether lreads has one of two elements
##' @return A boolean vector indicating whether read passes filter
##' @keywords internal
##' @export
filterByLength <- function(lreads, minlength=12, paired=FALSE) {
  loginfo("filterQuality.R/filterQuality: filterByLength...")
  
  z <- width(lreads[[1]]) >= minlength
  if (paired) {
    z <- z & width(lreads[[2]]) >= minlength
  }
  loginfo(paste("filterQuality.R/filterByLength: fraction of filtered short reads=", mean(!z)))
  loginfo("filterQuality.R/filterByLength: done")
  
  return(z)
}

##' Trim off low quality tail
##' 
##' The illuminsa manuals states:
##' If a read ends with a segment of mostly low quality (Q15 or below),
##' then all of the quality values in the segment are replaced with a value of 2(
##' encoded as the letter B in Illumina's text-based encoding of quality scores)...
##' This Q2 indicator does not predict a specific error rate, but rather indicates
##' that a specific final portion of the read should not be used in further analyses.
##'
##' For illumina 1.8 the special char is encoded as '#', which we chhose as default here.
##' For illumina 1.5 make sure to set the minqual to 'B'
##' 
##' @param lreads A list (usually a pair) of ShortReadQ object
##' @param minqual An ascii encoded quality score
##' @return A list of quality trimmed ShortReadQ objects
##' @keywords internal
##' @importMethodsFrom ShortRead trimEnds
##' @export
trimTailsByQuality <- function (lreads, minqual="#"){  
  loginfo(paste("preprocessReads.R/preprocessReadsChunk: Starting trimTailsByQuality ..."))
  total_reads <- length(lreads[[1]])
  if(total_reads > 0){
    lreads <- lapply(lreads, function(x)  trimEnds(x, a=minqual, left=FALSE))
    loginfo(paste("preprocessReads.R/preprocessReadsChunk: done"))
  }
  return(lreads)
}
