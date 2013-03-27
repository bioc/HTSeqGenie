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
    
  comb <- HTSeqGenie::isAboveQualityThresh(lreads[[1]], minquality=minquality,
                                           minfrac=minfrac, minlength=minlength)
  if (paired_ends) {
    comb <- comb & HTSeqGenie::isAboveQualityThresh(lreads[[2]], minquality=minquality,
                                                    minfrac=minfrac, minlength=minlength)
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
##' @param minlength The minimum read length. Default is no filter on read length.
##' @return A boolean vector indicating whether read is considered high quality.
##' @keywords internal
##' @export
##' @importMethodsFrom Biostrings unlist
##' @importMethodsFrom IRanges  unlist as.vector as.factor
isAboveQualityThresh <- function(reads, minquality, minfrac, minlength=NULL) {
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

  ## filter on read length
  if (!is.null(minlength)) z <- z & width(reads)>=minlength

  ## return filter
  z
}
