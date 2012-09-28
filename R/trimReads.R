##' Trim/truncate a set of reads
##'
##' @param lreads A list of ShortReadQ objects
##' @param trim_len The length reads will be truncated to 
##' @return A list of truncated ShortReadQ objects
##' @export
##' @keywords internal
trimReads <- function (lreads, trim_len){
  loginfo("trimReads.R/trimReads: trim reads...")

  ## trim reads
  lreads <- lapply(lreads, function(reads) truncateReads(reads, trim_len))

  loginfo("trimReads.R/trimReads: done")
  return(lreads)
}

##' Trim/truncate a set of reads
##'
##' @param reads A set of reads as ShortReadQ object
##' @param trim_len The length reads will be truncated to 
##' @return A truncated ShortReadQ object
##' @export
##' @keywords internal
truncateReads <- function(reads, trim_len) {
  if (length(reads)>0) {
    ## nucleotides
    trunc <- width(reads)
    trunc[trunc>trim_len] <- trim_len   
    treads <- subseq(sread(reads), 1, trunc)
    
    ## quals
    quals <- quality(reads)
    tquals <- subseq(quality(quals), 1, trunc)
    constructor <- get(class(quals))
    tquals <- constructor(tquals)
    
    ## rebuild ShortReadQ object
    ShortReadQ(sread=treads, quality=tquals, id=id(reads))
  } else return(reads)
}
