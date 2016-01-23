##' Trim/truncate a set of reads
##'
##' @param lreads A list of ShortReadQ objects
##' @param trim_len The length reads will be truncated to; default is NULL (no length truncation)
##' @param trim5 The number of nucleotides to trim from the 5'-end; default is 0
##' @return A list of truncated ShortReadQ objects
##' @export
##' @keywords internal
trimReads <- function (lreads, trim_len=NULL, trim5=0){
  loginfo("trimReads.R/trimReads: trim reads...")

  ## trim reads
  lreads <- lapply(lreads, function(reads) truncateReads(reads, trim_len, trim5))

  loginfo("trimReads.R/trimReads: done")
  return(lreads)
}

##' Trim/truncate a set of reads
##'
##' @param reads A set of reads as ShortReadQ object
##' @param trim_len The length reads will be truncated to; default is NULL (no length truncation)
##' @param trim5 The number of nucleotides to trim from the 5'-end; default is 0
##' @return A truncated ShortReadQ object
##' @export
##' @keywords internal 
##' @importMethodsFrom ShortRead sread ShortReadQ id
truncateReads <- function(reads, trim_len=NULL, trim5=0) {
  if (length(reads)>0) {
    ## nucleotides
    start <- rep(1+trim5, length(reads))
    len <- width(reads)
    ## first, use trim5: make reads of zero length if trim5 is longer than read length
    z <- start>len
    start[z] <- len[z]+1
    ## second, use trim_len to make reads of the desired trim_len read lenth
    if (!is.null(trim_len)) {
        z <- len>(trim_len+start-1)
        len[z] <- trim_len+start[z]-1
    }
    treads <- subseq(sread(reads), start, len)
    
    ## quals
    quals <- quality(reads)
    tquals <- subseq(quality(quals), start, len)
    constructor <- get(class(quals))
    tquals <- constructor(tquals)
    
    ## rebuild ShortReadQ object
    ShortReadQ(sread=treads, quality=tquals, id=id(reads))
  } else return(reads)
}
