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

##' Merge paired ends reads
##'
##' Merges overlapping paired end reads using the pear programm
##' via the peared package
##' 
##' @param reads A set of reads as ShortReadQ object
##' @param save_dir Save directory
##' @return A single ShortReadQ object
##' @export 
##' @importFrom peared pear
##' @keywords internal 
mergeReads <- function(reads, save_dir=NULL){

  pear.params <- getConfig("mergeReads.pear_param")
  if (missing(save_dir)) save_dir <- file.path(getConfig("save_dir"))

  tmp.dir <- makeDir(file.path(save_dir, "mergeReads"),
                     overwrite="erase")

  paths <- writeFastQFiles(reads, tmp.dir, "1.fq", "2.fq");

  args <- list(fq1=paths[[1]],
               fq2=paths[[2]],
               sample.id=getConfig("prepend_str"),
               out.dir=tmp.dir,
               temp.dir=tmp.dir)

  if(!is.null(pear.params)){
    args <- c(args, pear.param=pear.params)
  }
  
  result.paths <- do.call(pear, args)
  reads <- readFastq(result.paths[['assembled']])
  unlink(tmp.dir, recursive=TRUE)
  return(reads)
}
