##' Does a SAM flag indicate the first fragment
##'
##' Compute whether a SAM/BAM flag indicates a
##' first fragment. Method is not foolproof, as it ignores
##' a lot of SAM semantics. E.g the SAM spec says:
##' "If 0x1 is unset, no assumptions can be made about
##' 0x2, 0x8, 0x20, 0x40 and 0x80".
##' For our purpose this should be enough, but we should
##' keep an open eye for a more robust implementation in Rsamtools.
##'
##' @param flag A flag from the BAM/SAM file
##' @keywords internal
##' @return Logical
isFirstFragment <- function(flag) {
  int <- flag %% 256
  int < 8 * 16 ##per samtools flag specs
}

##' Returns the end number of an end from a paired-end read
##'
##' Given an integer from the BAM flag field, tells which end it is in
##' a read
##' @title Get Read End Number 
##' @param int an int from a SAM flag
##' @return {1,2}
##' @author Cory Barr
##' @keywords internal
##' @export
getEndNumber <- function(int) {
  isFirstFrag <- isFirstFragment(int)
  ifelse(isFirstFrag, 1L, 2L)
}

mergeBams <- function(inbams, outbam, sort=TRUE, index=TRUE) {
  if (length(inbams)>1) {
    ## merge
    mergeBam(inbams, outbam, overwrite=TRUE)
    
    ## sort outbam?
    if (sort) {
      sortedbam <-  paste(outbam, "_sorted", sep="")
      sortedbam <- sortBam(outbam, sortedbam)
      file.copy(sortedbam, outbam, overwrite=TRUE)
      unlink(sortedbam)
    }
  } else {
    file.copy(inbams[1], outbam)
  }
  
  ## index outbam?
  if (index) indexBam(outbam)

  ## return outbam
  invisible(outbam)
}

## Rsamtools does not expose samtools cat, so we have to use the cmd line
## check new htslib/Rsamtool when it gets released
catBams <- function(inbams, outbam){

  if(length(inbams)==1){
    ## renaming across file systems fails, so we have to copy instead.
    ## since bams are large we try renaming first, then copying
    file.rename(inbams[[1]], outbam) || file.copy(inbams[[1]], outbam)
  } else {
    if(system2("samtools", stderr=FALSE) == 0){
      ## this is mostly to make sure it runs on the test servers, which might not have samtools installed
      warning("Problems calling samtools. Falling back to Rsamtools merge, which might be slower")
      ## mergeBams already indexes, so we can return early 
      return(mergeBams(inbams, outbam))
    }
    args <- paste( c("cat", "-o", paste(outbam), inbams))
    
    retcode <- system2( "samtools", args)
    if(retcode != 0){
      stop( paste("samtools command [ ", paste(args, collapse=" "), "] failed."))
    }
  }
  indexBam(outbam)
  return(outbam)
}


##' Uniquely count number of reads in bam file
##'
##'
##' @title Uniquely count number of reads in bam file
##' @param bam Name of bam file
##' @return number of reads
##' @author Jens Reeder
##' @keywords internal
##' @importMethodsFrom Rsamtools ScanBamParam scanBam
##' @export
bamCountUniqueReads <- function(bam) {
  sb_param <- ScanBamParam(what=c("qname"))
  res <- scanBam(bam, param=sb_param)[[1]]
  nbreads <- length(unique(res$qname))
  nbreads
}

##' Compute record statistics from a bam file
##'
##' The statistics are additive over chunks/lanes.
##' @title Compute record statistics from a bam file
##' @param bam A character string containing an existing bam file
##' @return A numeric vector
##' @author Gregoire Pau
##' @keywords internal
##' @export
computeBamStats <- function(bam) {
  ## read bam
  res <- scanBam(bam, param=ScanBamParam(what=c("cigar", "strand", "mapq", "rname", "isize"), tag=c("NM")))[[1]]
  zna <- !is.na(res$rname)

  ## preprocess cigar strings
  events <- c("M", "N", "S", "I", "D")
  cigarletters <- paste(gsub("[0-9]","", res$cigar[zna]), collapse="")
  cigarletters <- table(factor(strsplit(cigarletters, "")[[1]], events))
  storage.mode(cigarletters) <- "numeric"
  
  ## preprocess isize
  target_length <- abs(res$isize[zna])
  cuts <- c(seq(1, 600, by=10),Inf)
  target_length <- table(cut(target_length, cuts, right=TRUE))
  names(target_length) <- sprintf("target_length.bin.%d", head(cuts, -1))

  ## compute nm, mapq and nchar.cigar
  nm <- as.numeric(res$tag$NM[zna])
  mapq <- as.numeric(res$mapq[zna])
  nchar.cigar <- nchar(res$cigar[zna])

  ## chromosomes
  chromosomes <- table(res$rname)
  storage.mode(chromosomes) <- "numeric"
  names(chromosomes) <- paste("chromosome.bin.", names(chromosomes), sep="")
  
  stats <- c(
             ## general
             total.nbrecords=length(res$rname),
             mapped.nbrecords=sum(zna),
             
             ## NM
             nm.sum=sum(nm),
             nm.sum2=sum(nm^2),
             nm.bin.0=sum(nm==0),
             nm.bin.1=sum(nm==1),
             nm.bin.2=sum(nm==2),
             nm.bin.3=sum(nm==3),
             nm.bin.4=sum(nm==4),
             nm.bin.5=sum(nm==5),
             nm.bin.6=sum(nm==6),
             nm.bin.7=sum(nm==7),
             nm.bin.8=sum(nm>=8),
             
             ## mapq
             mapq.sum=sum(mapq),
             mapq.sum2=sum(mapq^2),
             mapq.bin.36=sum(mapq<=36),
             mapq.bin.37=sum(mapq==37),
             mapq.bin.38=sum(mapq==38),
             mapq.bin.39=sum(mapq==39),
             mapq.bin.40=sum(mapq==40),
             mapq.bin.41=sum(mapq>=41),
             
             ## cigar
             cigar.length.sum=sum(nchar.cigar),
             cigar.length.sum2=sum(nchar.cigar^2),
             cigar.S.sum=as.numeric(cigarletters["S"]),
             cigar.N.sum=as.numeric(cigarletters["N"]), 
             cigar.ID.sum=as.numeric(cigarletters["I"]+cigarletters["D"]),
             cigar.withS.nbrecords=length(grep("S", res$cigar[zna])),
             cigar.withN.nbrecords=length(grep("N", res$cigar[zna])),
             cigar.withID.nbrecords=length(grep("[ID]", res$cigar[zna])),
             
             ## chromosomes
             chromosomes,

             ## target_length bins
             target_length
             )
  
  stats
}
