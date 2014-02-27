##' Read paired end BAM files with requested columns from the BAM
##'
##' This function returns a GRangesList containing information about
##' pair-end alignment data. It is essentially
##' Rsamtools::readGAlignments, but it uses Rsamtools::scanBam
##' to get information that cannot be obtained from
##' readBamGappedAlignments.
##' 
##' @title Read Paired End Bam Files
##' @param bam_file Path to a bam file
##' @param scan_bam_what additional paramters passed to scanBamParam()
##' @param remove.strandness A logical indicating whether read strands should be set to "*".
##' @return GRangesList
##' @author Cory Barr
##' @export
##' @keywords internal
##' @importMethodsFrom GenomicAlignments cigar 
##' @importFrom GenomicAlignments readGAlignments
readRNASeqEnds <- function(bam_file, scan_bam_what=NULL, remove.strandness=TRUE) {
  scan_bam_what0 <- union(scan_bam_what, "flag")
  sbp <- ScanBamParam(what=scan_bam_what0)
  ## readGappedAlignments() filters out unmapped reads (can happen in
  ## halfmapping_uniq files).
  gapped <- readGAlignments(bam_file, use.names=TRUE, param=sbp)
  ## Adding the number to the read name indicating which pair of the
  ## paired-end read it is.
  names(gapped) <- paste(names(gapped),
                         getEndNumber(values(gapped)$flag),
                         sep="#")
  grl <- grglist(gapped)

  ##adding strand, fusion, and junction info by unlisting then
  ##relisting. Preferred method for replacing parts of GRangesLists

  ## NOTE (H.P. - Nov 8, 2011): The fusion/junction info really describes
  ## the top-level elements of the GRangesList object 'grl', not its
  ## individual ranges. Therefore this info should be stored one level up
  ## i.e. in 'values(grl)', not in 'values(unlisted)'.
  ## In addition to being conceptually the right place for this info, that
  ## would allow the code to be much faster (no need to unlist/relist) and
  ## would also result in a smallest object (no need to repeat the same info
  ## for all the ranges belonging to the same top-level element).
  unlisted <- unlist(grl, use.names=FALSE)

  ## why? (H.P.) In any case, would be better to do this in 'gapped', before doing grl <- grglist(gapped)
  if (remove.strandness) strand(unlisted) <- "*"

  ## NOTE (H.P. - Nov 8, 2011): The junction col is redundant. You already
  ## know you have a junction when 'elementLengths(grl) >= 2' is TRUE. 
  junction_by_element <- grepl("N", cigar(gapped))
  values(unlisted)$junction <- rep.int(junction_by_element, elementLengths(grl))
  fusion_by_element <- grepl("H", cigar(gapped))
  values(unlisted)$fusion <- rep.int(fusion_by_element, elementLengths(grl))

  ## Handle user-specified columns.
  if (length(scan_bam_what) != 0L) {
    ## NOTE (H.P. - Nov 8, 2011): Like for the fusion/junction info,
    ## 'values(unlisted)' is the wrong place to store the user-specified
    ## columns.
    idx <- rep.int(seq_len(length(gapped)), elementLengths(grl))
    user_values <- values(gapped)[idx, scan_bam_what]
    values(unlisted) <- cbind(values(unlisted), user_values)
  }
  relist(unlisted, grl)
}
