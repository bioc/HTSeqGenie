##' Align reads against genome
##'
##' @title  Align reads against genome
##' @return Nothing
##' @author Gregoire Pau
##' @keywords internal
##' @export
alignReads <- function() {
  safeExecute({
    loginfo(paste("alignReads.R/alignReads: starting alignment..."))
    
    ## get the chunks
    chunkdirs <- getChunkDirs()
    listIterator.init(chunkdirs)
    
    ## schedule alignment for each chunk
    num_cores <- getConfig.integer("num_cores")
    nbthreads_perchunk <- getConfig.integer("alignReads.nbthreads_perchunk")
    nb.parallel.jobs <- floor(num_cores/nbthreads_perchunk)
    processChunks(listIterator.next, function(chunkdir, save_dir) {
      fp1 <- file.path(chunkdir, "bams", "processed.aligner_input_1.fastq")
      fp2 <- file.path(chunkdir, "bams", "processed.aligner_input_2.fastq")
      if (!file.exists(fp2)) fp2 <- NULL ## for single end reads
      alignReadsChunk(fp1=fp1, fp2=fp2, save_dir=save_dir)
    }, nb.parallel.jobs=nb.parallel.jobs)
    
    ## merge bams
    save_dir <- getConfig("save_dir")
    prepend_str <- getConfig("prepend_str")
    num_cores <- getConfig.integer("num_cores")
    mergeAlignReads(chunkdirs, save_dir, prepend_str, num_cores)

    ## write report
    writePreprocessAlignReport()
    
    loginfo(paste("alignReads.R/alignReads: done"))
  }, memtracer=getConfig.logical("debug.tracemem"))
}

##' Genomic alignment using gsnap.
##'
##' Aligns reads in fp1 and fp2 to genome specified via
##' global config variable alignReads.genome. Gsnap output
##' is converted into BAM files and sorted + indexed.
##' 
##' @title Genomic alignment
##' @param fp1 Path to FastQ file
##' @param fp2 Path to second FastQ file if paired end data, NULL if single ended
##' @param save_dir Save directory
##' @keywords internal
##' @return Nothing
alignReadsChunk <- function(fp1, fp2=NULL, save_dir=NULL) {
  ## init
  if (missing(save_dir)) save_dir <- file.path(getConfig("save_dir"))
  setUpDirs(save_dir, overwrite="overwrite")
  loginfo(paste("alignReads.R/alignReadsChunk: running gsnap..."))
          
  ## get config parameters
  paired_ends <- getConfig.logical("paired_ends")
  prepend_str <- getConfig("prepend_str")
  genome <- getConfig("alignReads.genome")

  ## build parameters
  gsnapParams <- buildGsnapParams()
  
  ## call gsnap
  bam_dir <- file.path(save_dir, "bams")
  output <- file.path(bam_dir, prepend_str)
  gsnapOutput <- gsnap(fp1, fp2, gsnapParams, output=output)
  
  ## convert gsnap output
  buildAnalyzedBam(bam_dir, prepend_str)
  
  ## create summary
  createSummaryAlignment(save_dir, prepend_str)
  
  loginfo("alignReads.R/alignReadsChunk: done")
}

parseAdditionalParameters <- function(ap) {
  ## split ap
  sap <- strsplit(ap, " ")[[1]]

  ## remove trailing spaces
  z <- sap==""
  if (any(z)) sap <- sap[!z]
  
  ## checks the absence of short-form parameters (starting with single hyphens)
  z <- grepl("^-[^-]", sap)
  if (any(z)) stop("alignReads.R/parseAdditionalParameters: cannot process short-form parameter(s) from [",
                   paste(sap[z], collapse=", "), "] in config parameter 'alignReads.additional_parameters'")

  ## checks that all parameters are long-form (starting with double hyphens and with =)
  z <- grepl("^--", sap) & grepl("=", sap)
  if (any(!z)) stop("alignReads.R/parseAdditionalParameters: invalid long-form parameter(s) syntax [",
                    paste(sap[!z], collapse=", "), "] in config parameter 'alignReads.additional_parameters'")

  ## get parameter names
  z <- regexpr("=", sap)
  apars <- substr(sap, z+1, nchar(sap))
  names(apars) <- substr(sap, 3, z-1)
  
  ## checks no parameter is duplicated
  z <- duplicated(names(apars))
  if (any(z)) stop("alignReads.R/parseAdditionalParameters: duplicate parameter(s) [",
                   paste(sap[z], collapse=", "), "] in config parameter 'alignReads.additional_parameters'")
  
  as.list(apars)
}

##' Build gsnap paramters from current config
##'
##' @title Build Aligner Command-Line Args
##' @return a GsnapParam object
##' @author Gregoire Pau
##' @keywords internal
##' @export
buildGsnapParams <- function() {
  ## get config parameters
  path.gsnap_genomes <- getConfig("path.gsnap_genomes")
  genome <- getConfig("alignReads.genome")
  max_mismatches <- getConfig.integer("alignReads.max_mismatches")
  quality_encoding <- getConfig("quality_encoding")
  snp_index <- getConfig("alignReads.snp_index")
  splice_index <- getConfig("alignReads.splice_index")
  nbthreads_perchunk <- getConfig.integer("alignReads.nbthreads_perchunk")
  sam_id <- getConfig("alignReads.sam_id")
  additional_parameters <- getConfig("alignReads.additional_parameters")
  
  ## build gmapGenome
  gmapGenome <- GmapGenome(genome, directory=path.gsnap_genomes, create=FALSE)

  ## convert quality_encoding to gsnap quality protocol
  quality_protocol <- switch(quality_encoding,
                             sanger="sanger",
                             solexa="illumina",
                             illumina1.3="illumina",
                             illumina1.5="illumina",
                             illumina1.8="sanger")
  
  ## parse additional parameters
  apars <- parseAdditionalParameters(additional_parameters)

  ## isolate the parameters that have to be passed to GsnapParam separately: suboptimal-levels, mode, npaths, novelsplicing, batch
  dpar <- list("suboptimal-levels"="0", "mode"="standard", "npaths"="100", "novelsplicing"="0", "batch"="2")
  
  ## add default parameters, if missing
  z <- names(dpar)%in%names(apars)
  if (any(!z)) apars <- c(apars, dpar[!z])

  ## build extrapars
  z <- names(apars)%in%names(dpar)
  extrapars <- c(apars[!z], "quality-protocol"=quality_protocol)
  if (!is.null(sam_id)) extrapars <- c(extrapars, "read-group-id"=sam_id)

  ## replace - by _ in parameter names
  names(extrapars) <- gsub("-", "_", names(extrapars))
                           
  ## check that extrapars does not contain GsnapParam arguments
  gpargs <- c("genome", "unique_only", "max_mismatches", "suboptimal_levels",
              "mode", "snps", "npaths", "quiet_if_excessive", "nofails",
              "split_output", "novelsplicing", "splicing",
              "nthreads", "part", "batch")
  z <- names(extrapars)%in%gpargs
  if (any(z)) stop("alignReads.R/buildGsnapParams: additional parameters cannot contain: ", paste(names(extrapars)[z], collapse=", "))
  
  ## build GsnapParam
  do.call(GsnapParam, c(list(genome=gmapGenome,
                             unique_only=FALSE,
                             max_mismatches=max_mismatches,
                             suboptimal_levels=apars[["suboptimal-levels"]],
                             mode=apars[["mode"]],
                             snps=snp_index,
                             npaths=apars[["npaths"]],
                             quiet_if_excessive=FALSE,
                             nofails=FALSE,
                             split_output=TRUE,
                             novelsplicing=ifelse(apars[["novelsplicing"]]=="1", TRUE, FALSE),
                             splicing=splice_index,
                             nthreads=nbthreads_perchunk,
                             part=NULL,
                             batch=apars[["batch"]]
                             ), extrapars))
}

## create the analyzed bam
buildAnalyzedBam <- function(bam_dir, prepend_str) {
  ## get uniq bams
  inbams <- grep("_uniq.*\\.bam$", dir(bam_dir, full.names=TRUE), value=TRUE)
  
  ## merge inbams
  outbam <- paste(prepend_str, ".analyzed.bam", sep="")
  outbam <- file.path(bam_dir, outbam)

  ## merge, sort and index bams
  mergeBamSortIndex(inbams, outbam)
}

## given a bam file, return the bam type
getBamType <- function(bams) {
  bamtypes <- unname(sapply(basename(bams), function(bam) {
    pieces <- unlist(strsplit(bam, "\\."))
    paste(tail(pieces, n=2), collapse=".")
  }))
  gsub("\\.bam$", "", bamtypes)
}

## merges similarly named BAM residing in different directories
mergeBAMsAcrossDirs <- function(indirs, outdir, prepend_str, nb.parallel.jobs=1) {
  loginfo("alignReads.R/mergeBAMsAcrossDirs: starting...")

  ## prepare bamindirs and bamoutdir
  bamindirs <- file.path(indirs, "bams")
  bamoutdir <- file.path(outdir, "bams")

  ## split bam files per bam type
  inbams <- dir(bamindirs, full.names=TRUE, pattern="\\.bam$")
  bamtypes <- getBamType(inbams)
  sinbams <- split(inbams, bamtypes)

  ## merge files
  outbams <- paste(prepend_str, names(sinbams), "bam", sep=".")
  outbams <- file.path(bamoutdir, outbams)
  z <- mclapply(seq_len(length(sinbams)), function(i) {
    inbams <- sinbams[[i]]
    outbam <- outbams[i]

    ## merge, sort and index bam
    mergeBamSortIndex(inbams, outbam)
  }, mc.cores=nb.parallel.jobs)
  
  loginfo("alignReads.R/mergeBAMsAcrossDirs: done")
}

createSummaryAlignment <- function(save_dir, prepend_str) {
  loginfo("alignReads.R/createSummaryAlignment: counting unique bam reads...")
  
  ## count unique bam reads
  bamdir <- file.path(save_dir, "bams")
  bams <- dir(bamdir, full.names=TRUE)
  bams <- grep("bam$", bams, value=TRUE)
  uniq <- sapply(bams, function(bam) bamCountUniqueReads(bam))
  names(uniq) <- getBamType(names(uniq))

  ## save summary_alignment
  df <- data.frame(name=names(uniq), value=uniq)
  saveWithID(df, "summary_alignment", id=prepend_str, save_dir=file.path(save_dir, "results"), format="tab")

  ## compute stats on analyzed.bam
  bam <- getObjectFilename(bamdir, "analyzed.bam$")
  bstats <- computeBamStats(bam)
  df <- data.frame(name=names(bstats), value=bstats)
  saveWithID(df, "summary_analyzed_bamstats", id=prepend_str, save_dir=file.path(save_dir, "results"), format="tab")
}

##' Merge BAMs and create summary alignment file 
##'
##' @title Merge after alignReads
##' @param indirs A character vector, indicating which directories have to be merged
##' @param outdir A character string indicating the output directory (which must exist)
##' @param prepend_str A character string, containing a prefix going to be appended on all output result files
##' @param num_cores Number of cores available for parallel processing (for the merge bam step)
##' @return Nothing
##' @author Gregoire Pau
##' @export
##' @keywords internal
mergeAlignReads <- function(indirs, outdir, prepend_str, num_cores) {
  ## merge bams
  mergeBAMsAcrossDirs(indirs, outdir, prepend_str, nb.parallel.jobs=num_cores)
  
  ## merge summary alignments
  mergeSummaryAlignment(indirs, outdir, prepend_str)
}

##' Merge summary alignments
##'
##' @title Merge summary alignments
##' @param indirs A character vector, indicating which directories have to be merged
##' @param outdir A character string indicating the output directory (which must exist)
##' @param prepend_str A character string, containing a prefix going to be appended on all output result files
##' @return Nothing
##' @author Gregoire Pau
##' @keywords internal
##' @export
##' @keywords internal
mergeSummaryAlignment <- function(indirs, outdir, prepend_str) {
  ## merge summary alignment
  summary_alignment <- lapply(indirs, getNumericVectorDataFromFile, object_name="summary_alignment")
  summary_alignment <- as.data.frame(do.call(rbind, summary_alignment))
  summary_alignment <- apply(summary_alignment, 2, sum)
  df <- data.frame(name=names(summary_alignment), value=unlist(summary_alignment))
  saveWithID(df, "summary_alignment", id=prepend_str, save_dir=file.path(outdir, "results"), format="tab")

  ## merge summary_analyzed_bamstats
  summary_analyzed_bamstats <- lapply(indirs, getNumericVectorDataFromFile, object_name="summary_analyzed_bamstats")
  summary_analyzed_bamstats <- as.data.frame(do.call(rbind, summary_analyzed_bamstats))
  summary_analyzed_bamstats <- apply(summary_analyzed_bamstats, 2, sum)
  df <- data.frame(name=names(summary_analyzed_bamstats), value=unlist(summary_analyzed_bamstats))
  saveWithID(df, "summary_analyzed_bamstats", id=prepend_str, save_dir=file.path(outdir, "results"), format="tab")
  
  summary_alignment
}
