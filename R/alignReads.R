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
##' @return List of alignment files in BAM format
##' @keywords internal
alignReadsChunk <- function(fp1, fp2=NULL, save_dir=NULL) {
  ## init
  if (missing(save_dir)) save_dir <- file.path(getConfig("save_dir"))
  setUpDirs(save_dir, overwrite="overwrite")
  loginfo(paste("alignReads.R/alignReadsChunk: running gsnap..."))
          
  ## get config parameters
  prepend_str <- getConfig("prepend_str")
  bam_dir <- file.path(save_dir, "bams")
  
  ## align reads
  gsnapParams <- buildGsnapParams()
  wrapGsnap(fp1, fp2, gsnapParams, bam_dir, prepend_str)

  ## build analyzed bam
  buildAnalyzedBam(bam_dir, prepend_str)
  
  ## create summary
  createSummaryAlignment(save_dir, prepend_str)
  loginfo("alignReads.R/alignReadsChunk: done")
}

buildGsnapParams <- function() {
  ## get config parameters
  quality_encoding <- getConfig("quality_encoding")
  path.gsnap_genomes <- getConfig("path.gsnap_genomes")
  genome <- getConfig("alignReads.genome")
  snp_index <- getConfig("alignReads.snp_index")
  splice_index <- getConfig("alignReads.splice_index") 
  max_mismatches <- getConfig.integer("alignReads.max_mismatches")
  static_parameters <- getConfig("alignReads.static_parameters")
  nbthreads_perchunk <- getConfig.integer("alignReads.nbthreads_perchunk")
  sam_id <- getConfig("alignReads.sam_id")
  
  ## convert quality_encoding to gsnap quality protocol
  quality_protocol <- switch(quality_encoding,
                             sanger="sanger",
                             solexa="illumina",
                             illumina1.3="illumina",
                             illumina1.5="illumina",
                             illumina1.8="sanger",
                             'GATK-rescaled'="sanger")

  ## build core params
  params <- paste("-D", path.gsnap_genomes,
                  "-t", nbthreads_perchunk,
                  "-d", genome,
                  paste("--quality-protocol", quality_protocol, sep="="),
                  static_parameters,
                  "-A sam")

  ## extra params
  if (!is.null(sam_id)) params <- paste(params, " --read-group-id=", sam_id, sep="")
  if (!is.null(snp_index)) params <- paste(params, "-v", snp_index)
  if (!is.null(splice_index)) params <- paste(params, "-s", splice_index)
  if (!is.null(max_mismatches)) params <- paste(params, "-m", max_mismatches)
  
  return(params)
}

## create the analyzed bam
buildAnalyzedBam <- function(bam_dir, prepend_str) {
  ## get config parameters
  paired_ends <- getConfig.logical("paired_ends")
  alignReads.analyzedBam <- getConfig("alignReads.analyzedBam")

  ## prepare analyzed_bamregexp
  if (alignReads.analyzedBam=="uniq") analyzed_bamregexp <- "_uniq.*\\.bam$"
  else {
    if (paired_ends) analyzed_bamregexp <- "concordant_uniq.*\\.bam$"
    else analyzed_bamregexp <- "_uniq.*\\.bam$"
  }
  
  ## get uniq bams
  bam_files <- dir(bam_dir, full.names=TRUE)
  inbams <- grep(analyzed_bamregexp, bam_files, value=TRUE)
  
  ## merge inbams
  outbam <- paste(prepend_str, ".analyzed.bam", sep="")
  outbam <- file.path(bam_dir, outbam)
  mergeBams(inbams, outbam, sort=FALSE)
}

## given a bam file, return the bam type
getBamType <- function(bams) {
  bamtypes <- unname(sapply(basename(bams), function(bam) {
    pieces <- unlist(strsplit(bam, "\\."))
    paste(tail(pieces, n=2), collapse=".")
  }))
  gsub("\\.bam$", "", bamtypes)
}

## merges similarly named BAM residing in different directories and
## creates the BAI index files.
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
  z <- mclapply(seq_len(length(sinbams)),
                function(i) {
                  inbams <- sinbams[[i]]
                  outbam <- outbams[i]
                  mergeBams(inbams, outbam, sort=FALSE)
                },
                mc.cores=nb.parallel.jobs)
  
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
##' @keywords internal
##' @export
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
