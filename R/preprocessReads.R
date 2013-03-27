##' Pipeline preprocessing
##'
##' The preprocessing for our NGS pipelines
##' consists of :
##'  - quality filtering
##'  - check for adapter contamination
##'  - filtering of rRNA reads
##'  - read trimming
##'  - shortRead report generation of surviving reads
##'
##' These steps are mostly controlled by the global config.
##' @return A named vector containg the path to the preproceesed
##'         FastQ files and a few other statistics
##' @keywords internal
##' @export
preprocessReads <- function() {
  safeExecute({
    loginfo("preprocessReads.R/preprocessReads: starting...")
    
    ## initialise streamer
    FastQStreamer.init(input_file=getConfig("input_file"), input_file2=getConfig("input_file2"),
                       chunk_size=getConfig.integer("chunk_size"),
                       subsample_nbreads=getConfig.integer("subsample_nbreads"),
                       max_nbchunks=getConfig.integer("max_nbchunks"))
    
    ## preprocess reads using preprocessReadsChunk
    processChunks(FastQStreamer.getReads, preprocessReadsChunk)
    
    ## release streamer
    FastQStreamer.release()
    
    ## merge chunks
    chunkdirs <- getChunkDirs()
    save_dir <- getConfig("save_dir")
    prepend_str <- getConfig("prepend_str")
    summary_preprocess <- mergePreprocessReads(chunkdirs, save_dir, prepend_str)
    
    loginfo("preprocessReads.R/preprocessReads: done")
    invisible(summary_preprocess)
  }, memtracer=getConfig.logical("debug.tracemem"))
}

preprocessReadsChunk <- function(lreads, save_dir=NULL) {
  ## init
  if (missing(save_dir)) save_dir <- file.path(getConfig("save_dir"))
  setUpDirs(save_dir, overwrite="overwrite")
  
  ## determine read length and num mismatches
  total_reads <- length(lreads[[1]])
  if (total_reads>0) read_length <- width(sread(lreads[[1]])[1])
  else read_length <- NA
  summary_preprocess <- list(total_reads=total_reads, read_length=read_length)
  loginfo(paste("preprocessReads.R/preprocessReadsChunk: processing nbreads=", total_reads, "width=", read_length))

  ## compute input read min/max length
  z <- unlist(lapply(lreads, function(z) width(sread(z))))
  if (length(z)>0) summary_preprocess <- c(summary_preprocess, input_min_read_length=min(z), input_max_read_length=max(z))
  else summary_preprocess <- c(summary_preprocess, input_min_read_length=NA, input_max_read_length=NA)
  
  ## trim reads
  if (getConfig("trimReads.do")) {
    trim_len <- getConfig.integer("trimReads.length")
    lreads <- trimReads(lreads, trim_len)
  }
  
  ## filter reads
  if (getConfig.logical("filterQuality.do")) {
    lreads <- filterQuality(lreads)
    summary_preprocess$highqual_reads <- length(lreads[[1]])
  } else summary_preprocess$highqual_reads <- 0

  ## detect adapter contamination
  if(getConfig.logical("detectAdapterContam.do")){
    adapter_contaminated <- detectAdapterContam(lreads, save_dir=save_dir)
    summary_preprocess$adapter_contam <- sum(adapter_contaminated)
  } else summary_preprocess$adapter_contam <- 0
  
  ## detect rRNA contamination
  if (getConfig.logical("detectRRNA.do")) {
    is_rRNA_contam <- detectRRNA(lreads, save_dir=save_dir)
    lreads <- lapply(lreads, function(read) read[!is_rRNA_contam])    
    summary_preprocess$rRNA_contam_reads <- sum(is_rRNA_contam)
  } else summary_preprocess$rRNA_contam_reads <- 0

  ## compute processed read min/max length
  z <- unlist(lapply(lreads, function(z) width(sread(z))))
  if (length(z)>0) summary_preprocess <- c(summary_preprocess, processed_min_read_length=min(z), processed_max_read_length=max(z))
  else summary_preprocess <- c(summary_preprocess, processed_min_read_length=NA, processed_max_read_length=NA)
    
  ## write filtered and deduped fastq to run aligner on
  summary_preprocess$processed_reads <- length(lreads[[1]])
  writeFastQFiles(lreads, dir=file.path(save_dir, "bams"), filename1="processed.aligner_input_1.fastq",
                  filename2="processed.aligner_input_2.fastq")

  ## save summary_preprocess
  prepend_str <- getConfig("prepend_str")
  df <- data.frame(name=names(summary_preprocess), value=unlist(summary_preprocess))
  saveWithID(df, "summary_preprocess", id=prepend_str, save_dir=file.path(save_dir, "results"), format="tab")
  
  loginfo(paste("preprocessReads.R/preprocessReadsChunk: done"))
  return(save_dir)
}

##' Merge detectAdapterContam, merge preprocessed reads, create summary preprocess, build shortReadReport, remove processed 
##' 
##' @title Merge after preprocessReads
##' @param indirs A character vector, indicating which directories have to be merged
##' @param outdir  A character string indicating the output directory (which must exist)
##' @param prepend_str A character string, containing a prefix going to be appended on all output result files
##' @return Nothing
##' @author Gregoire Pau
##' @keywords internal
##' @export
mergePreprocessReads <- function(indirs, outdir, prepend_str) {
  ## get config parameters
  paired_ends <- getConfig.logical("paired_ends")
  detectAdapterContam.do <- getConfig.logical("detectAdapterContam.do")
  shortReadReport.do <- getConfig.logical("shortReadReport.do")
  remove_processedfastq <- getConfig.logical("remove_processedfastq")
  
  ## merge detectAdapterContam
  if (detectAdapterContam.do) mergeDetectAdapterContam(indirs, outdir, prepend_str, paired_ends)

  ## merge preprocess summary
  summary_preprocess <- mergeSummaryPreprocess(indirs, outdir, prepend_str)
  
  if (shortReadReport.do) {
    ## merge preprocessed reads
    mergePreprocessedFastQ(indirs, outdir, paired_ends)
    
    ## build ShortRead reports
    buildShortReadReports(outdir, paired_ends)
  }
  
  ## remove preprocessed reads
  if (remove_processedfastq) {
    pfilenames <- c("processed.aligner_input_1.fastq", "processed.aligner_input_2.fastq")
    pfilenames <- file.path(outdir, "bams", pfilenames)
    for (pfilename in pfilenames) {
      if (file.exists(pfilename)) unlink(pfilename)
    }
  }
  
  return(summary_preprocess)
}

mergePreprocessedFastQ <- function (indirs, outdir, paired_ends) {
  ## set the files to merge
  tfilenames <- c("processed.aligner_input_1.fastq")
  if (paired_ends) tfilenames <- c(tfilenames, "processed.aligner_input_2.fastq")

  for (tfilename in tfilenames) {
    ## guess filename from tfilename
    infiles <- sapply(indirs, function(indir) getObjectFilename(file.path(indir, "bams"), tfilename))

    ## merge files
    outfile <- file.path(outdir, "bams", tfilename)
    loginfo(paste("preprocessReads.R/mergePreprocessedReads: merging file=", outfile, "..."))
    system(paste("cat", paste(infiles, collapse=" "), ">", outfile))
  }
}

mergeSummaryPreprocess <- function(indirs, outdir, prepend_str) {
  ## read summary_preprocess
  summary_preprocess <- lapply(indirs, getNumericVectorDataFromFile, object_name="summary_preprocess")
  summary_preprocess <- as.data.frame(do.call(rbind, summary_preprocess))
  
  ## get read length
  read_length <- unique(na.omit(summary_preprocess$read_length))[1]
 
  ## merge
  summary_preprocess <- c(total_reads=sum(summary_preprocess$total_reads), 
                          highqual_reads=sum(summary_preprocess$highqual_reads),
                          adapter_contam=sum(summary_preprocess$adapter_contam),
                          read_length=read_length,
                          rRNA_contam_reads=sum(summary_preprocess$rRNA_contam_reads),
                          processed_reads=sum(summary_preprocess$processed_reads),
                          input_min_read_length=min(summary_preprocess$input_min_read_length, na.rm=TRUE),
                          input_max_read_length=max(summary_preprocess$input_max_read_length, na.rm=TRUE),
                          processed_min_read_length=min(summary_preprocess$processed_min_read_length, na.rm=TRUE),
                          processed_max_read_length=max(summary_preprocess$processed_max_read_length, na.rm=TRUE)
                          )
  summary_preprocess <- as.list(summary_preprocess)
  
  ## log merged summary_preprocess
  msg <- paste(paste(names(summary_preprocess), summary_preprocess, sep="="), collapse=" ")
  loginfo(paste("preprocessReads.R/mergeSummaryPreprocess:", msg))

  ## save merged summary_preprocess
  df <- data.frame(name=names(summary_preprocess), value=unlist(summary_preprocess), stringsAsFactors=FALSE)
  saveWithID(df, "summary_preprocess", id=prepend_str, save_dir=file.path(outdir, "results"), format="tab")
  
  invisible(summary_preprocess)
}

buildShortReadReports <- function(save_dir, paired_ends) {
  ## set the files to merge
  objects <- list(c(filename="processed.aligner_input_1.fastq", report_dir="shortReadReport_1"))
  if (paired_ends) {
    objects <- c(objects, list(c(filename="processed.aligner_input_2.fastq", report_dir="shortReadReport_2")))
  }
  
  for (object in objects) {
    filename <- file.path(save_dir, "bams", object["filename"])
    report_dir <- file.path(save_dir, "reports", object["report_dir"])
    if (file.exists(filename)) {
      loginfo(paste("preprocessReads.R/buildShortReadReports: generating report_dir=", report_dir, "..."))
      subsample_nbreads <-  getConfig.integer("shortReadReport.subsample_nbreads")
      if (is.null(subsample_nbreads)) {
        ## read all reads in memory
        reads <- readFastq(filename)
      } else {
        ## take a random subset of reads
        ## save and set fixed random seed, to allow reproducible behaviors (the seed is passed to FastqSampler)
        z = runif(1) ; originalseed <- .Random.seed ; set.seed(1)
      
        ## get reads
        fqs <- FastqSampler(filename, n=subsample_nbreads)
        reads <- yield(fqs)
        close(fqs)
        
        ## restore seed
        set.seed(originalseed)
      }
      
      ## generate report
      qao <- qa(reads, lane=filename)
      report(qao, dest=report_dir, type="html")
    } else {
      logwarn(paste("preprocessReads.R/buildShortReadReports: I can't generate the ShortReadReport since '", filename,
                    "' is not present", sep=""))
    }
  }
}
