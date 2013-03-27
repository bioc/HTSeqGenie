##' Run the NGS analysis pipeline
##'
##' This function starts the pipeline. It first preprocesses the input FASTQ reads, align them, count the read overlaps with genomic features and
##' compute the coverage. See the vignette for details.
##'
##' @title  Run the NGS analysis pipeline
##' @param ... A list of parameters. See the vignette for details.
##' @return The path to the NGS output directory.
##' @author Jens Reeder, Gregoire Pau
##' @seealso TP53Genome, TP53GenomicFeatures
##' @examples \dontrun{
##' ## build genome and genomic features
##' tp53Genome <- TP53Genome()
##' tp53GenomicFeatures <- TP53GenomicFeatures()
##'  
##'  ## get the FASTQ files
##' fastq1 <- system.file("extdata/H1993_TP53_subset2500_1.fastq.gz", package="HTSeqGenie")
##' fastq2 <- system.file("extdata/H1993_TP53_subset2500_2.fastq.gz", package="HTSeqGenie")
##' 
##' ## run the pipeline
##' save_dir <- runPipeline(
##'     ## input
##'     input_file=fastq1,
##'     input_file2=fastq2,
##'     paired_ends=TRUE,
##'     quality_encoding="illumina1.8",
##'     
##'     ## output
##'     save_dir="test",
##'     prepend_str="test",
##'     overwrite_save_dir="erase",
##'     
##'     ## aligner
##'     path.gsnap_genomes=path(directory(tp53Genome)),
##'     alignReads.genome=genome(tp53Genome),
##'     alignReads.additional_parameters="--indel-penalty=1 --novelsplicing=1 --distant-splice-penalty=1",
##'     
##'     ## gene model
##'     path.genomic_features=dirname(tp53GenomicFeatures),
##'     countGenomicFeatures.gfeatures=basename(tp53GenomicFeatures)
##'     )
##'}
##' 
##' 
##' @export
runPipeline <- function(...) {
  runPipelineConfig(config_update=list(...))
}


##' Run the NGS analysis pipeline from a configuration file
##'
##' This is the launcher function for all
##' pipeline runs. It will do some
##' preprocessing steps, then aligns the reads,
##' counts overlap with genomic Features such as genes, exons etc
##' and applies a variant caller.
##'
##' @title  Run the NGS analysis pipeline
##' @param config_filename Path to a pipeline configuration file
##' @param config_update A list of name value pairs that will update the config parameters
##' @return Nothing
##' @author Jens Reeder, Gregoire Pau
##' @export
runPipelineConfig <- function(config_filename, config_update) {
  initPipelineFromConfig(config_filename, config_update)

  preprocessReads()
  alignReads()

  if (getConfig.logical("markDuplicates.do")) {
    markDups()
  }
  if (getConfig.logical("countGenomicFeatures.do")) {
    countGenomicFeatures()
  }
  if (getConfig.logical("coverage.do")) calculateCoverage()

  if (getConfig.logical("analyzeVariants.do")) analyzeVariants()

  removeChunkDir()

  loginfo("Pipeline run successful.")
  invisible(getConfig("save_dir"))
}

##' Remove chunk directories
##'
##' A pipeline run processes the data in small
##' chunks, which are eventually combined into the
##' final result. Afterwards, this function can be called to remove
##' the temporary results per chunk.
##' @title Remove chunk directories
##' @return Nothing 
##' @author Jens Reeder
##' @export
##' @keywords internal
removeChunkDir <- function() {
  ## get configuration parameters
  save_dir <- getConfig("save_dir")
  fst.chunkdir <- getChunkDirs()[1]

  ## copy the first chunk log
  fin <- file.path(fst.chunkdir, "logs", "progress.log")
  fout <- file.path(save_dir, "logs", "progress-chunk1.log")
  file.copy(fin, fout)
  
  if (getConfig.logical("remove_chunkdir")) {
    ## remove chunkdir
    if (is.null(getConfig('chunk_dir'))){
      stop("Can't remove chunk dirs, since config variable is not set")
    }
    unlink(getConfig('chunk_dir'), recursive=TRUE)
  }
}
