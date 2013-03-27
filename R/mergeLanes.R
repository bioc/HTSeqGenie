##' Merge input lanes built by the NGS pipeline
##'
##' @title Merge input lanes
##' @param indirs A character vector of directory paths containing NGS pipeline output
##' @param outdir A character string pointing to a non-existing output directory
##' @param prepend_str A character string, containing a prefix going to be appended on all output result files
##' @param num_cores Number of cores available for parallel processing (for the merge bam step)
##' @param config_update List of name value pairs that will update the config parameters
##' @param preMergeChecks.do A logical, indicating whether to perform pre merge checks
##' @param ignoreConfigParameters A character vector containing the configuration parameters that are not required to be identical
##' @return Nothing
##' @author greg
##' @keywords internal
##' @export
mergeLanes <- function(indirs, outdir, prepend_str, num_cores, config_update, preMergeChecks.do=TRUE, ignoreConfigParameters) {
  ## set up logger
  addHandler(writeToConsole, level="DEBUG")

  ## check ignoreConfigParameters
  if (missing(ignoreConfigParameters)) {
    ignoreConfigParameters <- c("input_file", "input_file2", "subsample_nbreads", "quality_encoding",
                                "chunk_size", "max_nbchunks",
                                "save_dir", "overwrite_save_dir", "prepend_str", "remove_processedfastq",
                                "remove_chunkdir",
                                "num_cores",
                                "debug.level", "debug.tracemem",
                                "trimReads.do", "trimReads.length",
                                "shortReadReport.do", "shortReadReport.subsample_nbreads",
                                "alignReads.sam_id", "alignReads.nbthreads_perchunk",
                                
                                ## additions in 3.6.0
                                "tmp_dir", "chunk_dir", "countGenomicFeatures.do"
                                )
  }
  
  ## check input configs and build output config
  config <- checkInputConfigs(indirs, ignoreConfigParameters=ignoreConfigParameters)
  loadConfig()
  updateConfig(config)
  updateConfig(list(save_dir=outdir, prepend_str=prepend_str, debug.level="DEBUG", overwrite_save_dir="never", num_cores=num_cores))
  if (!missing(config_update)) updateConfig(config_update)
 
  ## set up outdir and logger
  initDirs()
  initLogger()

  ## heavy lifting
  safeExecute({
    ## check lanes are mergeable
    if (preMergeChecks.do) preMergeChecks(indirs)
  
    ## merge lanes
    doMergeLanes(indirs)

    ## check output
    postMergeChecks(indirs, outdir)
    
  }, memtracer=TRUE)
}

doMergeLanes <- function(indirs) {
  ## get parameters
  outdir <- getConfig("save_dir")
  prepend_str <- getConfig("prepend_str")
  num_cores <- getConfig.integer("num_cores")
  loginfo("mergeLanes.R/doMergeLanes: starting...")
  
  mergePreprocessReads(indirs, outdir, prepend_str)
  mergeAlignReads(indirs, outdir, prepend_str, num_cores)
  writePreprocessAlignReport()
  if (getConfig.logical("countGenomicFeatures.do")) {
    safeExecute({mergeCountGenomicFeatures(indirs, outdir, prepend_str)}, memtracer=FALSE)
  }
  if (getConfig.logical("coverage.do")) {
    safeExecute({mergeCoverage(indirs, outdir, prepend_str)}, memtracer=FALSE)
  }
  if(getConfig.logical("analyzeVariants.do")){
    safeExecute({analyzeVariants()}, memtracer=FALSE)
  }
  
  loginfo("mergeLanes.R/doMergeLanes: merge lanes successful.")
}

preMergeChecks <- function(indirs) {
  loginfo("mergeLanes.R/preMergeChecks: checking indirs...")

  summaries <- parseSummaries(indirs, "summary_preprocess")
  max_read_lengths <- summaries[,'processed_max_read_length']
  min_read_lengths <- summaries[,'processed_min_read_length']

  if (length(unique(max_read_lengths))!=1) {
    stop("mergeLanes.R/preMergeChecks: max read lengths of input dirs are not identical [",
         paste(paste(indirs, ": ", max_read_lengths, sep=""), collapse=", "), "]")
  }

  if (length(unique(min_read_lengths))!=1) {
    stop("mergeLanes.R/preMergeChecks: min read lengths of input dirs are not identical [",
         paste(paste(indirs, ": ", min_read_lengths, sep=""), collapse=", "), "]")
  }
  
  if(getConfig.logical("analyzeVariants.do")){
    .IDverify(indirs)
  } else {
    logwarn("Variant concordance between samples has not been checked, as variants were not called for inputs!");
  }

  loginfo("mergeLanes.R/preMergeChecks: done")
}

.IDverify <- function(indirs, iivar_min=0.8){
  if (length(indirs)>1) {
    variants <- sapply(indirs, function(indir) getObjectFilename(file.path(indir, "results"), "filtered_variants_granges"))
    iivar <- rep(NA, length(indirs)-1)
    var1 <- get(load(variants[1]))
    nvar1 <- names(variants[1])
    for (i in 1:length(iivar)) {
      var2 <- get(load(variants[i+1]))
      nvar2 <- names(variants[i+1])
      iiv <- VariantTools:::checkVariantConcordance(var1, var2)
      loginfo(paste("mergeLanes.R/IDVerify: variantConcordance between [", nvar1, "] and [", nvar2, "] is a=", iiv[1], "b=", iiv[2]))
      iivar[i] <- iiv$fraction
      var1 <- var2
      nvar1 <- nvar2
    }
    if (any(is.na(iivar)) || min(iivar)<iivar_min) {
      stop("mergeLanes.R/IDverify: ID verify failed [", paste(iivar, collapse=", "), "], iivar_min=", iivar_min)
    }
  }
}

postMergeChecks <- function(indirs, outdir) {
}

checkInputConfigs <- function(indirs, ignoreConfigParameters=NULL) {
  ## load all configs
  config_filenames <- file.path(indirs, "logs", "config.txt")
  configs <- lapply(config_filenames, function(config_filename) {
    if (!file.exists(config_filename)) stop("mergeConfigs: cannot open config file=", config_filename)
    else unlist(parseDCF(config_filename))
  })
  names(configs) <- config_filenames

  ## check versions
  versions <- lapply(configs, function (config) config["version"])
  if (length(unique(versions))>1) {
    msg <- paste(paste(indirs, " (", versions, ")", sep=""), collapse=", ")
    stop("mergeConfigs: cannot merge samples analyzed with different versions: ", msg)
  }

  ## do all input configs contain the same parameters?
  pconfigs <- lapply(configs, names)
  allp <- unique(unlist(pconfigs))
  for (i in 1:length(pconfigs)) {
    z <- allp%in%intersect(pconfigs[[i]], allp)
    config_filename <- names(pconfigs)[i]
    if (any(!z)) stop("mergeConfigs: config parameters [", paste(allp[!z], collapse=", "), "] are missing from config file=", config_filename)
  }

  ## merge configs
  configs <- do.call(rbind, configs)
  oconfig <- sapply(colnames(configs), function(parameter) {
    mvalue <- unique(configs[,parameter])
    if (length(mvalue)>1) {
      mmvalue <- paste("<merged>", paste(mvalue, collapse=" "), "</merged>", sep="")
      if (parameter%in%ignoreConfigParameters) mvalue <- mmvalue
      else stop("mergeLanes.R/mergeConfig: parameter '", parameter, "' has different values [", paste(mvalue, collapse=", "), "]")
    }
    mvalue
  })
  
  oconfig
}

mergeCountGenomicFeatures <- function(indirs, outdir, prepend_str) {
  merge_nbfeatures <- mergeCounts(indirs, outdir, prepend_str)
  mergeSummaryCounts(indirs, outdir, prepend_str, merge_nbfeatures)   
  writeGenomicFeaturesReport()
}
