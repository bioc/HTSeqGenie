##' Check configuration
##'
##' Performs all configuration checks
##' 
##' @return Nothing. Individual checks will throw error instead.
##' @keywords internal
##' @export
checkConfig <- function() {
  ## check config parts
  checkConfig.template()
  checkConfig.input()
  checkConfig.output()
  checkConfig.system()
  checkConfig.debug()
  checkConfig.trimReads()
  checkConfig.filterQuality()
  checkConfig.detectAdapterContam()
  checkConfig.detectRRNA()
  checkConfig.shortReadReport()
  checkConfig.alignReads()
  checkConfig.countGenomicFeatures()
  checkConfig.coverage()
  checkConfig.markDuplicates()
  checkConfig.realignIndels()
  checkConfig.analyzeVariants()
  checkConfig.noextraparameters()
  
  ## tag current version
  updateConfig(list(version=sessionInfo()$otherPkgs$HTSeqGenie$Version))
  
  ## everything OK
  invisible(TRUE)
}

## read template configuration file (and recursively if needed)
checkConfig.template <- function() {
  template_config <- getConfig("template_config")
  
  if (!is.null(template_config)) {
    ff <- findTemplate(template_config)
    ## load ff
    if (!is.null(ff)) {
      cat(paste("checkConfig.R/checkConfig.template: loading template config=", ff), "\n")
      updateConfig(list(template_config="")) ## remove template_config parameter from current config
      current_config <- getConfig() ## save current config
      loadConfig(ff) ## load template_config
      checkConfig.template() ## update config recursively
      updateConfig(current_config)      
    } else stop("checkConfig.R/checkConfig.template: cannot open file indicated by config parameter template_config=", getConfig("template_config"))
  }
}


findTemplate <- function(template_config){
  potential.config.dirs <- list("", "./", ## allows absolute and relative paths
                                strsplit(Sys.getenv("HTSEQGENIE_CONFIG"), ':')[[1]],
                                ## this will pull in inst/config while interactive
                                ## even if we are in a wrapper package that also defines a inst/config
                                getPackageFile("config"),
                                ## fallback is always the installed version of HTSeqGenie
                                system.file("config", package="HTSeqGenie"))

  paths <- file.path(unlist(potential.config.dirs), template_config)
  ff <- paths[file.exists(paths)][1]
  if (is.na(ff)) ff <- NULL
  return(ff)
}

checkConfig.input <- function() {
  ## input file
  file1 <- getConfig("input_file", stop.ifempty=TRUE)
  if (regexpr(" ", file1)>0) stop("'input_file' contains a whitespace: ", file1)
  if(! file.exists(file1)) stop("Input file 1 does not exist:", file1)
  
  ## input file 2
  paired_ends <- getConfig.logical("paired_ends", stop.ifempty=TRUE)
  if (paired_ends){
    file2 <- getConfig("input_file2", stop.ifempty=TRUE)
    if (regexpr(" ", file2)>0) stop("'input_file2' contains a whitespace: ", file2)
    if(! file.exists(file2)) stop("Input file 2 does not exist:", file2)
    if(file2==file1) stop("Input files are identical. Check your config file!")
  }
  else {
    p <- getConfig("input_file2")
    if (!is.null(p)) stop("'paired_ends' is FALSE but input_file2=", p)
  }

  ## detect quality
  cquals <- detectQualityInFASTQFile(file1)
  if (is.null(cquals)) stop("cannot open/detect quality of FASTQ file=", file1)
  
  ## set quality_encoding if unspecified
  quality_encoding <- getConfig("quality_encoding")
  if (is.null(quality_encoding)) {
    cat("possible qualities of filename=", file1, " are: ", paste(cquals, collapse=", "), "\n", sep="")
    quality_encoding <- cquals[1]
    cat("quality_encoding is not set! setting quality_encoding to", quality_encoding, "\n")
    updateConfig(list(quality_encoding=quality_encoding))
  }
  
  ## check if quality_encoding is compatible with the detected quality
  if (!(quality_encoding %in% cquals)) {
    stop("config parameter 'quality_encoding' is set to [", quality_encoding, "] but possible qualities of filename=", file1, " are [", paste(cquals, collapse=", "), "]\n")
  }
  
  getConfig.integer("chunk_size", stop.ifempty=TRUE) 
  getConfig.integer("subsample_nbreads")
  getConfig.integer("max_nbchunks")
}

checkConfig.output <- function() {
  ## save_dir
  save_dir <- getConfig("save_dir", stop.ifempty=TRUE)
  if (regexpr(" ", save_dir)>0) stop("'save_dir' contains a whitespace: ", save_dir)
  
  overwrite_save_dir <- getConfig("overwrite_save_dir", stop.ifempty=TRUE)
  if (!(overwrite_save_dir %in% c("never", "erase","overwrite"))) {
    stop("config parameter 'overwrite_save_dir' must be either 'never', 'erase' or 'overwrite'")
  }
  if (overwrite_save_dir=="never" && file.exists(save_dir)) stop("I won't overwrite data present in save_dir=", save_dir)

  ## prepend_str
  prepend_str <- getConfig("prepend_str", stop.ifempty=TRUE)
  if (regexpr(" ", prepend_str)>0) stop("'prepend_str' contains a whitespace: ", prepend_str)

  getConfig.logical("remove_processedfastq", stop.ifempty=TRUE)
  getConfig.logical("remove_chunkdir", stop.ifempty=TRUE)

  ## tmp_dir
  tmp.dir = getConfig('tmp_dir')
  if(!(is.null(tmp.dir))){
    if(!file.exists(tmp.dir)){
      stop("tmp_dir does not exists")
    }
  }
  ## chunk.dir is auto set from tmp_dir or save_dir
  if(is.null(getConfig('chunk_dir'))){
    setChunkDir()
  }
  getConfig("chunk_dir", stop.ifempty=TRUE)
}

checkConfig.system <- function() {
  ## check num_cores
  p <- "num_cores"
  a <- getConfig.integer(p)
  if (a<=0) stop("config parameter '", p, "' must be a non-zero positive integer")
}

checkConfig.debug <- function() {
  getConfig.logical("debug.tracemem", stop.ifempty=TRUE)
  level <- getConfig("debug.level", stop.ifempty=TRUE)
  if (!(level %in% c("ERROR", "WARN","INFO", "DEBUG"))) {
    stop("config parameter 'debug.level' must be either ERROR, WARN, INFO or DEBUG")
  }
}

checkConfig.trimReads <- function() {
  trimReads.do <- getConfig.logical("trimReads.do", stop.ifempty=TRUE)
  if (trimReads.do) {
    p <- "trimReads.length"
    a <- getConfig.integer(p, stop.ifempty=TRUE)
    if (a<=0) stop("config parameter '", p, "' must be a non-zero positive integer")
  }
}

checkConfig.filterQuality <- function() {
  filterQuality.do <- getConfig.logical("filterQuality.do", stop.ifempty=TRUE)
  if (filterQuality.do) {
    p <- "filterQuality.minQuality"
    a <- getConfig.integer(p, stop.ifempty=TRUE)
    if (a<0 || a>40) stop("config parameter '", p, "' must be an integer between 0 and 40")
    
    p <- "filterQuality.minFrac"
    a <- getConfig.numeric(p, stop.ifempty=TRUE)
    if (a<0 || a>1) stop("config parameter '", p, "' must be a number between 0 and 1")

    ## filterQuality.minLength
    filterQuality.minLength <- getConfig.numeric("filterQuality.minLength")
    if (!is.null(filterQuality.minLength) && filterQuality.minLength<0) {
      stop("checkConfig.R/checkConfig.filterQuality: 'filterQuality.minLength' must be a positive number")
    }
  }
}

checkConfig.detectAdapterContam <- function() {
  getConfig.logical("detectAdapterContam.do", stop.ifempty=TRUE)
  getConfig.logical("detectAdapterContam.force_paired_end_adapter", stop.ifempty=TRUE)
}

checkConfig.detectRRNA <- function() {
  detectRRNA.do <- getConfig.logical("detectRRNA.do", stop.ifempty=TRUE)
  if (detectRRNA.do) {    
    ## check if the genome exists
    p <- "detectRRNA.rrna_genome"
    a <- getConfig("path.gsnap_genomes", stop.ifempty=TRUE)
    genome <- getConfig(p, stop.ifempty=TRUE)
    dirname <- file.path(a, genome)
    if (!file.exists(dirname)) stop("cannot find the gsnap genome indicated by '", p, "' at=", dirname) 
  }
}

checkConfig.shortReadReport <- function() {
  shortReadReport.do <- getConfig.logical("shortReadReport.do", stop.ifempty=TRUE)
  if (shortReadReport.do) {
    a <- getConfig.integer("shortReadReport.subsample_nbreads")
    if (!is.null(a) && a<0) stop("shortReadReport.subsample_nbreads must be a positive integer")
  }
}

checkConfig.alignReads <- function() {
  ## alignReads.genome
  genome <- getConfig("alignReads.genome", stop.ifempty=TRUE)
  
  ## check if the genome exists
  p <- "alignReads.genome"
  a <- getConfig("path.gsnap_genomes", stop.ifempty=TRUE)
  genome <- getConfig(p, stop.ifempty=TRUE)
  dirname <- file.path(a, genome)
  if (!file.exists(dirname)) stop("cannot find the gsnap genome indicated by '", p, "' at=", dirname) 
  
  ## alignReads.max_mismatches
  max_mismatches <- getConfig.integer("alignReads.max_mismatches")
  if (!is.null(max_mismatches) && max_mismatches<0) stop("config parameter 'alignReads.max_mismatches' must be a positive integer")
  
  ## alignReads.sam_id
  sam_id <- getConfig("alignReads.sam_id")
  if (!is.null(sam_id)) {
    if (regexpr(" ", sam_id)>0) stop("'alignReads.sam_id' contains a whitespace: ", sam_id)
  }

  ## alignReads.nbthreads_perchunk
  num_cores <- getConfig.integer("num_cores", stop.ifempty=TRUE)
  nbthreads_perchunk <- getConfig.integer("alignReads.nbthreads_perchunk")
  if (is.null(nbthreads_perchunk)) {
    nbthreads_perchunk <- min(4, num_cores) ## set to min(4, num_cores) if empty
    updateConfig(list(alignReads.nbthreads_perchunk=nbthreads_perchunk)) 
  }
  if (nbthreads_perchunk<1 || nbthreads_perchunk>num_cores) stop("config parameter 'alignReads.nbthreads_perchunk' must be within (1;num_cores)")
  
  ## optional
  getConfig("alignReads.snp_index")
  getConfig("alignReads.splice_index")
  getConfig("alignReads.static_parameters")

  ## alignReads.analyzedBam
  alignReads.analyzedBam <- getConfig("alignReads.analyzedBam", stop.ifempty=TRUE)
  if (!(alignReads.analyzedBam %in% c("uniq", "concordant_uniq"))) {
    stop("config parameter 'alignReads.analyzedBam' must be either 'uniq' or 'concordant_uniq'")
  }

  ## alignReads.use_gmapR_gsnap
  alignReads.use_gmapR_gsnap <- getConfig.logical("alignReads.use_gmapR_gsnap", stop.ifempty=TRUE)
  if (!alignReads.use_gmapR_gsnap) {
    isgsnap <- system("gsnap", ignore.stdout=TRUE, ignore.stderr=TRUE)!=127
    if (!isgsnap) stop("'alignReads.use_gmapR_gsnap' is FALSE but 'gsnap' is not found in PATH")
  }

  invisible()
}

checkConfig.noextraparameters <- function() {
  ## read default config
  dconfig <- parseDCF(getPackageFile(file.path("config", "default-config.txt")))

  ## only parameters present in the default configuration file are allowed
  allowed_parameters <- names(dconfig)
  current_parameters <- names(getConfig())
  z <- setdiff(current_parameters, allowed_parameters)
  if (length(z)>0) {
    stop("unsupported config parameter(s): ", paste(z, collapse=", "))
  }

  ## check that no parameter is duplicated
  current_parameters <- names(getConfig())
  if (any(duplicated(current_parameters))) {
    stop("the following config parameters are duplicated: ", paste(current_parameters[duplicated(current_parameters)], collapse=", "))
  }
    
  invisible()
}

checkConfig.countGenomicFeatures <- function () {
  ## countGenomicFeatures.do
  countGenomicFeatures.do <- getConfig.logical("countGenomicFeatures.do", stop.ifempty=TRUE)
  if (countGenomicFeatures.do) {
    gpath <- getConfig("path.genomic_features", stop.ifempty=TRUE)
    gfeatures <- getConfig("countGenomicFeatures.gfeatures", stop.ifempty=TRUE)
  
    ## test existence of path.genomic_features/countGenomicFeatures.gfeatures
    filename <- file.path(gpath, gfeatures)
    if (!file.exists(filename)) stop("checkConfig.R/checkConfig.countGenomicFeatures: cannot open path.genomic_features/countGenomicFeatures.gfeatures=", filename)
  }
}

checkConfig.coverage <- function() {
  ## coverage.do
  coverage.do <- getConfig.logical("coverage.do", stop.ifempty=TRUE)
  if (coverage.do) {
    ## coverage.extendReads
    coverage.extendReads <- getConfig.logical("coverage.extendReads", stop.ifempty=TRUE)
    
    ## coverage.fragmentLength
    coverage.fragmentLength <- getConfig.numeric("coverage.fragmentLength")
    if (!is.null(coverage.fragmentLength) && !coverage.extendReads) {
      stop("checkConfig.R/checkConfig.coverage: 'coverage.extendReads' is FALSE but 'coverage.fragmentLength' is not empty")
    }
    
    ## coverage.maxFragmentLength
    coverage.maxFragmentLength <- getConfig.numeric("coverage.maxFragmentLength")
    if (!is.null(coverage.maxFragmentLength) && coverage.maxFragmentLength<0) {
      stop("checkConfig.R/checkConfig.coverage: 'coverage.maxFragmentLength' must be a positive number")
    }
  }
}

checkConfig.markDuplicates <- function() {
  markDuplicates.do <- getConfig.logical("markDuplicates.do", stop.ifempty=TRUE)
  if (markDuplicates.do) {    
    picard.path <- getConfig("path.picard_tools", stop.ifempty=TRUE)
    ## check if the tool is callable
    if (!checkPicardJar("MarkDuplicates", path=picard.path))
      stop("MarkDuplicates.jar not found or not executable")
  }
}

checkConfig.analyzeVariants <- function() {

  if (getConfig.logical("analyzeVariants.do", stop.ifempty=TRUE) ){  
    ## analyzeVariants.method
    analyzeVariants.method <- getConfig("analyzeVariants.method", stop.ifempty=TRUE)
    if (length(analyzeVariants.method)==0) stop("'analyzeVariants.method' must contain a character string")

    knownMethods <- c("VariantTools", "GATK")
    z <- analyzeVariants.method%in%knownMethods
    if (!z) stop("'analyzeVariants.method' must contain the name of a method among: ", paste(knownMethods, collapse=", "))
    if ("GATK" %in% analyzeVariants.method) {
      checkConfig.GATK()
      if(getConfig.logical("gatk.filter_repeats")){
        rmask.file <- getConfig("analyzeVariants.rep_mask")
        if( is.null(rmask.file) || !file.exists(rmask.file)){
          stop(paste("Repeat mask file either not set or not found at", rmask.file))
        }
      }
    }
    if ("VariantTools" %in% analyzeVariants.method) {
      ## analyzeVariants.bqual
      analyzeVariants.bqual <- getConfig.numeric("analyzeVariants.bqual", stop.ifempty=TRUE)
      if(analyzeVariants.bqual < 0 || analyzeVariants.bqual > 41)
        stop("analyzeVariants.bqual must be an integer between 0 and 41")     
      
      ## analyzeVariants.rep_mask
      rmask.file <- getConfig("analyzeVariants.rep_mask")
      if(!(is.null(rmask.file))){
        if(!file.exists(rmask.file)){
          stop(paste("Repeat mask file not found at", rmask.file))
        }
      }
      ## analyzeVariants.dbsnp
      dbsnp <- getConfig("analyzeVariants.dbsnp", stop.ifempty=FALSE)
      if(!(is.null(dbsnp))){
        if(!file.exists(dbsnp)){
          stop(paste("dbsnp Rdata file not found at", dbsnp))
        }
      }
      ## analyzeVariants.callingFilters
      calling.filters <- getConfig.vector("analyzeVariants.callingFilters")
      all.filters = VariantCallingFilters()
      if (! all(calling.filters %in% names(all.filters))) {
        stop("Some of the VariantCallingFilters are not provided by VariantTools!")
      }
      ## analyzeVariants.postFilters
      post.filters <- getConfig.vector("analyzeVariants.postFilters")
      all.filters = VariantPostFilters()
      if (! all(post.filters %in% names(all.filters))) {
        stop("Some of the VariantPostFilters are not provided by VariantTools!")
      }
    }
    ## analyzeVariants.indels
    getConfig.logical("analyzeVariants.indels", stop.ifempty=TRUE)
  }
  invisible(TRUE)
}


## check for GATK and presence of a genome file
checkConfig.GATK <- function() {
  gatk.path <- getConfig("path.gatk", stop.ifempty=TRUE)
  if ( checkGATKJar(gatk.path) != TRUE ){
    stop("Can't run indelRealigner without GATK! Check your path.gatk settings and make sure you have java in your PATH")
  }
  gatk.genomes <- getConfig("path.gatk_genomes", stop.ifempty=TRUE)
  genome <- getConfig("alignReads.genome", stop.ifempty=TRUE)
  genome.path <- file.path(gatk.genomes, paste0(genome, ".fa"))
  if(!file.exists(genome.path)){
    stop(paste("GATK genome not found at", genome.path))
  }
}

checkConfig.realignIndels <- function() {
  if (getConfig.logical("realignIndels.do", stop.ifempty=TRUE) ){  
    checkConfig.GATK()
  }
}
