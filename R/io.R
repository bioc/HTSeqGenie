##' Read FastQ input files
##'
##' Uses the global config to find input files
##' @return Reads as list of ShortRead objects
##' @export
##' @keywords internal
readInputFiles <- function() {
  loginfo("io.R/readInputFiles: reading FASTQ file(s)...")

  ## get config parameters
  paired_ends <- getConfig.logical("paired_ends")
  filename1 <- getConfig("input_file")
  if (paired_ends) filename2 <- getConfig("input_file2")

  ## check file existence
  if (!file.exists(filename1)) stop("cannot open file indicated by config parameter input_file=", filename1)
  if (paired_ends) if (!file.exists(filename2)) stop("cannot open file indicated by config parameter input_file2=", filename2)

  ## read files
  lreads <- list(readFastq(filename1))
  loginfo(paste("io.R/readInputFiles: read", length(lreads[[1]]), "reads from input_file=", filename1))
  if (paired_ends) {
    lreads <- c(lreads, list(readFastq(filename2)))
    loginfo(paste("io.R/readInputFiles: read", length(lreads[[2]]), "reads from input_file2=", filename2))
  }
  
  loginfo("io.R/readInputFiles: done")
  return(lreads)
}

##' Write reads to file
##'
##' @param lreads List of reads as ShortRead objects
##' @param dir Save directory
##' @param filename1 Name of file 1
##' @param filename2 Name of file 2
##' @return Named list of filepaths
##' @export
##' @keywords internal
##' @importMethodsFrom ShortRead writeFastq
writeFastQFiles <- function (lreads, dir, filename1, filename2) {
  ## get config parameters
  paired_ends <- getConfig.logical("paired_ends")

  ## write forward reads
  filepaths <- list()
  filepaths$fastq_for_aligner_1 <- file.path(dir, filename1)
  loginfo(paste("io.R/writeFastQFiles: writing filename=", filepaths$fastq_for_aligner_1))
  if (length(lreads[[1]])>0) writeFastq(lreads[[1]], file=filepaths$fastq_for_aligner_1, lane="", compress=FALSE)
  else cat("", file=filepaths$fastq_for_aligner_1) ## empty FASTQ file 

  ## write reverse reads
  if (paired_ends) {
    filepaths$fastq_for_aligner_2 <- file.path(dir, filename2)
    loginfo(paste("io.R/writeFastQFiles: writing filename=", filepaths$fastq_for_aligner_2))
    if (length(lreads[[2]])>0) writeFastq(lreads[[2]], file=filepaths$fastq_for_aligner_2, lane="", compress=FALSE)
    else cat("", file=filepaths$fastq_for_aligner_2) ## empty FASTQ file 
  }

  return(filepaths)
}

##' Count reads in Fastq file
##'
##' @title Count reads in Fastq file
##' @param filename Name of FastQ file
##' @return Number of reads
##' @author Gregoire Pau
##' @export
##' @importMethodsFrom ShortRead FastqStreamer
##' @keywords internal
getNumberOfReadsInFASTQFile <- function(filename) {
  if (length(grep("\\.gz$", filename))>0) con <- gzfile(filename)
  else con <- file(filename)
  fqs <- FastqStreamer(con, n=1e6)
  n <- 0
  repeat {
    nreads <- length(safe.yield(fqs))
    if (nreads==0) break
    else n <- n + nreads
  }
  close(fqs)
  n
}

##' Open a streaming connection to a FastQ file
##'
##' Only one FastQStreamer object can be open at any time.
##' 
##' @title Open a streaming connection to a FastQ file
##' @param input_file Path to a FastQ file
##' @param input_file2 Optional path to a FastQ file. Default is NULL.
##' @param chunk_size Number of reads per chunk 
##' @param subsample_nbreads Optional number of reads to subsample (deterministic) from the input files. Default is NULL.
##' @param max_nbchunks Optional maximal number of chunks to read
##' @seealso FastQStreamer.getReads
##' @return Nothing. 
##' @author Gregoire Pau
##' @export
##' @keywords internal
FastQStreamer.init <- function(input_file, input_file2=NULL, chunk_size, subsample_nbreads=NULL, max_nbchunks=NULL) {
  ## check file existences
  if (!file.exists(input_file)) stop("cannot open file indicated by parameter input_file=", input_file)
  if (!is.null(input_file2)) {
    if (!file.exists(input_file2)) stop("cannot open file indicated by parameter input_file2=", input_file2)
  }
  
  ## build subsampling_filter if subsample_nbreads is not NULL
  subsampling_filter <- NULL
  if (!is.null(subsample_nbreads)) {    
    loginfo(paste("io.R/FastQStreamer.init: counting number of reads in file=", input_file))
    totnbreads <- getNumberOfReadsInFASTQFile(input_file)
    
    if (subsample_nbreads<totnbreads) {
      ## save and set fixed random seed, to allow reproducible behaviors
      runif(1) ; originalseed <- .Random.seed ; set.seed(1)

      ## build subsampling filter
      subsampling_filter <- rep(FALSE, totnbreads)
      subsampling_filter[sample(seq_len(totnbreads), subsample_nbreads)] <- TRUE
      
      ## restore seed
      set.seed(originalseed)
      loginfo(paste("io.R/FastQStreamer.init: subsampling_filter set (subsampled reads=", subsample_nbreads,
                    ", totnbreads=", totnbreads, ")"))
    } else {
      loginfo(paste("io.R/FastQStreamer.init: the requested number of subsampled reads=", subsample_nbreads,
                    "is higher than the total number of reads=", totnbreads, ": I won't subsample"))
    }
  } 
  
  ## initialise streamers
  filenames <- list(input_file)
  if (!is.null(input_file2)) filenames <- c(filenames, input_file2)
  lfqs <- lapply(filenames, function (filename) {
    if (length(grep("\\.gz$", filename))>0) con1 <- gzfile(filename)
    else con1 <- file(filename)
    fqs <- FastqStreamer(con1, n=chunk_size)
    loginfo(paste("io.R/FastQStreamer.init: initialised FastQ streamer for filename=", filename))
    fqs
  })
  
  FastQStreamer.lfqs <<- lfqs
  FastQStreamer.chunkid <<- 0
  FastQStreamer.subsampling_filter <<- subsampling_filter
  
  if (is.null(max_nbchunks)) FastQStreamer.max_nbchunks <<- Inf
  else FastQStreamer.max_nbchunks <<- as.integer(max_nbchunks)
}

##' Get FastQ reads from the FastQ streamer
##'
##' @title Get FastQ reads from the FastQ streamer
##' @return A list of ShortRead object containing reads. NULL if there are no more reads to read.
##' @author Gregoire Pau
##' @seealso FastQStreamer.init
##' @export
##' @keywords internal
FastQStreamer.getReads <- function() {
  ## check
  if (!exists("FastQStreamer.lfqs")) {
    FastQStreamer.lfqs <- NULL ## to avoid the 'no visible binding' R CMD check warning
    FastQStreamer.chunkid <- NULL ## to avoid the 'no visible binding' R CMD check warning
    FastQStreamer.max_nbchunks <- NULL ## to avoid the 'no visible binding' R CMD check warning
    FastQStreamer.subsampling_filter <- NULL ## to avoid the 'no visible binding' R CMD check warning
    stop("io.R/FastQStreamer.getReads: FastQStreamer.init() has not been called")
  }
  
  ## get reads
  lfqs <- FastQStreamer.lfqs
  paired_ends <- length(lfqs)==2
  lreads <- list(safe.yield(lfqs[[1]]))
  if (paired_ends) {
    lreads <- c(lreads, list(safe.yield(lfqs[[2]])))
    if (length(lreads[[1]])!=length(lreads[[2]]))
      stop("io.R/FastQStreamer.getReads: input files must have the same number of reads")
  }

  ## increment chunkid
  FastQStreamer.chunkid <<- FastQStreamer.chunkid + 1
  if (length(lreads[[1]])==0 || FastQStreamer.chunkid>FastQStreamer.max_nbchunks) lreads <- NULL
  else {
    ## subsample, if required
    if (!is.null(FastQStreamer.subsampling_filter)) {
      len <- length(lreads[[1]])
      if (len>0) {
        z <- FastQStreamer.subsampling_filter[1:len]
        ## if (sum(z)>0) {
          lreads <- lapply(lreads, function(reads) reads[z])
          FastQStreamer.subsampling_filter <<- FastQStreamer.subsampling_filter[-(1:len)]
        ## } else {
       ##   stop("io.R/FastQStreamer.getReads: (known bug) the chunk is empty after subsampling. I won't process an empty chunk! Please disable or increase 'subsample_nbreads', or 'increase chunk_size'")
        ##}
      }
    }
  }
  
  lreads
}

##' Close the FastQStreamer 
##'
##' @title Close the FastQStreamer 
##' @return Nothing
##' @seealso FastQStreamer.init
##' @author Gregoire Pau
##' @export
##' @keywords internal
FastQStreamer.release <- function() {
  ## check
  if (!exists("FastQStreamer.lfqs")) {
    FastQStreamer.lfqs <- NULL ## to avoid the 'no visible binding' R CMD check warning
    stop("io.R/FastQStreamer.release: FastQStreamer.init() has not been called")
  }
  
  invisible(lapply(FastQStreamer.lfqs, close))
}

##' Save an R object
##' 
##' Exists so objects can be serialized and reloaded with the a unique
##' identifier in the symbol. Stores the data object with a new nane
## 'created by appending 'orig_name' and 'id'
##' @param data The data to store
##' @param orig_name The original name of the data
##' @param id A meaningful id the is prepended to the stored objects name
##' @param save_dir The directory where the data should be saved in
##' @param compress Save the data compressed or not
##' @param format Choice of 'RData' or 'tab'(ular)
##' @return Name of the stored file
##' @export
##' @keywords internal
saveWithID <- function(data, orig_name, id, save_dir, compress=TRUE, format="RData") {
  new_name <- paste(id, orig_name, sep=".")
  
  if (format=="RData") {
    ## RData format
    assign(new_name, data)
    filename <- paste(new_name, ".RData", sep="")
    filename <- file.path(save_dir, filename)
    loginfo(paste("io.R/saveWithID: saving file=", filename))
    save(list=new_name, file=filename, compress=compress)
  } else if (format=="tab"){
    ## tab-separated format
    filename <- paste(new_name, ".tab", sep="")
    filename <- file.path(save_dir, filename)
    loginfo(paste("io.R/saveWithID: saving file=", filename))
    write.table(data, filename, row.names=FALSE, sep="\t", quote=FALSE)
  }
  else{
    stop("saveWithID: unknown format: ", format)
  }
  return(filename)
}

##' Make a directory after performing an existence check
##'
##' Throws an exception if file or directory with same name exist
##' and overwrite is TRUE.
##' @param dir Name of directory to create
##' @param overwrite A character string: never (default), erase, overwrite
##' @return Path to created directory
##' @export
##' @keywords internal
makeDir <- function(dir, overwrite="never"){
  if (file.exists(dir)) {
    if (overwrite=="never") stop("io.R/makeDir: I won't overwrite data present in dir=", dir)
    if (overwrite=="erase") safeUnlink(dir)
  }
  dir.create(dir, recursive=TRUE, showWarnings=FALSE)
  return(dir)
}

##' Create a random directory with prefix in R temp dir
##'
##' Especially for testing code it is very helpful to have
##' a temp directory with a defined prefix, so one knows which test
##' produced which directory.
##' @param prefix A string that will preceed the directory name
##' @param dir Directory where the random dir will be created under.
##'            Defaults to tempdir()
##' @return Name of temporary directory
##' @keywords internal
createTmpDir <- function(prefix=NULL, dir=tempdir()) {
  tmpname <- NULL
  while(is.null(tmpname)){
    tmpname <- try(
                   {
                     tmp <- file.path(dir, paste(c(prefix, sample(letters,8)), collapse=""))
                     makeDir(tmp, overwrite="never")
                     return(tmp)
                   }, silent=TRUE)
     if (class(tmpname)=="try-error") {
       tmpname <- NULL ## restarts while
     }
  }
  return(tmpname)
}

##' Creates a directory with all needed subdirectories for pipeline
##' outputs
##'
##' @title Create output directory and subdirectories for sequencing pipeline analysis outputs
##' @param save_dir path to the directory that will contain all
##' needed subdirectories
##' @param overwrite A character string: never (default), erase, overwrite
##' @return Nothing. Called for its side effects
##' @author Cory Barr, Jens Reeder
##' @export
##' @keywords internal
setUpDirs <- function(save_dir, overwrite="never") {
  makeDir(save_dir, overwrite=overwrite)
  dirs <- c("bams", "logs", "results", "reports", "reports/images")
  dirs <- file.path(save_dir, dirs)
  invisible(sapply(dirs, makeDir, overwrite=overwrite))
}

##' Safely load a R data file
##'
##' Attempts to load a file given by object_name.
##' Bails out if none or more than one files match the object name.
##' 
##' @param dir_path Save dir of a pipeline run
##' @param object_name object name, can be a regexp
##' @return loaded object
##' @export
##' @keywords internal
safeGetObject <- function(dir_path, object_name){
  filename <- getObjectFilename(file.path(dir_path,"results"), object_name)
  obj <- try({
    get(load(filename))
  }, silent=TRUE)

  if (class(obj)=="try-error") {
    stop("error: cannot load ", object_name, "in", dir_path, "\n", sep="")
  }
  return(obj)
}

##' Load tabular data from the NGS pipeline result directory
##' @param save_dir A character string containing an NGS pipeline output directory.
##' @param object_name A character string ontaining the regular expression matching a filename in dir_path
##' @return A data frame.
##' @export
getTabDataFromFile <- function(save_dir, object_name) {
  filename <- getObjectFilename(file.path(save_dir, "results"), paste(".", object_name, "\\.tab$", sep=""))
  read.table(filename, sep="\t", header=TRUE, stringsAsFactors=FALSE)
}

##' Load data as numerical values
##'
##' @param dir_path Save dir of a pipeline run
##' @param object_name Object name
##' @return loaded data as table of numbers
##' @author Jens Reeder
##' @export
##' @keywords internal
getNumericVectorDataFromFile <- function(dir_path, object_name) {
  df <- getTabDataFromFile(dir_path, object_name)
  v <- as.numeric(df$value)
  names(v) <- df$name
  return(v)
}

##' Detect quality protocol from a FASTQ file
##'
##' @param filename Path to a FASTQ or gzipped-FASTQ file
##' @param nreads Number of reads to test quality on. Default is 5000.
##' @return A character vector containing the compatible qualities. NULL if none.
##' @author Jens Reeder
##' @export
##' @importMethodsFrom IRanges as.vector
##' @keywords internal
detectQualityInFASTQFile <- function(filename, nreads=5000) {
  cquals <- try({
    if (length(grep("\\.gz$", filename))>0) con <- gzfile(filename)
    else con <- file(filename)
    
    ## read qualities
    fqs <- FastqStreamer(con, n=nreads)
    reads <- safe.yield(fqs)
    z <- paste(as.vector(quality(quality(reads))), collapse="")
    quals <- as.numeric(charToRaw(z))
    close(con)
    
    ## check
    if (length(quals)==0) stop()
      
    ## build a dictionary of known qualities
    ## "illumina1.5" upper bound is now 105 instead of 104, to accomodate Phred-qualities of 41 (instead of 40)
    ## some bam files from the TCGA project come with rescaled phred qualities from 1-50
    ## to accomodate for them we invent a new quality range "GATK-rescaled"
    knownquals <- data.frame(qual=c("solexa", "illumina1.3", "illumina1.5", "illumina1.8", "GATK-rescaled", "sanger"),
                             min=c(59,  64,   66, 33, 33, 33), 
                             max=c(104, 104, 105, 74, 83, 73), stringsAsFactors=FALSE)
    
    ## check compatible qualities
    rquals <- range(quals)
    unlist(lapply(1:nrow(knownquals), function(i) {
      if (rquals[1]>=knownquals$min[i] && rquals[2]<=knownquals$max[i]) knownquals$qual[i]
      else NULL
    }))
  }, silent=TRUE)
  
  if (class(cquals)=="try-error") stop("io.R/detectQualityInFASTQFile: cannot read FASTQ reads in filename=", filename, cquals)
  else cquals
}

##' Symlink-safe file/directory delete function
##'
##' Unlike unlink(), safeUnlink() does not follow symlink directories for deletion.
##' 
##' @title safeUnlink
##' @param path A character string indicating which file/directory to delete.
##' @return Nothing
##' @author Gregoire Pau
##' @export
##' @keywords internal
safeUnlink <- function(path) {
  cmd <- paste("rm -rf", path, "2>&1")
  z <- suppressWarnings(system(cmd, intern=TRUE))
  if (length(z)>0) stop("io.R/safeUnlink: cannot delete file=", path, "; ", z)
}

##' Get a filename given a directory and the object name
##'
##' @title Get a filename given a directory and the object name
##' @param dir_path A character string containing a dir path
##' @param object_name A character string containing the regular expression matching a filename in dir_path
##' @return A character vector containing an existing filename, stops if 0 or more than 1
##' @author Gregoire Pau
##' @export
##' @keywords internal
getObjectFilename <- function(dir_path, object_name) {
  files <- dir(dir_path, full.names=TRUE)
  filename <- grep(paste(object_name, sep=""), files, value=TRUE)
  if (length(filename)!=1) {
    stop("io.R/getObjectFilename: found non-unique (0 or multiple) filename matching '", object_name,
         "' in directory ", dir_path)
  }
  return(filename)
}

##' Get the filename of the variant file
##'
##' Depending on the variant caller used and the version of
##' VariantAnotation used to create the file a file might have the
##' ending vcf.gz, vcf.bgz. To function hides all this mess.
##' @title Get a vcf filename given a HTSeqGenie directory
##' @param dir_path A character string containing a dir path
##' @return A character vector containing an existing filename, stops if 0 or more than 1
##' @author Jens Reeder
##' @export
##' @keywords internal
findVariantFile <- function(save_dir){
  vcf <- getObjectFilename(file.path(save_dir, "results"),
                           "variants.vcf.b?gz$")
  return(vcf)
}

## this path is ciontructed over and over again, so we factor it out here
## similar functions for other files will follow
getAnalyzedBamFile <- function(){
   file <- file.path(getConfig('save_dir'), 'bams',
                     paste(getConfig('prepend_str'),
                           'analyzed.bam', sep="."))
   return(file)
 }


##' Overloaded yield(...) method catching truncated exceptions for FastqStreamer
##'
##' @title Overloaded yield(...) method catching truncated exceptions for FastqStreamer
##' @param fqs An instance from the FastqSampler or FastqStreamer class.
##' @return Same as FastqStreamer::yield
##' @author Gregoire Pau
##' @export
##' @keywords internal
##' @importMethodsFrom ShortRead yield
safe.yield <- function(fqs) {
  tryCatch(yield(fqs), IncompleteFinalRecord=function(x) stop("io.R/safe.yield: input gzipped file is corrupted/truncated. Aborting.\n"))
}

## file.rename only warns, but we like it to fail instead
safe.file.rename <- function(from, to){
  success <- file.rename(from, to)
  if(!all(success)){
    stop(paste("Renaming file", from, "to", to, "failed"))
  }
}

##' Read single/paired end BAM files with requested columns from the BAM
##'
##' @title Read single/paired End Bam Files
##' @param filename Path to a bam file
##' @param paired_ends A logical indicating whether the reads are paired
##' @param remove.strandness A logical indicating whether read strands should be set to "*".
##' @return GRangesList
##' @author Cory Barr
##' @export
##' @keywords internal
##' @importMethodsFrom Rsamtools ScanBamParam 
##' @importFrom Rsamtools scanBamFlag
##' @importFrom GenomicAlignments readGAlignments
readRNASeqEnds <- function(filename, paired_ends, remove.strandness=TRUE) {
  ## bam file
  scf <- scanBamFlag(isNotPassingQualityControls=FALSE, isDuplicate=FALSE)
  sbp <- ScanBamParam(flag=scf, what="groupid")
  if (paired_ends) bamfile <- BamFile(filename, asMates=TRUE)
  else bamfile <- BamFile(filename)

  ## read bam
  reads <- readGAlignmentsFromBam(bamfile, param=sbp, use.names=FALSE, with.which_label=FALSE)
  if (remove.strandness) strand(reads) <- "*"
  groupid <- values(reads)$groupid
  reads <- grglist(reads) ## use cigar information to split gapped reads into granges
  
  ## pair reads
  if (paired_ends) {
    z <- rep(groupid, elementLengths(reads))
    reads <- unlist(reads)
    reads <- split(reads, z)
  }

  ## return reads
  return(reads)
}
