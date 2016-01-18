##' @importMethodsFrom Rsamtools asBam 
wrapGsnap <- function(fp1, fp2, gsnapParams, outputdir, prepend_str, remove.nomapping=FALSE) {
  ## get config parameters
  nbthreads_perchunk <- getConfig.integer("alignReads.nbthreads_perchunk")
  use_gmapR_gsnap <- getConfig.logical("alignReads.use_gmapR_gsnap")
                               
  ## build gsnap command
  args <- buildGsnapArgs(fp1, fp2, gsnapParams, output=file.path(outputdir, prepend_str))

  ## empty input files? create a temporary fake FASTQ file to not fail gsnap
  sizefp1 <- file.info(fp1)$size
  if (sizefp1==0) {
    fakeread <- "@fakeread\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
    cat(fakeread, file=fp1)
    if (!is.null(fp2)) cat(fakeread, file=fp2)
  }

  ## call gsnap command
  if (use_gmapR_gsnap) {
    loginfo(paste("alignReads.R/alignReadsChunk: calling gsnap (from gmapR) with args=", args))
    ret <- try(suppressWarnings(gmapR:::..gsnap(args)), silent=TRUE)
    if (class(ret)=="try-error" || !any(grepl("^Processed \\d+ queries", ret, perl=TRUE))) {
      stop("alignReads.R/alignReadsChunk: gsnap (from gmapR) failed with args='", args, "'")
    }
  }
  else {
    loginfo(paste("alignReads.R/alignReadsChunk: calling gsnap (from PATH) with args=", args))
    ret <- suppressWarnings(system(paste("gsnap", args), intern=TRUE))
    if (!any(grepl("^Processed \\d+ queries", ret, perl=TRUE))) {
      stop("alignReads.R/alignReadsChunk: gsnap (from PATH) failed with args='", args, "'")
    }
  }

  ## output sams
  ## gsnap sam files have no ending, so are hard to recognize.
  ## We find them by excluding files that do not look like sam files. This is fragile.
  outsams <- grep(paste("/", prepend_str, sep=""), dir(outputdir, full.names=TRUE), value=TRUE)
  outsams <- grep("ba[mi]$", outsams, invert=TRUE, value=TRUE)
  outsams <- grep("fastq$", outsams, invert=TRUE, value=TRUE)

  ## remove nomapping?
  if (remove.nomapping) {
    z <- grep("nomapping", outsams)
    if (length(z)>0) {
      unlink(outsams[z])
      outsams <- outsams[-z]
    }
  }

  ## empty input files? restore empty files and remove the fakeread
  if (sizefp1==0) {
    cat("", file=fp1)
    if (!is.null(fp2)) cat("", file=fp2)
    
    ## remove the fake read in outsams
    for (sam in outsams) {
      a <- readLines(sam)
      z <- grep("^fakeread", a)
      if (length(z)>0) writeLines(a[-z], sam)
    }
  }

  ## convert sams to bams, using multiple threads
  bams <- as.character(mclapply(outsams, mc.cores=nbthreads_perchunk,
                                function(sam) {
                                  bam <- asBam(sam, sam)
                                  unlink(sam)
                                  return(bam)
                                }))
  if(length(bams) != length(outsams) ||
     !all(file.exists(bams))){
    stop("Error during sam to bam conversions. Some bam files are missing.")
  }

  invisible(bams)
}

buildGsnapArgs <- function(fp1, fp2=NULL, gsnapParams, output) {
  ## core args
  args <- paste(gsnapParams, "--split-output", output, fp1)
  if (!is.null(fp2)) args <- paste(args, "-a paired", fp2)

  ## redirect stderr to stdout
  args <- paste(args, "2>&1")
  args
}
