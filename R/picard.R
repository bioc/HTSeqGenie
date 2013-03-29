##' Generic function to call all picard command line java tools
##'
##' @title picard
##' @param tool Name of the Picard Tool, e.g. MarkDuplicates
##' @param ... Arguments forwarded to the picard tool
##' @param path full path to the picard tool jar file.
##' @return Nothing
##' @author Jens Reeder, Michael Lawrence
##' @keywords internal
picard <- function(tool, ..., path = getOption("picard.path"))
{
  args <- list(...)
  unnamed.args <- nchar(names(args)) == 0L
  quoted.argnames <- sapply(as.list(substitute(list(...)))[-1], deparse)
  names(args)[unnamed.args] <- quoted.argnames[unnamed.args]
  jar <- paste0(tool, ".jar")
  if (!is.null(path))
    jar <- file.path(path, jar)

  ## Force some settings that make sense on our big multicore machines
  full_args <- paste( "-XX:ParallelGCThreads=4 -Xmx4g -jar", jar,
                     paste(names(args), args, sep = "=", collapse = " "))
  
  retcode <- system2("java", full_args, stderr=FALSE)

  if(retcode == 127){
    stop(paste("Could not call picard tool. Make sure to specify the proper path.",
               "Command:", paste("java", full_args, sep=" "), sep="\n"))
  }
  if(retcode != 0){
    stop(paste("Picard command [", paste("java", full_args, sep=" "), "]failed for an unknown reason"))
  }
}

##' Mark duplicates in bam
##'
##' Use MarkDuplicates from PicardTools to
##' mark duplicate alignments in bam file.
##' @title markDuplicates
##' @param bamfile Name of input bam file
##' @param outfile Name of output bam file
##' @param path Full path to MarkDuplicates jar 
##' @return Path to output bam file
##' @author Jens Reeder
##' @export
markDuplicates <- function(bamfile, outfile=NULL, path=NULL) {
  if (is.null(outfile)){
    ## put it right next to input file
    outfile <- paste(gsub("\\.bam$", "", bamfile) , "marked_dups", "bam", sep='.')
  }
  if (!file.exists(dirname(outfile))){
    dir.create(dirname(outfile), recursive = TRUE)
  }
  metrics_file <- paste( gsub("\\.bam$", "", outfile), "picard_metrics.txt", sep='.')
  
  picard("MarkDuplicates", INPUT=bamfile, OUTPUT=outfile,
         METRICS_FILE=metrics_file,
         ASSUME_SORTED = TRUE,
         PROGRAM_RECORD_ID='null', ## so we see the gmap annotation
         path=path
         )
  invisible(outfile)
}

##' Check for a jar file from picard tools
##'
##' Call a tool from picard and see if it responds.
##' @title checkPicardJar
##' @param toolname Name of the Picard Tool, e.g. MarkDuplicates
##' @param path Path to folder containing picard jars
##' @return TRUE if tool can be called, FALSE otherwise
##' @author Jens Reeder
##' @keywords internal
##' @export
checkPicardJar <- function(toolname, path=getOption("picard.path")){

  if(is.null(path)) return(FALSE)
  jar.file <- file.path(path, paste(toolname,"jar", sep="."))
  if (!file.exists(jar.file))
    return(FALSE)

  ## java call comes back with return code 1, which confuses system2
  retval <- system2("java", paste("-jar", jar.file, "--version"),
                    stderr=FALSE, stdout=FALSE)
  if (retval==127){ #127 signals could not execute command at all
    return(FALSE)
  } else {
    return(TRUE)
  }
}
     
##' Mark duplicates in pipeline context
##'
##' High level function call to mark duplicates
##' in the analyzed.ba file of a pipelin run.
##' @title markDups
##' @return Nothing
##' @author Jens Reeder
##' @export
markDups <- function(){
  loginfo("picard.R/markDups: starting to mark duplicates")    

  save_dir <- getConfig("save_dir")
  picard.path <- getConfig('path.picard_tools')
  analyzed.bam <- getBams(save_dir)$"analyzed"

  safeExecute({
    markDups.bam <- markDuplicates(analyzed.bam, path=picard.path)
    file.rename(markDups.bam, analyzed.bam)
    indexBam(analyzed.bam)
  }, memtracer=getConfig.logical("debug.tracemem"))


  loginfo("picard.R/markDups: finished marking duplicates")
}
