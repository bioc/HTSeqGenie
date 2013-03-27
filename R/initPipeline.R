##' Init pipeline environment
##'
##' @title Init pipeline environment
##' @param config_filename Name of config file
##' @param config_update List of name value pairs that will update the config parameters
##' @return Nothing
##' @author Jens Reeder
##' @keywords internal
##' @export
initPipelineFromConfig <- function(config_filename, config_update) {
  loadConfig(config_filename)
  if (!missing(config_update)) updateConfig(config_update)
  checkConfig()
  Sys.umask("0022") ## file permissions are now -rw-r--r-- and dir permissions are drwxr-xr-x
  initDirs()
  initLogger()
}

##' Init Pipeline environment from previous run
##'
##' Loads the config file from a previous run
##' stored in [save_dir]/logs/config.txt
##' @title Init Pipeline environment from previous run
##' @param save_dir Save dir of a previous pipeline run
##' @param config_update List of name value pairs that will update the config parameters
##' @return Nothing
##' @author Gregoire Pau
##' @keywords internal
##' @export
initPipelineFromSaveDir <- function(save_dir, config_update) {
  config_filename <- file.path(save_dir, "logs", "config.txt")
  loadConfig(config_filename)
  updateConfig(list(save_dir=save_dir))
  if (!missing(config_update)) updateConfig(config_update)
  Sys.umask("0022")  ## file permissions are now -rw-r--r-- and dir permissions are drwxr-xr-x
  initLogger()
}

##' Runs the preprocesing steps of the pipeline
##'
##' 
##' @title Run the preprocesing steps of the pipeline
##' @param config_filename Path to configuration file
##' @param config_update List of name value pairs that will update the config parameters
##' @return Nothing
##' @author Jens Reeder
##' @keywords internal
##' @export
runPreprocessReads <- function(config_filename, config_update) {
  initPipelineFromConfig(config_filename, config_update)
  preprocessReads()
}

##' Runs the read alignment step of the pipeline
##'
##' 
##' @title Runs the read alignment step of the pipeline
##' @param config_filename Path to configuration file
##' @param config_update List of name value pairs that will update the config parameters
##' @return Nothing 
##' @author Jens Reeder
##' @keywords internal
##' @export
runAlignment <- function(config_filename, config_update) {
  initPipelineFromConfig(config_filename, config_update)
  preprocessReads()
  alignReads()
}

##' Set up NGS output dir (using save_dir from getConfig)
##'
##' @title Set up NGS output dir 
##' @return Nothing
##' @author Gregoire Pau
##' @keywords internal
##' @export
initDirs <- function(){
  ## setup master output directory
  save_dir <- getConfig("save_dir")
  setUpDirs(save_dir, overwrite=getConfig("overwrite_save_dir"))

  ## write audit and config files
  writeAudit(file.path(save_dir, "logs", "audit.txt"))
  writeConfig(file.path(save_dir, "logs", "config.txt"))

  ## Return the config explicitely, as this one should be used for all
  ## subsequenct loadConfig() calls
  return(file.path(save_dir, "logs", "config.txt"))
}

##' Init loggers (output dir log, using save_dir from getConfig, and console log)
##'
##' @title Init loggers
##' @return Nothing
##' @author Gregoire Pau
##' @keywords internal
##' @export
initLogger <- function(){ 
  ## set up logger
  debug_level <- getConfig("debug.level")
  initLog(getConfig("save_dir"), debug_level)
  addHandler(writeToConsole, level=debug_level)
}

##' Initialize the logger
##'
##' Setup logging file in {save_dir}/progress.log
##' and log sessionInfo and configuration 
##' 
##' @param save_dir Save dir of a pipeline run
##' @param debug_level One of INFO, WARN, ERROR, FATAL
##' @return Log file name
##' @keywords internal
##' @export
initLog <- function(save_dir, debug_level="INFO") {
  ## set log handlers
  logReset()
  logfilename <- file.path(save_dir, "logs", "progress.log")
  addHandler(writeToFile, file=logfilename, level=debug_level)
  setLevel(debug_level)

  ## hi!
  loginfo("runPipeline.R/initLog: Hi!")
  logfilename
}

##' Write Session information
##'
##' 
##' @title Write Session information
##' @param filename  Optional name of file. If missing, prints session information on the standard output.
##' @return Nothing
##' @author Gregoire Pau
##' @keywords internal
##' @export
writeAudit <- function(filename) {
  ## open file
  if (!missing(filename)) faudit <- file(filename, open="wt")
  else faudit <- ""

  ## horizontal bar
  hbar <- paste(c(rep("-", 140), "\n"), collapse="")

  ## system
  cat(hbar, "system info:\n", file=faudit, sep="")
  cat("date=", date(), "\n", file=faudit, sep="")
  cat("hostname=", Sys.info()['nodename'], "\n", file=faudit, sep="")
  cat("pipeline.version=", sessionInfo()$otherPkgs$HTSeqGenie$Version, file=faudit, "\n")
  cat("\n", file=faudit)

  ## session info
  cat(hbar, "session info:\n", file=faudit, sep="")
  ssi <- paste(capture.output(sessionInfo()), collapse="\n")
  cat(ssi, "\n", file=faudit, sep="")
  cat("\n", file=faudit)

  ## gsnap
  cat(hbar, "gsnap --version:\n", file=faudit, sep="")
  if (system("gsnap", ignore.stderr=TRUE)!=127) {
    cat(paste(system("gsnap --version  2>&1", intern=TRUE), collapse="\n"), "\n", file=faudit)
  } else {
    cat("gsnap is not present in the PATH\n", file=faudit)
  }
  cat("\n", file=faudit)

  ## samtools
  cat(hbar, "samtools:\n", file=faudit, sep="")
  if (system("samtools", ignore.stderr=TRUE)!=127) {
    out <- suppressWarnings(system("samtools 2>&1", intern=TRUE))
    cat(paste(out, collapse="\n"), "\n", file=faudit)
  } else {
    cat("samtools is not present in the PATH\n", file=faudit)
  }
  cat("\n", file=faudit)
  
  ## PATH
  cat(hbar, "PATH:\n", file=faudit, sep="")
  path <- Sys.getenv("PATH")
  path <- gsub(":", ":\n", path)
  cat(path, "\n", file=faudit)
  cat("\n", file=faudit)
  
  ## close audit
  if (!missing(filename)) close(faudit)
}
