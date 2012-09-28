##' Load configuration file
##'
##' Loads the indicated configuration file. Creates and installs
##' a global variable that should be accesssed only via getConfig().
##' 
##' @param filename Path to configuration file
##' @return Nothing. Called for its side effect, which is setting
##'         the global config variable.
##' @keywords internal
##' @export
loadConfig <- function(filename) {
  ## load configuration file
  if (missing(filename)){
    config <- list(template_config="default-config.txt")
  }
  else{
    config <- parseDCF(filename)
  }
  assign("config", config, 1)
}

##' Update the existing config
##'
##' @param tconfig List of configuration name value pairs
##' @return Nothing.
##' @keywords internal
##' @export
updateConfig <- function(tconfig) {
  if (!exists("config")) stop("'config' has not been loaded")
  if (length(tconfig)>0) {
    config[names(tconfig)] <- as.character(tconfig)
    assign("config", config, 1)
  }
}

##' Get a configuration parameter
##'
##' @param p Name of parameter
##' @param stop.ifempty throw error if value is not set, otherwise returns NULL
##' @return If parameter is missing, return the config list
##' otherwise return the value of the parameter name as a character string
##' throws an exception if the parameter is not present in the config
##' @keywords internal
##' @export
getConfig <- function(p, stop.ifempty=FALSE) {
  if (!exists("config")) {
    config <- NULL ## to avoid the 'no visible binding' R CMD check warning
    stop("'config' has not been loaded")
  }
  
  if (missing(p)) return(config)
  else {
    if (isConfig(p)) {
      a <- config[[p]]
      if (a=="") a <- NULL
    }
    else stop("config parameter '", p, "' must be present")
  }
  
  if (stop.ifempty && is.null(a)) stop("config parameter '", p, "' must be specified")
  a
}

##' Test the presence of the parameter in the current config
##' @param parameter Name of parameter
##' @return TRUE if present, FALSE otherwise
##' @keywords internal
##' @export
isConfig <- function(parameter) {
  if (!exists("config")) {
    config <- NULL ## to avoid the 'no visible binding' R CMD check warning
    stop("'config' has not been loaded")
  }
  
  z <- config[[parameter]]
  !is.null(z)
}

##' Check if a config parameter is a numeric
##'
##' Throws exception if value can't be cast into numeric
##' 
##' @param p Name of parameter
##' @param ... Extra params passed to getConfig
##' @return Value of parameter as numeric
##' @keywords internal
##' @export
getConfig.numeric <- function(p, ...) {
  a <- getConfig(p, ...)
  if (!is.null(a)) {
    a <- suppressWarnings(as.numeric(a))
    if (is.na(a)) stop("config parameter '", p, "' must be a numeric")
  }
  a
}

##' Check if a config parameter is an integer
##'
##' Throws exception if value is no integer
##' 
##' @param p Name of parameter
##' @param tol Tolerance that controls how far a value can be from the next integer.
##' @param ... Additinal parameters passed to getConfig() 
##' @return Value of parameter as integer
##' @keywords internal
##' @export
getConfig.integer <- function(p, tol=1e-8, ...) {
  a <- getConfig.numeric(p, ...)
  if (!is.null(a)) {
    if (abs(as.integer(a)-a)>tol) stop("config parameter '", p, "' must be an integer")
    a <- as.integer(a)
  }
  a
}

##' Check if a config parameter has a logical value
##'
##' Throws exception if value is not logical
##' 
##' @param p Name of parameter
##' @param ... extra params passed to getConfig
##' @return Logical value of parameter
##' @keywords internal
##' @export
getConfig.logical <- function(p, ...) {
  a <- getConfig(p, ...)
  if (!is.null(a)) {
    a <- as.logical(a)
    if (is.na(a)) stop("config parameter '", p, "' must be a boolean value")
  }
  a
}

##' Read and parse a configuration file
##' 
##' From a file like
##' x1: y1
##' x2: y2
##' extract field, using the rules:
##'   - split on ':'
##'   - first element of split id name of parameter, second is value
##'   - trailing whitespaces (tabs and spaces) are removed
##'   - comments (text flow starting with #) are removed
##' @param filename File name
##' @return Named list
##' @keywords internal
##' @export
parseDCF <- function(filename) {
  if (!file.exists(filename)) stop('config.R/parseDCF: cannot open file "',filename,'"',sep='')

  ## read file
  tt <- paste(readLines(filename, warn=FALSE), '')

  ## remove comments
  tt <- gsub("#.*$", "", tt)
  
  ## split first :
  z <- regexpr(':',tt)
  z[z<0] <- nchar(tt)[z<0]+1
  aa <- matrix('', nrow=length(tt), ncol=2)
  aa[,1] <- substr(tt, 1, z-1)
  aa[,2] <- substr(tt, z+1, nchar(tt))
  tt <- apply(aa, 1, as.list)

  ## remove trailing whitespaces
  tt[lapply(tt, length)<2] <- NULL
  tt <- lapply(tt, lapply, function (z) {gsub("^[\t ]*|[\t ]*$", "", z)})

  ## rebuild list
  names(tt) <- sapply(tt,'[[', 1)
  tt <- sapply(tt, '[[', 2)
  tt <- as.list(tt)
  
  ## remove unnamed parameters
  tt <- tt[names(tt)!=""]
 
  tt
}

##' Write a config file
##'
##' Writes the currently active configuration to file
##'
##' @param config.filename Optional name of output file. If missing, print the config file on the standard output.
##' @return Name of saved file
##' @keywords internal
##' @export
writeConfig <- function(config.filename) { 
  ## open file
  if (!missing(config.filename)) fconfig <- file(config.filename, open="wt")
  else fconfig <- ""
    
  ## current configuration file
  conf <- getConfig()
  conf <- paste(names(conf),": ", as.character(conf), sep="", collapse="\n")
  cat(conf, file=fconfig, sep="")
  cat("\n", file=fconfig)

  ## close config
  if (!missing(config.filename)) {
    close(fconfig)
    return(config.filename)
  }
  else invisible(NULL)
}

##' Build a configuration file based on a list of parameters
##'
##' @title Build a configuration file based on a list of parameters
##' @param config_filename The path of a configuration filename.
##' @param ... A list of named value pairs.
##' @return Nothing.
##' @seealso runPipeline
##' @author Gregoire Pau
##' @keywords internal
##' @export
buildConfig <- function(config_filename, ...) {
  fconfig <- file(config_filename, open="wt")
  config <- list(...)
  config <- paste(names(config),": ", as.character(config), sep="", collapse="\n")
  cat(config, file=fconfig, sep="")
  cat("\n", file=fconfig)
  close(fconfig)
}
