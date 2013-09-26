## API based on the "logging" package from Mario Frasca
loglevels <- setNames(c(10, 20, 30, 40), c("DEBUG", "INFO", "WARN", "ERROR"))

##' Log debug (with a try statement)
##'
##' @title Log debug using the logging package
##' @param ... Arguments passed to logging::logdebug
##' @return Nothing
##' @author Gregoire Pau
##' @keywords internal
##' @export
logdebug <- function(msg) try(loglevel(msg, level="DEBUG"), silent=TRUE)

##' Log info (with a try statement)
##'
##' @title Log info using the logging package
##' @param ... Arguments passed to logging::loginfo
##' @return Nothing
##' @author Gregoire Pau
##' @keywords internal
##' @export
loginfo <- function(msg) try(loglevel(msg, level="INFO"), silent=TRUE)

##' Log warning (with a try statement)
##'
##' @title Log warning using the logging package
##' @param ... Arguments passed to logging::logwarn
##' @return Nothing
##' @author Gregoire Pau
##' @keywords internal
##' @export
logwarn <- function(msg) try(loglevel(msg, level="WARN"), silent=TRUE)

##' Log error (with a try statement)
##'
##' @title Log info using the logging package
##' @param ... Arguments passed to logging::loginfo
##' @return Nothing
##' @author Gregoire Pau
##' @keywords internal
##' @export
logerror <- function(msg) try(loglevel(msg, level="ERROR"), silent=TRUE)

loglevel <- function(msg, level) {
  timestamp <- sprintf("%s", Sys.time())
  msg <- paste(timestamp, " ", level, "::", msg, sep="")
  llevel <- loglevels[level]
  glevel <- loglevels[logging.loglevel]
  z <- llevel>=glevel
  if (!is.na(z) && z) {
    sapply(logging.handlers, function(handler) {
      handler(msg)
    })
  }
  invisible(TRUE)
}

writeToConsole <- function (msg) {
  cat(paste(msg, "\n", sep = ""))
}

writeToFile <- function (msg) {
  cat(paste(msg, "\n", sep = ""), file=logging.file, append=TRUE)
}

logReset <- function() {
  logging.handlers <<- list()
  logging.loglevel <<- "DEBUG"
  logging.file <<- NULL
}

setLevel <- function(level) {
  logging.loglevel <<- level
}

addHandler <- function(handler, file, level) {
  logging.handlers <<- c(logging.handlers, handler)
  if (!missing(file)) logging.file <<-file
}
