.onLoad <- function(libname, pkgname) {
  ## disable ShortRead OpenMP (which crashes FastqStreamer when invoked inside a mcparallel() statement)
  .Call(ShortRead:::.set_omp_threads, 1L)
}
.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.22")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}
