.onLoad <- function(libname, pkgname) {
  ## disable ShortRead OpenMP (which crashes FastqStreamer when invoked inside a mcparallel() statement)
  .Call(ShortRead:::.set_omp_threads, 1L)
}
