test.sclapply <- function() {
  ## test 1: single job
  listIterator.init(1)
  z <- unlist(sclapply(listIterator.next, function(x, ...) x^2, max.parallel.jobs=3))
  checkEquals(z, 1)

  ## test 2: standard
  listIterator.init(1:100)
  z <- unlist(sclapply(listIterator.next, function(x, ...) x^2, max.parallel.jobs=3))
  checkEquals(z, (1:100)^2)

  ## test 3: one node
  listIterator.init(1:100)
  z <- unlist(sclapply(listIterator.next, function(x, ...) x^2, max.parallel.jobs=1))
  checkEquals(z, (1:100)^2)

  ## test 4: different running times
  listIterator.init(1:23)
  z <- unlist(sclapply(listIterator.next, function(x, ...) {Sys.sleep(runif(1)*0.5) ; x^2}, max.parallel.jobs=5))
  checkEquals(z, (1:23)^2)

  ## test 5: exception
  listIterator.init(1:23)
  f <- function(x, ...) { if (x==4) stop("illegal operation") ; x^2}
  z <- try(sclapply(listIterator.next, f, max.parallel.jobs=3), silent=TRUE)
  checkEquals(class(z), "try-error")

  ## test 6: stop.onfail=FALSE
  listIterator.init(1:23)
  f <- function(x, ...) { if (x==4) stop("illegal operation") ; x^2}
  z <- sclapply(listIterator.next, f, max.parallel.jobs=3, stop.onfail=FALSE)
  checkEquals(class(z[[4]]), "try-error")
  checkEquals(unlist(z[-4]), ((1:23)^2)[-4])

  ## test 7: tracer timer
  tracefun <- function(type, ...) {}
  listIterator.init(1:5)
  z <- unlist(sclapply(listIterator.next, function(x, ...) {Sys.sleep(runif(1)*2) ; x^2}, max.parallel.jobs=3,
                       tracefun=tracefun, tracefun.period=1))
  checkEquals(z, (1:5)^2)
  
  ## test 8: tracer exception
  listIterator.init(1:23)
  f <- function(x, ...) { if (x==4) stop("illegal operation") ; x^2}
  z <- sclapply(listIterator.next, f, max.parallel.jobs=3, stop.onfail=FALSE,
                                 tracefun=tracefun, tracefun.period=1)
  checkEquals(class(z[[4]]), "try-error")
  checkEquals(unlist(z[-4]), ((1:23)^2)[-4])
}

test.tryKeepTraceback <- function() {
  ## simple pass
  z <- tryKeepTraceback("ok")
  checkEquals(z, "ok")

  ## simple fail
  z <- tryKeepTraceback(stop("ko"))
  checkEquals(class(z), "try-error", "tryKeepTraceback() did not fail")
  checkEquals(grep("ko", getTraceback(z)), 1, "getTraceback() did not return the message 'ko'")
}
