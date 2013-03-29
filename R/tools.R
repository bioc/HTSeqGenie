##' Reload package source code
##'
##' When developing code this function can be used to quickly
##' reload all of the packages code, without installing it.
##' @param dirname  Directore with files to source
##' @return Nothing
##' @export
##' @keywords internal
resource <- function(dirname=".") {
  dirname <- file.path(dirname, "R")
  rfiles <- dir(dirname, pattern="*.R$", full.names=TRUE)
  z = sapply(rfiles, source)
  if (length(rfiles)==0) cat("did not source anything in dirname=", dirname, "\n")
  else cat("sourced:", paste(rfiles, sep=", "), "\n")
}

##' Get a package file
##'
##' Magically get package files from the inst directory,
##' which will be in different location, depending on whether
##' we run in:
##'  - local mode: if interactive() is TRUE
##'  - package mode: if interactive() is FALSE
##' @param filename Name of package file
##' @param package Name of the package the file is coming from
##' @param mustWork Boolean, will stop the code if set tot TRUE and file not found
##' otherwise returns Nothing.
##' @return relative path to requested file
##' @export
##' @keywords internal
getPackageFile <- function(filename, package="HTSeqGenie", mustWork=TRUE) {
  fp <- file.path("inst", filename)
  if (interactive() && file.exists(fp)) return(fp)
  else system.file(filename, package=package, mustWork=mustWork)
}

##' Returns memory usage in bytes
##'
##' For debugging.
##' @return Memory usage in bytes
##' @export
##' @keywords internal
getMemoryUsage <- function() {
  mem <- try(system("free", intern=TRUE, ignore.stderr=TRUE), silent=TRUE)
  if (class(mem)=="try-error") res <- c(NA, NA, NA, NA, NA, NA)
  else {
    ## parse "free" command output
    mem <- strsplit(mem[2], " ")[[1]]
    mem <- mem[mem!=""]
    res <- as.numeric(mem[2:7])*1024
  }
  names(res) <- c("total", "used", "free", "shared", "buffers", "cached")
  res
}

##' Process chunk in the pipeline framework
##' 
##' High-level pipeline-specific version of sclapply, with chunk loggers and safeExecute
##' @title Process chunk in the pipeline framework
##' @param inext A function (without argument) returning an object to process; NULL if none left; this function is run in the main thread
##' @param fun Function to process the object returned by inext; this function is run in children threadfunction to apply to a chunk
##' @param nb.parallel.jobs number of parallel jobs
##' @return Nothing
##' @author Gregoire Pau
##' @export
##' @keywords internal
processChunks <- function(inext, fun, nb.parallel.jobs) {
  logdebug("tools.R/processChunks: starting...")
  
  ## get config parameters
  if (missing(nb.parallel.jobs)) nb.parallel.jobs <- getConfig.integer("num_cores")

  chunk.dir <- getConfig('chunk_dir')
  ## tracer
  tracefun.starttime <- list()
  tracefun <- function(type, ...) {
    with(list(...), {
      if (type=="timer") {
        logdebug(paste("tools.R/processChunks: waiting for chunkid=[",
                       paste(which(sjobs=="running"), collapse=", "), "] ..."))
      }
      else {
        logfile <- file.path(chunk.dir, sprintf("chunk_%06d", chunkid),
                             "logs", "progress.log")
        if (type=="start") {
          logdebug(paste("tools.R/processChunks: starting chunkid=", chunkid,
                         "; see logfile=", logfile))
          tracefun.starttime[[chunkid]] <<- proc.time()["elapsed"]
          Sys.sleep(2)
        }
        if (type=="done") {
          endtime <- proc.time()["elapsed"]
          elapsed.time <- (endtime - tracefun.starttime[[chunkid]])/60
          logdebug(paste("tools.R/processChunks: done with chunkid=", chunkid, "; elapsed.time=",
                         round(elapsed.time, 3), "minutes"))
        }
        if (type=="error") {
          logerror(paste("tools.R/processChunks: error while processing chunkid=", chunkid,
                         "; see logfile=", logfile))
        }
      }
    })
  }
  
  ## decorate fun to set up logger and catch exception with traceback
  funlog <- function(..., chunkid) {
    save_dir <- file.path(chunk.dir, sprintf("chunk_%06d", chunkid))
    setUpDirs(save_dir, overwrite="overwrite")
    writeConfig(file.path(save_dir, "logs", "config.txt"))
    initLog(save_dir, getConfig("debug.level"))
    safeExecute(fun(..., save_dir=save_dir), memtracer=FALSE)
  }

  ## clean up memory before firing new threads, this should decrease a lot of memory usage
  gc(reset=TRUE)
  gc()
  
  ## process chunks
  res <- sclapply(inext=inext, fun=funlog, max.parallel.jobs=nb.parallel.jobs, stop.onfail=TRUE,
           tracefun=tracefun, tracefun.period=60)

  ## all jobs failed? => stop
  if (all(sapply(res, class)=="try-error")) {
    stop("tools.R/processChunks: all jobs failed")
  }
  
  logdebug("tools.R/processChunks: done")
}

##' Scheduled parallel processing
##'
##' @param inext A function (without argument) returning an object to process; NULL if none left; this function is run in the main thread
##' @param fun Function to process the object returned by inext; this function is run in children thread
##' @param max.parallel.jobs Number of jobs to start in parallel 
##' @param ... Further arguments passed to fun
##' @param stop.onfail Throw error if one 
##' @param tracefun Callback function that will be executed in a separate thread 
##' @param tracefun.period Time intervall between calls to tracefun
##' @return Return value of applied function
##' @export
##' @keywords internal
sclapply <- function(inext, fun, max.parallel.jobs, ..., stop.onfail=TRUE, tracefun=NULL, tracefun.period=60) {
  ## initialise scheduler
  sjobs <- character(0) ## jobs (state)
  rjobs <- list() ## jobs (result)
  pnodes <- vector(mode="list", max.parallel.jobs) ## nodes (process)
  jnodes <- rep(NA, max.parallel.jobs) ## nodes (job id)

  ## cleanup procedure, based on mclapply
  cleanup <- function() {
    z <- pnodes[!sapply(pnodes, is.null)]
    if (length(z)>0) mccollect(z, wait=FALSE, timeout=4)
    z <- pnodes[!sapply(pnodes, is.null)]
    if (length(z)>0) mccollect(z, wait=FALSE, timeout=4)
    z <- pnodes[!sapply(pnodes, is.null)]
    if (length(z)>0) try(parallel:::mckill(z, tools::SIGTERM), silent=TRUE)
    z <- pnodes[!sapply(pnodes, is.null)]
    if (length(z)>0) mccollect(z, wait=FALSE, timeout=4)
    z <- pnodes[!sapply(pnodes, is.null)]
    if (length(z)>0) mccollect(z, wait=FALSE, timeout=4)
  }
  on.exit(cleanup())

  ## start scheduler
  collect.timeout <- 2 ## 2 seconds wait between each iteration
  inextdata <- NULL
  i <- 0
  repeat {
    ## is there a new job to process?
    if (is.null(inextdata)) inextdata <- inext()
    
    ## are all the jobs done?
    if (is.null(inextdata)) {
      if (length(sjobs)==0) break ## no jobs have been run
      if (all(sjobs=="done")) break
    }
    
    ## fire fun(inextdata) on node i
    repeat {
      ## periodic trace
      if (!is.null(tracefun)) {
        current <- proc.time()["elapsed"]
        if (!exists("last.tracefun") || (current-last.tracefun)>tracefun.period) {
          tracefun(type="timer", sjobs=sjobs, jnodes=jnodes)
          last.tracefun <- current
        }
      }
      
      i <- (i %% length(pnodes))+1
      process <- pnodes[[i]]
      if (!is.null(process)) {
        status <- mccollect(process, wait=FALSE, timeout=collect.timeout) ## wait collect.timeout seconds        
        if (is.null(status)) {
          ## node busy
          fire <- FALSE
        }
        else {
          ## node done: save results and fire a new job
          mccollect(process) ## kill job
          ## did an error occur?
          if (class(status[[1]])=="try-error") { 
            if (!is.null(tracefun)) {
              tracefun("error", sjobs=sjobs, jnodes=jnodes, chunkid=jnodes[i], error=status[[1]])
            }
            if (stop.onfail) {
              stop(paste("tools.R/sclapply: error in chunkid=", jnodes[i], ": ", status[[1]], sep=""))
            }
          }
          if (!is.null(tracefun)) tracefun("done", chunkid=jnodes[i])
          rjobs[jnodes[i]] <- status
          sjobs[jnodes[i]] <- "done"
          fire <- TRUE
        }
      } else {
        ## virgin node
        fire <- TRUE
      }
      
      ## fire a new job
      if (fire && !is.null(inextdata)) {
        jnodes[i] <- length(sjobs)+1
        sjobs[jnodes[i]] <- "running"
        if (!is.null(tracefun)) tracefun("start", chunkid=jnodes[i])
        pnodes[[i]] <- mcparallel(fun(inextdata, chunkid=jnodes[i], ...), mc.set.seed=TRUE)
        inextdata <- NULL
        break
      }

      ## no more job
      if (is.null(inextdata)) break
    }
  }
  
  rjobs
}

##' Show memory usage
##'
##' For debugging purposes only.
##' Show memory usage if config variable
##' @return Nothing
##' @export
##' @keywords internal
traceMem <- function() {
  wmem <- getMemoryUsage()
  wmem <- (wmem["used"] - wmem["cached"])/1e9
  msg <- sprintf("tools.R/traceMem: wired.mem=%f GiB \n", wmem)
  logdebug(msg)
}

##' Wrapper around try-catch
##'
##' @param expr Expression to evaluate
##' @return Result of expression or error if thrown
##' @export
##' @keywords internal
tryKeepTraceback <- function(expr) {
  withRestarts(withCallingHandlers(expr,
                                   error=function(e) {
                                     error <- e
                                     calls <- sys.calls()
                                     invokeRestart("myAbort", e, calls)
                                   }),
               myAbort=function(e, calls){
                 err <- list(error=e, calls=calls)
                 class(err) <- c("try-error")
                 return(err)
               })
}

##' Get traceback from tryKeepTraceback()
##'
##' @param mto An object of the try-error class
##' @return Traceback as a string
##' @export
##' @keywords internal
getTraceback <- function(mto) {
  ## get trace
  trace <- lapply(mto$calls, function(call) {
    z = format(call)[1]
    attr(z, "srcref") <- attr(call, "srcref")
    z
  })
  
  ## reverse trace
  trace <- trace[length(trace):1]

  ## skip error handler calls
  skipFunctions <- c("withRestarts", "withCallingHandlers", "invokeRestart", "withOneRestart",
                     "tryKeepTraceback", "doWithOneRestart", ".handleSimpleError", "simpleError",
                     "doTryCatch", "tryCatchOne", "tryCatchList", "tryCatch", "sendMaster")
  utrace <- unlist(trace)
  z <- unique(unlist(sapply(1:length(skipFunctions), function(i) grep(skipFunctions[i], utrace))))
  trace <- trace[-z]

  ## build trace using traceback()
  capture.output(traceback(trace))
}

##' Execute function in try catch with trace function
##'
##' Requires the logger to be set
##' @param expr Expression to safely execute 
##' @param memtracer A boolean, to enable/disable a periodic memory tracer. Default is TRUE.
##' @return Nothing
##' @author Gregoire Pau
##' @export
##' @keywords internal
safeExecute <- function(expr, memtracer=TRUE) {
  ## start memory tracer
  if (memtracer) pmemtracer <- mcparallel(repeat {Sys.sleep(60) ; traceMem()})

  ## execute expr in a child thread (to save memory from the)
  z <- mccollect(mcparallel(tryKeepTraceback(expr)))
  if (length(z)>0) z <- z[[1]]
  else z <- NULL
  
  ## kill memory tracer
  if (memtracer) parallel:::mckill(pmemtracer, SIGTERM)
  
  if (class(z)=="try-error") {
    logerror("tools.R/safeExecute: caught exception:")
    logerror(z$e)
    logerror("tools.R/safeExecute: traceback:")
    sapply(getTraceback(z), logerror)
    stop(z$e)
  } else invisible(z)
}

##' Create a iterator on a list
##'
##' Only one listIterator object can be open at any time.
##' 
##' @title Create a iterator on a list
##' @param x A list.
##' @return Nothing
##' @author Gregoire Pau
##' @export
##' @keywords internal
listIterator.init <- function(x) {
  listIterator.n <<- 0
  listIterator.x <<- as.list(x)
}

##' Get reads from the listIterator
##'
##' @title Get reads from the listIterator
##' @return An object. NULL if there are no more objects in the listIterator.
##' @author Gregoire Pau
##' @seealso listIterator.init
##' @export
##' @keywords internal
listIterator.next <- function() {
  if (!exists("listIterator.n")) {
    listIterator.n <- NULL ## to avoid the 'no visible binding' R CMD check warning
    listIterator.x <- NULL ## to avoid the 'no visible binding' R CMD check warning
    stop("tools.R/listIterator.next: listIterator.init() has not been called")
  }
  
  listIterator.n <<- listIterator.n + 1
  if (listIterator.n<=length(listIterator.x)) return(listIterator.x[[listIterator.n]])
  else return(NULL)
}

##' Get the list of chunk directories
##'
##' @title Get the list of chunk directories
##' @return List of chunk directories
##' @author Gregoire Pau
##' @export
##' @keywords internal
getChunkDirs <- function() {
  chunkdirs <- dir(getConfig('chunk_dir'), full.names=TRUE)
  sort(chunkdirs[grep("chunk_", chunkdirs)])
}

##' Set the base directory for the chunks
##' 
##' @title Set the base directory for the chunks
##' @return path to chunk dir
##' @author Jens Reeder
##' @keywords internal
setChunkDir <- function() {

  tmp.dir <- getConfig("tmp_dir")
  
  if(is.null(tmp.dir)){
    chunk.dir <- file.path(getConfig("save_dir"), "chunks")
  }
  else {
    chunk.dir <- file.path(tmp.dir, paste(c("chunks_", sample(letters,8)), collapse=""))
    makeDir(chunk.dir)
  }
  updateConfig(list(chunk_dir=chunk.dir)) 
  return(chunk.dir)
}

##' Get bam files of a pipeline run
##'
##' @title Get bam files of a pipeline run
##' @param save_dir Save directory of a pipeline run
##' @return named list of bam files
##' @author Gregoire Pau
##' @export
##' @keywords internal
getBams <- function(save_dir) {
  bfiles <- dir(file.path(save_dir, "bams"), full.names=TRUE)
  bfiles <- bfiles[grep(".bam$", bfiles)] ## keep only files ending by .bam
  btypes <- strsplit(bfiles, "\\.") ## split on dots
  btypes <- sapply(btypes, function(z) z[length(z)-1]) ## keep the next to last
  bfiles <- as.list(bfiles)
  names(bfiles) <- btypes
  bfiles
}
