test.createTmpDir <- function(){
  tmp1 <- HTSeqGenie:::createTmpDir()
  checkTrue(file.exists(tmp1), "createTmpDir() creates temporary directory")

  tmp2 <- HTSeqGenie:::createTmpDir()
  checkTrue(file.exists(tmp1))
  checkTrue(tmp1 != tmp2, "createTmpDir does return uniq dir names")
  unlink(tmp2, recursive=TRUE)
  unlink(tmp1, recursive=TRUE)

  tmp3 <- HTSeqGenie:::createTmpDir(prefix='test123_')
  checkTrue(file.exists(tmp3), "createTmpDir() works with a prefix")
  unlink(tmp3, recursive=TRUE)
}
  
test.FastQStreamer.getReads.pefq <- function() {
  ## paired-end fastq
  ## setup test framework
  config.update <- list(chunk_size=10)
  setupTestFramework(config.filename="test-data/test_config.txt", config.update=config.update,
                     testname="test.FastQStreamer.getReads.pefq")

  ## init FastQStreamer
  FastQStreamer.init(input_file=getConfig("input_file"), input_file2=getConfig("input_file2"),
                     chunk_size=getConfig.integer("chunk_size"))

  ## read 10 reads
  lreads <- FastQStreamer.getReads()
  checkEquals(length(lreads[[2]]), 10, "FastQStreamer.getReads() doesn't return 10 reads")
  
  ## release FastQStreamer
  FastQStreamer.release()
}

test.FastQStreamer.getReads.pefq.subsample <- function() {
  ## paired-end fastq with subsample
  ## setup test framework
  config.update <- list(chunk_size=100, subsample_nbreads=17)
  setupTestFramework(config.filename="test-data/test_config.txt", config.update=config.update,
                     testname="test.FastQStreamer.getReads.pefq.subsample")

  ## init FastQStreamer
  FastQStreamer.init(input_file=getConfig("input_file"), input_file2=getConfig("input_file2"),
                     chunk_size=getConfig.integer("chunk_size"),
                     subsample_nbreads=getConfig.integer("subsample_nbreads"))

  ## read 17 reads
  lreads <- FastQStreamer.getReads()
  checkEquals(length(lreads[[2]]), 17, "FastQStreamer.getReads() doesn't return 17 reads")
  
  ## release FastQStreamer
  FastQStreamer.release()
  lreads
}

test.FastQStreamer.subsampler.isdeterministic <- function() {
  ## check if the subsampler is deterministic
  lreads1 <- test.FastQStreamer.getReads.pefq.subsample()
  lreads2 <- test.FastQStreamer.getReads.pefq.subsample()
  
  ## convert to character strings
  reads1 <- as.character(sread(lreads1[[1]]))
  reads2 <- as.character(sread(lreads2[[1]]))

  ## check
  checkEquals(reads1, reads2, "FastQStreamer is not deterministic")
}

test.FastQStreamer.getReads.segz <- function() {
  ## single-end fastq.gz
  ## setup test framework
  config.update <- list(input_file="inst/test-data/reads.fastq.gz", chunk_size=10, quality_encoding="illumina1.8")
  setupTestFramework(config.filename="test-data/test_config.txt", config.update=config.update,
                     testname="test.FastQStreamer.getReads.segz")

  ## init FastQStreamer
  FastQStreamer.init(input_file=getConfig("input_file"), input_file2=getConfig("input_file2"),
                     chunk_size=getConfig.integer("chunk_size"))

  ## read 10 reads
  lreads <- FastQStreamer.getReads()
  checkEquals(length(lreads[[1]]), 10, "FastQStreamer.getReads() doesn't return 10 reads")

  ## release FastQStreamer
  FastQStreamer.release()
}

test.writeAudit <- function() {
  ## setup test framework
  config.filename <- "test-data/test_config.txt"
  setupTestFramework(config.filename=config.filename, testname="test.writeAudit")
  writeAudit("audit.txt")
}

test.safeUnlink <- function() {
  tmpdir <- tempfile()
  dir.create(tmpdir)
  currentdir <- getwd()
  setwd(tmpdir)

  ## create blob/tk and sblob linking to blob
  system("mkdir blob ; touch blob/tk ; ln -s blob sblob")

  ## safeUnlink should not delete tk!
  safeUnlink("sblob")
  checkTrue(!file.exists("sblob"), "safeUnlink did not work!")
  checkTrue(file.exists("blob/tk"), "safeUnlink did not work!")

  setwd(currentdir)
}
