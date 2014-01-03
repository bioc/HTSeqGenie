test.findTemplate <- function() {  
  checkTrue(file.exists(HTSeqGenie:::findTemplate("default-config.txt")), "findTemplate() finds template in package")

  tfile <- "template_test_config.txt"
  file(file.path(tempdir(), tfile), "w")
 
  ## make sure we can find files in absolute paths
  abs.path <- file.path(tempdir(), tfile)
  checkTrue(file.exists(HTSeqGenie:::findTemplate(abs.path)), "findTemplate() finds template in absolute path")

  ## test finding in current dir 
  cwd <- getwd()
  setwd(tempdir())
  checkTrue(file.exists(HTSeqGenie:::findTemplate(tfile)), "findTemplate() finds template in current dir")
  setwd(cwd)

  Sys.setenv(HTSEQGENIE_CONFIG=tempdir())
  checkTrue(file.exists(HTSeqGenie:::findTemplate(tfile)), "findTemplate() finds template from HTSEQGENIE_CONFIG env")

  ## create env with multiple values
  Sys.setenv(HTSEQGENIE_CONFIG=paste(c("bla",tempdir()), collapse=":"))
  checkTrue(file.exists(HTSeqGenie:::findTemplate(tfile)), "findTemplate() finds template from HTSEQGENIE_CONFIG env with multiple entries")

  ## clean up
  unlink(file.path(tempdir(), tfile))
}

test.checkConfig.analyzeVariants <- function() {

  setupTestFramework(config.filename="test-data/test_config.txt",
                     config.update=list(),
                     testname="test.checkConfig.analyzeVariants")

  checkTrue(HTSeqGenie:::checkConfig.analyzeVariants(),
            "checkConfig.analyzeVariants() passes valid config")

  updateConfig(list(analyzeVariants.method = "bogus"))
  checkException(HTSeqGenie:::checkConfig.analyzeVariants(),
            "checkConfig.analyzeVariants() detects bogus method", silent=TRUE)

  updateConfig(list(analyzeVariants.method = 'VariantTools',
                     analyzeVariants.callingFilters='nonRef, nonNRef'))
  checkTrue(HTSeqGenie:::checkConfig.analyzeVariants(),
            "checkConfig.analyzeVariants() checks callingFilter config")

  updateConfig(list(analyzeVariants.callingFilters='nonRef, XRef'))
  checkException(HTSeqGenie:::checkConfig.analyzeVariants(),
                 "checkConfig.analyzeVariants() detects filters not provided by VT",
                 silent=TRUE)
  updateConfig(list(analyzeVariants.callingFilters=''))
  checkTrue(HTSeqGenie:::checkConfig.analyzeVariants(),
                 "checkConfig.analyzeVariants() handles empty vector configs")
}   
