test.parseDCF <- function() {
  ## write test config file
  tt <- list(
    "parameter1: 12.5", 
    " parameter2: 2 # comment1",
    " ",
    "",
    "# # comment2",
    "  parameter3:yuk  ",
    "parameter4 \t:\tblob ",
    "parameter5 :"       
    )
  configfile <- tempfile()
  cat(paste(tt, collapse="\n"), file=configfile)

  ## parse and check config file
  config = parseDCF(configfile)
  checkEquals(length(config), 5) ## 5 parameters
  checkEquals(config$parameter1, "12.5") ## no trailing whitespaces in "12.5"
  checkEquals(as.numeric(config$parameter1), 12.5)
  checkEquals(as.numeric(config$parameter2), 2) ## no trailing whitespaces in "parameter2"
  checkEquals(config$parameter3, "yuk") ## no trailing whitespaces in "yuk"
  checkEquals(config$parameter4, "blob") ## no trailing whitespaces in "blob"
  checkEquals(config$parameter5, "") ## empty parameter
}

test.checkConfig <- function() {
  
  conf_file <- "test-data/test_config.txt"
  ## read and parse default config  
  setupTestFramework(conf_file)
  checkTrue(checkConfig(),"checkConfig() returns true on valid config")

  ## setupTestFramework already does the checkConfig for us, so we have to wrap
  ## it into the try()
  
  ## test with unspecified parameter
  observed <- try(
                  setupTestFramework(conf_file, config.update=list(paired_ends="")),
                  silent=TRUE)
  checkEquals(class(observed), "try-error", "checkConfig() detects missing required param")

  ## test with specified, but empty parameter,
  ## this should behave the same as a missing param, as the space is removed by loadConfig()
  observed <- try(
                  setupTestFramework(conf_file, config.update=list(paired_ends=" ")),
                  silent=TRUE)
  checkEquals(class(observed), "try-error", "checkConfig() detects empty required param")

  ## test with bogus parameter
  observed <- try(
                  setupTestFramework(conf_file, config.update=list(blob=TRUE)),
                  silent=TRUE)
  checkEquals(class(observed), "try-error", "checkConfig() detects bogus parameter")

  ## test with non-existent template config file
  observed <- try(
                  setupTestFramework(conf_file, config.update=list(template_config="_blob.txt")),
                  silent=TRUE)
  checkEquals(class(observed), "try-error", "checkConfig() detects non-existent template config file")
}

test.checkConfig.alignReads <- function() {
  
  conf_file <- "test-data/test_config.txt"

  ## test a valid config
  update.list <- list(num_cores='4',
                      alignReads.nbthreads_perchunk='4')

  setupTestFramework(conf_file, config.update=update.list)
  ## NULL is good, as errors will raise a stop()
  checkEquals(HTSeqGenie:::checkConfig.alignReads(), NULL,
              "checkConfig.alignReads() returns true on valid config")

  ## nbthreads undefined triggers magic
  update.list <- list(num_cores=2,
                      alignReads.nbthreads_perchunk='')

  setupTestFramework(conf_file, config.update=update.list)
  checkEquals(HTSeqGenie:::checkConfig.alignReads(), NULL,
              "checkConfig.alignReads() returns true on valid config")
  checkEquals(getConfig.integer("alignReads.nbthreads_perchunk"), 2,
              "checkConfig.alignReads() sets alignReads.nbthreads_perchunk if undefined")

  ## now check invalid ones
  ## setupTestFramework already calls checkConfig which in turn calls checkConfig.alignReads
  ## for us, so we can't call it directly
  ## Since a failure will throw an error, we wrap it into the try()  

  ## num_cores < nbthreads
  update.list <- list(num_cores=4,
                      alignReads.nbthreads_perchunk=5)
  observed <- try(
                  setupTestFramework(conf_file, config.update=update.list),
                  silent=TRUE)
  checkEquals(class(observed), "try-error", "checkConfig() detects incompatible alignReads.nbthreads_perchunk setting")

  ## nbthreads smaller than 1
  update.list <- list(num_cores=4,
                      alignReads.nbthreads_perchunk=0)
  observed <- try(
                  setupTestFramework(conf_file, config.update=update.list),
                  silent=TRUE)
  checkEquals(class(observed), "try-error", "checkConfig() detects incompatible alignReads.nbthreads_perchunk setting")
}

test.loadConfig <- function() {
  ## when called without args, loads the default
  loadConfig()
  
  ## check that at least one param is set
  checkEquals(getConfig('template_config'), "default-config.txt",
              "loadConfig() loads default config")
  
  ## check that at least one param is set
  loadConfig(getPackageFile("test-data/test_config_single_end.txt"))
  checkTrue(!getConfig.logical('paired_ends'),
            "loadConfig() updates default config with given file")
}

test.updateConfig <- function() {
  loadConfig()

  updateConfig(list(template_config="foo.bar"))
  checkEquals(getConfig("template_config"),"foo.bar",
              "updateConfig() overwrites default config")
}

test.getConfig <- function() {
  loadConfig()
  updateConfig(list(blob="2", foo="", ming="FALSE",
                    multi='a, b, c '))

  ## simple checks
  checkEquals(getConfig("blob"), "2")
  checkEquals(getConfig("foo"), NULL)
  checkEquals(class(try(getConfig("bar"), silent=TRUE)), "try-error")
  checkEquals(getConfig.logical("foo"), NULL)
  checkEquals(getConfig.logical("ming"), FALSE)
  checkEquals(class(try(getConfig.logical("blob"), silent=TRUE)), "try-error")
  checkEquals(getConfig.numeric("foo"), NULL)
  checkEquals(getConfig.integer("blob"), 2)
  checkEquals(getConfig.integer("foo"), NULL)
  checkEquals(class(try(getConfig.integer("ming"), silent=TRUE)), "try-error")
  checkEquals(getConfig("ming", stop.ifempty=TRUE), "FALSE")
  checkEquals(class(try(getConfig("foo", stop.ifempty=TRUE), silent=TRUE)), "try-error")
  checkEquals(getConfig.vector('multi'), c('a','b','c'))
}
