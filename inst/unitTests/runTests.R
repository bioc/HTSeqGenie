## load library
library("HTSeqGenie")
library("RUnit")
library(gmapR) # to get the TP53genome

## reload sources
resource()

## run test suite
testSuite <- defineTestSuite(name=paste("HTSeqGenie unit testing"),
                             dirs=getPackageFile("unitTests", "HTSeqGenie"))
tests <- runTestSuite(testSuite)
outputFilename <- file.path(getwd(), "HTSeqGenie-runit.txt")
printTextProtocol(tests, fileName=outputFilename)
err <- getErrors(tests)

## conclude
if( (err$nFail + err$nErr) > 0) {
  cat(paste(readLines(outputFilename), collapse="\n"), "\n\n")
  stop("runTests.R: failed to run unit tests, see ", outputFilename," for details")
} else {
  cat(paste(readLines(outputFilename), collapse="\n"), "\n\n")
  cat("runTests.R: OK !\n")
}
