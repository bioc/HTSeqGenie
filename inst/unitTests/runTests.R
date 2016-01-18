## load library
library("RUnit")
library("HTSeqGenie")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(gmapR) # to get the TP53genome

## setting this env will let the GATK tests allow to run
gatk <- Sys.getenv("GATK_PATH")
if(file.exists(gatk)) options(gatk.path=gatk)

picard <- Sys.getenv("PICARD_PATH")
if(file.exists(picard)) options(picard.path=picard)

## reload sources
resource()

## run test suite
testSuite <- defineTestSuite(name=paste("HTSeqGenie unit testing"),
                             dirs=getPackageFile("unitTests", "HTSeqGenie"))
tests <- runTestSuite(testSuite)
outputFilename <- file.path(getwd(), "HTSeqGenie-runit.txt")
printTextProtocol(tests, fileName=outputFilename)
printHTMLProtocol(tests, fileName=file.path(getwd(), "HTSeqGenie-runit.html"))

err <- getErrors(tests)

## show warnings
warnings()

## conclude
if( (err$nFail + err$nErr) > 0) {
  cat(paste(readLines(outputFilename), collapse="\n"), "\n\n")
  stop("runTests.R: failed to run unit tests, see ", outputFilename," for details")
} else {
  cat(paste(readLines(outputFilename), collapse="\n"), "\n\n")
  cat("runTests.R: OK !\n")
}
