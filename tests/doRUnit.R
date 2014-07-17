library("HTSeqGenie")

## setting this env will let the GATK tests allow to run
gatk <- Sys.getenv("GATK_PATH")
if(file.exists(gatk)) options(gatk.path=gatk)

picard <- Sys.getenv("PICARD_PATH")
if(file.exists(picard)) options(picard.path=picard)

source(getPackageFile("unitTests/runTests.R"))
