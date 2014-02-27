library("HTSeqGenie")

## setting this env will let the GATK tests allow to run
gatk <- Sys.getenv("GATK_PATH")
if(file.exists(gatk)) options(gatk.path=gatk)

source(getPackageFile("unitTests/runTests.R"))
