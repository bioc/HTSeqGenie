test.storeCoverage <- function() {
  ## set up test
  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename, testname="test.storeCoverage")

  ## compute coverage
  bam <- getPackageFile("test-data/variant_calling/test_set.merged.uniques.bam")
  file <- storeCoverage(bam)

  observed <- get(load(file))
  
  checkEquals(length(observed), 84,
              "storeCoverage() reports coverage for all 84 genomic fragments")

  checkTrue(file.exists(file.path(save.dir, "results", "test_pe.coverage.RData")),
            "storeCoverage() writes coverage file")

  checkTrue(file.exists(file.path(save.dir, "results", "test_pe.coverage.bw")),
            "storeCoverage() writes bigwig coverage file")
}

test.storeCoverage.extension <- function() {
  ## set up test
  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename, testname="test.storeCoverage")

  ## compute coverage with different options
  bam <- getPackageFile("test-data/variant_calling/test_set.merged.uniques.bam")
  cov0p <- get(load(storeCoverage(bam, extendReads=FALSE, paired_ends=TRUE)))
  covxp <- get(load(storeCoverage(bam, extendReads=TRUE, paired_ends=TRUE)))
  cov0s <- get(load(storeCoverage(bam, extendReads=FALSE, paired_ends=FALSE)))
  covxs <- get(load(storeCoverage(bam, extendReads=TRUE, paired_ends=FALSE)))
  covxs150 <- get(load(storeCoverage(bam, extendReads=TRUE, paired_ends=FALSE, fragmentLength=150)))
}
