test.computeCoverage <- function() {
  ## set up test
  bam <- getPackageFile("test-data/variant_calling/test_set.merged.uniques.bam")

  ## compute coverage with different options
  cov1 <- computeCoverage(bam, paired_ends=TRUE, extendReads=FALSE, fragmentLength=NULL)
  cov2 <- computeCoverage(bam, paired_ends=TRUE, extendReads=TRUE, fragmentLength=NULL)
  cov3 <- computeCoverage(bam, paired_ends=FALSE, extendReads=FALSE, fragmentLength=NULL)
  cov4 <- computeCoverage(bam, paired_ends=FALSE, extendReads=TRUE, fragmentLength=NULL)
  cov5 <- computeCoverage(bam, paired_ends=FALSE, extendReads=TRUE, fragmentLength=150)
  cov6 <- computeCoverage(bam, paired_ends=TRUE, extendReads=TRUE, maxFragmentLength=1e4)
  
  ## check length OK
  checkEquals(length(cov1), 84,
              "computeCoverage() reports coverage for all 84 genomic fragments")

  ## check correct values
  checkEquals(hashCoverage(cov1), 890968)
  checkEquals(hashCoverage(cov2), 895658)
  checkEquals(hashCoverage(cov3), 890968)
  checkEquals(hashCoverage(cov4), 1302696)
  checkEquals(hashCoverage(cov5), 1760400)
  checkEquals(hashCoverage(cov6), 867753)

  ## check correct fragmentLengths
  checkEquals(attr(cov4, "fragmentLength"), 111)
  checkEquals(attr(cov5, "fragmentLength"), 150)       
}
