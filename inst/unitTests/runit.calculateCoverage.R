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

test.isSparse <- function() {
  ## obvious cases
  x <- RleList('1' = c(rep(0,10), rep(1,100), rep(0,10)),
               '2' = c(rep(0,100), rep(1,100), rep(2,100)))
  checkTrue(HTSeqGenie:::isSparse(x))

  y <- RleList('1' = c(1:10))
  checkTrue(!HTSeqGenie:::isSparse(y))

  ## 11 runs, total length 100, thus not sparse
  y <- RleList('1' = c(c(1:10), rep(1,90)))
  checkTrue(!HTSeqGenie:::isSparse(y))

  ## 10 runs, total length 100, thus not sparse (edge case)
  y <- RleList('1' = c(c(1:9), rep(1,91)))
  checkTrue(!HTSeqGenie:::isSparse(y))

  ## 10 runs, total length 100, thus sparse (edge case)
  y <- RleList('1' = c(c(1:9), rep(1,92)))
  checkTrue(HTSeqGenie:::isSparse(y))
}

test.mergeCoverage <- function (){
  bam <- getPackageFile("test-data/variant_calling/tp53_test.bam")
  cov1 <- computeCoverage(bam, paired_ends=TRUE, extendReads=FALSE, fragmentLength=NULL)

  tmp <- HTSeqGenie:::createTmpDir()  
  path1 <- file.path(tmp, "S1", "results", "coverage.RData")
  dir.create(dirname(path1), recursive=TRUE) 
  save(cov1, file=path1)

  outpath <- file.path(tmp, "merged")
  dir.create(file.path(outpath,"results"), recursive=TRUE)
  
  mergeCoverage(c(file.path(tmp, c("S1", "S1"))),
                outdir=outpath, prepend_str="bla")

  observed <- get(load(file.path(outpath, "results/bla.coverage.RData")))
  checkEquals(hashCoverage(observed), 846378)
}

test.mergeCoverage.sparse <- function (){

  ## our coveages are always SimpleRleList,
  ## the compressed ones are not able to store the genome in 32 bit vector
  cov1 <- as(RleList('1' = RleList('1' = c(1:10))), "SimpleRleList")
  
  tmp <- HTSeqGenie:::createTmpDir()  
  path1 <- file.path(tmp, "S1", "results", "coverage.RData")
  dir.create(dirname(path1), recursive=TRUE)
  save(cov1, file=path1)
  
  outpath <- file.path(tmp, "merged")
  dir.create(file.path(outpath,"results"), recursive=TRUE)
  
  mergeCoverage(c(file.path(tmp, c("S1", "S1", "S1"))),
                outdir=outpath, prepend_str="bla")
  
  observed <- get(load(file.path(outpath, "results/bla.coverage.RData")))

  checkTrue(all((cov1*3L)==observed))
  ## They are not equal, as observed has a seqinfo slot
  ## checkEquals(observed, cov1*3L)

  unlink(file.path(outpath, "results/bla.coverage.RData"))
}
