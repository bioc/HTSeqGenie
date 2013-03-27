test.wrap.callVariants <- function(){
  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 config.update=list(num_cores=1),
                                 testname="test.wrap.callVariants")

  ## create a mock summary file
  summary.preprocess <- list(processed_min_read_length=100,
                             processed_max_read_length=100)
  
  prepend_str <- getConfig("prepend_str")
  df <- data.frame(name=names(summary.preprocess), value=unlist(summary.preprocess))
  saveWithID(df, "summary_preprocess", id=prepend_str, save_dir=file.path(save.dir, "results"), format="tab")
  
  bam <- getPackageFile("test-data/variant_calling/tp53_test.bam")
  observed <- wrap.callVariants(bam)
  ## only check if we actually get something, more detailed checks should be done in VariantTools

  checkEquals(length(observed$filtered.variants), 62,
              "run.callVariants() reports correct number of variants")
  checkTrue(file.exists(file.path(save.dir, "results", "test_pe.filtered_variants_granges.RData")),
            "wrap.callVariants writes filtered variants file")
  checkTrue(file.exists(file.path(save.dir, "results", "test_pe.raw_variants_granges.RData")),
            "wrap.callVariants writes raw variants file")                      
}

test.wrap.callVariants.parallel <- function(){
  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 config.update=list(num_cores=4),
                                 testname="test.wrap.callVariants")
  ## create a mock summary file
  summary.preprocess <- list(processed_min_read_length=100,
                             processed_max_read_length=100)
  
  prepend_str <- getConfig("prepend_str")
  df <- data.frame(name=names(summary.preprocess), value=unlist(summary.preprocess))
  saveWithID(df, "summary_preprocess", id=prepend_str, save_dir=file.path(save.dir, "results"), format="tab")

  bam <- getPackageFile("test-data/variant_calling/tp53_test.bam")
  observed <- wrap.callVariants(bam)

  ## only check if we actually get something, more detailed checks should be done in VariantTools
  checkEquals(length(observed$filtered.variants), 62,
              "wrap.callVariants() runs in parallel")
}

test.writeVCF <-  function(){
  config.filename <- "test-data/test_config.txt"
  save.dir <- setupTestFramework(config.filename=config.filename, testname="test.writeVCF")
  gr = GRanges(
    ranges=IRanges(start=1, end=2), seqnames="1",
    alt=DNAStringSet('g'),
    count=2, count.ref='0', count.total='2',
    ref='g', location='test')

  writeVCF(gr)
  checkTrue(file.exists(file.path(save.dir, "results",
                                  "test_pe_variants.vcf.gz")),
            "writeVCF() write vcf file")
}

test.writeVCF.empty <-  function(){
  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 testname="test.writeVCF")
  gr = GRanges()
  checkException(writeVCF(gr),
                 "writeVCF() throws exception for empty input",
                 silent=TRUE)
}

test.plotMutationMatrix <-  function(){
  config.filename <- "test-data/test_config.txt"
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 testname="test.plotMutationMatrix")
  gr = GRanges(
    ranges = IRanges(1:2, width=c(1,1)),
    seqnames = c("a", "b"),
    alt = DNAStringSet(c('C','G')),
    ref = DNAStringSet(c('G','G')),
    count = c(1,2) )
  
  obs.matrix <- plotMutationMatrix(gr)
  checkEquals(obs.matrix[1,1], 0, "plotMutationMatrix() returns correct matrix")
  checkEquals(obs.matrix[3,2], 1, "plotMutationMatrix() returns correct matrix")
  checkEquals(obs.matrix[3,3], 1, "plotMutationMatrix() returns correct matrix")

  checkTrue(file.exists(file.path(save.dir, "reports", "images",
                                  "mutation_matrix.png")),
            "plotMutationMatrix() writes png file")
}

test.plotMutationMatrix.oneVariantOnly <-  function(){
  config.filename <- "test-data/test_config.txt"
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 testname="test.plotMutationMatrix.oneVariantOnly")
  gr = GRanges(
    ranges = IRanges(1:1, width=c(1)),
    seqnames = c("a"),
    alt = DNAStringSet(c('C')),
    ref = DNAStringSet(c('G')),
    count = c(1) )
  
  obs.matrix <- plotMutationMatrix(gr)
  checkEquals(obs.matrix[1,1], 0, "plotMutationMatrix() returns correct matrix")
  checkEquals(obs.matrix[3,2], 1, "plotMutationMatrix() returns correct matrix")

  checkTrue(file.exists(file.path(save.dir, "reports", "images",
                                  "mutation_matrix.png")),
            "plotMutationMatrix() writes png file")
}

test.processVariants <- function(){
  config.filename <- "test-data/test_config.txt"
  save.dir <- setupTestFramework(config.filename=config.filename, testname="test.processVariants")

  gr = GRanges(
    ranges = IRanges(1:2, width=c(1,1)),
    seqnames = c("a", "b"),
    alt = DNAStringSet(c('C','G')),
    ref = DNAStringSet(c('G','G')),
    count = c(1,2),
    count.pos = c(1,1),
    count.neg = c(0,1)
    )

  mock.variants <- list(raw.variants = gr,
                        filtered.variants = gr)
           
  reports.dir <- file.path(save.dir, "reports")
  observed.html <- processVariants(mock.variants, reports.dir)
  checkTrue(file.exists(observed.html), "processVariants() write html")
}
  
test.writeVariantSummary <- function(){
  config.filename <- "test-data/test_config.txt"
  save.dir <- setupTestFramework(config.filename=config.filename, testname="test.writeVariantSummary")

  gr = GRanges(
    ranges = IRanges(1:2, width=c(1,1)),
    seqnames = c("a", "b"),
    alt = DNAStringSet(c('C','G')),
    ref = DNAStringSet(c('G','G')),
    count = c(1,2),
    count.pos = c(1,1),
    count.neg = c(0,1)
    )
  mock.variants <- list(raw.variants = gr,
                        filtered.variants = gr)
  
  observed.summary <- HTSeqGenie:::.writeVariantSummary(mock.variants, 0.4)
  checkEquals(observed.summary$raw_variants_count, 2, ".writeVariantSummary reports variant count")
  checkEquals(observed.summary$filtered_variants_count, 2, ".writeVariantSummary reports variant count")

  results.dir <- file.path(save.dir, "results")
  summary.file = dir(pattern="*summary_variants.tab", path=results.dir, full.names=TRUE)
  checkTrue(file.exists(summary.file), ".writeVariantSummary() writes summary file")
}
