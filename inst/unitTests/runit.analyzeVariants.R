test.wrap.callVariants <- function(){
  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 config.update=list(num_cores=1),
                                 testname="test.wrap.callVariants")
  
  bam <- getPackageFile("test-data/variant_calling/tp53_test.bam")
  observed <- wrap.callVariants(bam)
  
  ## only check if we actually get something, more detailed checks should be done in VariantTools
  checkEquals(length(observed$filtered.variants), 61,
              "run.callVariants() reports correct number of variants")
  checkTrue(file.exists(file.path(save.dir, "results", "test_pe.filtered_variants.RData")),
            "wrap.callVariants writes filtered variants file")
  checkTrue(file.exists(file.path(save.dir, "results", "test_pe.raw_variants.RData")),
            "wrap.callVariants writes raw variants file")                      
}

test.wrap.callVariants.parallel <- function(){
  config.filename <- "test-data/test_config.txt"
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 config.update=list(num_cores=4),
                                 testname="test.wrap.callVariants.parallel")

  bam <- getPackageFile("test-data/variant_calling/tp53_test.bam")
  observed <- wrap.callVariants(bam)

  ## only check if we actually get something, more detailed checks should be done in VariantTools
  checkEquals(length(observed$filtered.variants), 61,
              "wrap.callVariants() runs in parallel")
}

test.wrap.callVariants.rmsk_dbsnp <- function(){
  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 config.update=list(
                                   analyzeVariants.rep_mask = getPackageFile("test-data/tp53_rmsk.bed"),
                                   analyzeVariants.dbsnp = getPackageFile("test-data/tp53_dbsnp.Rdata")
                                   ),
                                 testname="test.wrap.callVariants")
  
  bam <- getPackageFile("test-data/variant_calling/tp53_test.bam")
  observed <- wrap.callVariants(bam)

  ## only check if we actually get something, more detailed checks should be done in VariantTools
  checkEquals(length(observed$filtered.variants), 61,
              "run.callVariants() reports correct number of variants")
  checkTrue(file.exists(file.path(save.dir, "results", "test_pe.filtered_variants.RData")),
            "wrap.callVariants writes filtered variants file")
  checkTrue(file.exists(file.path(save.dir, "results", "test_pe.raw_variants.RData")),
            "wrap.callVariants writes raw variants file")                      
}

test.wrap.callVariants.filters <- function(){
  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 config.update=list(
                                   analyzeVariants.callingFilters='nonRef,readCount',
                                   analyzeVariants.postFilters=''),
                                 testname="test.wrap.callVariants.filters")

  bam <- getPackageFile("test-data/variant_calling/tp53_test.bam")
  observed <- wrap.callVariants(bam)
  checkEquals(names(VariantAnnotation::hardFilters(observed$filtered.variants)), c('nonRef', 'readCount'),
              "we can turn off some Filters for variant calling")

  ## now test if we can turn off all filters
  updateConfig(list(analyzeVariants.callingFilters=''))
  observed <- wrap.callVariants(bam)

  checkEquals(names(VariantAnnotation::hardFilters(observed$filtered.variants)), c(),
              "we can turn off all filkters for variant calling")
}

#de-activated until VariantAnnotatioon can handle simple VRanges
notest.writeVCF.vcfStat <-  function(){
  config.filename <- "test-data/test_config.txt"
  save.dir <- setupTestFramework(config.filename=config.filename, testname="test.writeVCF")
  vr = VRanges(
    ranges=IRanges(start=1, end=1),
    seqnames=Rle('1'),
    alt='G',
    ref='A',
    altCount=2,
    refCount='0',
    totalCount='2',
    location='test',
    sampleNames = 'foo'
    )

  vcf.filename <- writeVCF(vr)
  checkTrue(file.exists(file.path(save.dir, "results",
                                  "test_pe.variants.vcf.gz")),
            "writeVCF() write vcf file")

  ## check vcfStat
  vcfstat <- vcfStat(vcf.filename)
  checkEquals(vcfstat[["nb.variants"]], 1)
  checkEquals(vcfstat[["nb.snvs.G>A"]], 0)
  checkEquals(vcfstat[["nb.snvs.A>G"]], 1)
}

test.writeVCF.NULL <-  function(){
  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 testname="test.writeVCF")
  gr = GRanges()
  checkEquals(writeVCF(gr), NULL, 
                 "writeVCF() returns NULL")
}

