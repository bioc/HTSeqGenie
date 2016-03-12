test.wrap.callVariants <- function(){
  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 config.update=list(num_cores=1),
                                 testname="test.wrap.callVariants")
  
  bam <- getPackageFile("test-data/variant_calling/tp53_test.bam")
  observed <- wrap.callVariants(bam)
  
  ## only check if we actually get something, more detailed checks should be done in VariantTools
  checkEquals(length(observed$filtered.variants), 77L,
              "run.callVariants() reports correct number of variants")
  checkEquals(sum(VariantAnnotation::isIndel(observed$filtered.variants)), 14,
              "run.callVariants() reports correct number of indels")
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
  checkEquals(length(observed$filtered.variants), 77L,
              "wrap.callVariants() runs in parallel")
}

notest.wrap.callVariants.rmsk_dbsnp <- function(){

  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 config.update=list(
                                   analyzeVariants.rep_mask = getPackageFile("test-data/tp53_rmsk.bed"),
                                   analyzeVariants.dbsnp = getPackageFile("test-data/tp53_dbsnp.Rdata"),
                                   analyzeVariants.postFilters = 'avgNborCount'
                                   ),
                                 testname="test.wrap.callVariants")
  
  bam <- getPackageFile("test-data/variant_calling/tp53_test.bam")
  observed <- wrap.callVariants(bam)
  
  ## only check if we actually get something, more detailed checks should be done in VariantTools
  checkTrue('avgNborCount' %in% names(VariantAnnotation::hardFilters(observed$filtered.variants)),
            "AvgNborCount filter is activated")
  ##AvgNbor filter removes one snp
  checkEquals(length(observed$filtered.variants), 74,
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

test.wrap.callVariants.which <- function(){
  ## build analyzeVariant position file
  gr <- TP53Which()
  whichFile <- tempfile()
  saveRDS(gr, whichFile)

  ## run variant calling
  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 config.update=list(
                                   analyzeVariants.positions=whichFile
                                   ),
                                 testname="test.wrap.callVariants.which")
  
  bam <- getPackageFile("test-data/variant_calling/tp53_test.bam")
  observed <- wrap.callVariants(bam)
 
  checkEquals(length(observed$filtered.variants), 5,
              "run.callVariants() reports correct number of variants")
}

test.writeVCF.vcfStat <-  function(){
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

  filename <- file.path(save.dir, "results",
                        "test_pe.variants.vcf")
  vcf.filename <- writeVCF(vr, filename=filename)
  checkTrue(file.exists(vcf.filename),
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
  checkEquals(writeVCF(gr, file.path(save.dir,"bla.vcf.gz")), NULL, 
                 "writeVCF() returns NULL")
}

##Ideally I'd like to test the genotype function on a single variant,
## but VT does not seem to like it, so I give up for now.
not.yet.test.callGenotypes <- function() {
  config.filename <- "test-data/test_config.txt"
  save.dir <- setupTestFramework(config.filename=config.filename, testname="test.callGenotypes")

  genome.name <- gmapR:::geneGenomeName('TP53')
 
  vr = VRanges(
    ranges=IRanges(start=50, end=50),
    seqnames=Rle(c('TP53')),
    seqinfo = Seqinfo("TP53", 2025772, NA, genome.name),
    alt='G',
    ref='A',
    altDepth=c(5),
    refDepth=c(1),
    totalDepth=c(6),
    location='TP53',
    sampleNames = 'foo'
    )
 
  cov <- as(RleList('TP53' = RleList('1' = c(1:2025772))), "SimpleRleList")
  HTSeqGenie:::saveCoverage(cov, save.dir, "test_pe")

  observed <- HTSeqGenie:::.callGenotypes(vr)
}

test.callVariantsVariantTools.genotype <- function(){
  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename,
                                 config.update=list(num_cores=1,
                                   analyzeVariants.genotype=TRUE),
                                 testname="test.genotype")
  
  bam <- getPackageFile("test-data/variant_calling/tp53_test.bam")
  cov <- computeCoverage(bam)
  HTSeqGenie:::saveCoverage(cov, save.dir, "test_pe")

  ## When analyzeVariants.genotype=TRUE, this function
  ## creates a genotypes.vcf.gz file as side effect, for which we test 
  observed <- HTSeqGenie:::callVariantsVariantTools(bam)
  geno.file <- sub("variants", "genotypes", observed)

  checkTrue(file.exists(geno.file), "The genotypes.vcf file is created")
  checkTrue( length(geno(readVcf(geno.file, genome='bla'))$GT) > 1,
            "and it has some genotypes annotated")
}

test.annotateVariants <- function(){
  if (HTSeqGenie:::checkVEP()){
    vcf <- getPackageFile("test-data/braf_v600.vcf")
    config.filename <- "test-data/test_config.txt"  
    save.dir <- setupTestFramework(config.filename=config.filename,
                                   config.update=list(num_cores=2,
                                     analyzeVariants.vep_options="--vcf --sift b --polyphen b --protein --check_existing --check_alleles --symbol --canonical --ccds --maf_1kg --maf_esp --cache --force_overwrite --no_stats"),
                                   testname="test.annotateVariants")
    
    anno.file <- annotateVariants(vcf)
    checkTrue(file.exists(anno.file))
    checkEquals(1, length(info(readVcf(anno.file, genome="GRCh38"))$CSQ),
                "vcf contains consequence info")
  } else {
    DEACTIVATED("Skipped annotateVariants() test")
  }
}
