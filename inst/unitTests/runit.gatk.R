test.checkGATKJar <- function(){

  checkTrue(!checkGATKJar(path="foo"),
            "checkGATKJar() returns false as expected")

  path <- getOption("gatk.path")
  if(!is.null(path)){
    ## Note this tests relies on having the gatk.path
    ## correctly set in the options.
    ## In general this will no be the caSE for automated tests,
    ## but it is good enough for manual tests.
    checkTrue(checkGATKJar(path),
              'checkGATKJar() returns true')
  }
  else {
    DEACTIVATED("checkGATKJar() test needs gatk.path option set")
  }
}

test.gatk <- function() {
  if (checkGATKJar()){
    checkEquals(gatk(method="UnifiedGenotyper", args="--help"), 0,
              "gatk() executes jar file")

    checkException(gatk(method="Foo"), silent=TRUE, "gatk() detects error and stops")
  }
  else{
    DEACTIVATED("gatk() tests need gatk.path option set")
  }
}

test.callVariantsGATK <- function() {
  if (checkGATKJar()){
    config.filename <- "test-data/test_config.txt"

    buildTP53FastaGenome()

    save.dir <- setupTestFramework(config.filename=config.filename,
                                   config.update=list(num_cores=1,
                                     analyzeVariants.method="GATK",
                                     path.gatk=getOption('gatk.path'),
                                     path.gatk_genomes=tempdir()),
                                   testname="test.callVariantsGATK")
    
    bam.file <- getPackageFile("test-data/variant_calling/tp53_test.bam")

    HTSeqGenie:::callVariantsGATK(bam.file)
    checkTrue(file.exists(file.path(save.dir, "results", "test_pe.variants.vcf.gz")),
              "callVariantsGATK writes vcf.gz file")
    checkTrue(file.exists(file.path(save.dir, "results", "test_pe.variants.vcf.gz.tbi")),
              "callVariantsGATK writes vcf index file")
  }
  else{
    DEACTIVATED("callVariantsGATK() tests need gatk.path option set")
  }
}

test.callVariantsGATK.withFiltering <- function() {
  if (checkGATKJar()){
    config.filename <- "test-data/test_config.txt"

    buildTP53FastaGenome()

    save.dir <- setupTestFramework(config.filename=config.filename,
                                   config.update=list(num_cores=1,
                                     analyzeVariants.method="GATK",
                                     path.gatk=getOption('gatk.path'),
                                     path.gatk_genomes=tempdir(),
                                     gatk.filter_repeats=TRUE,
                                     analyzeVariants.rep_mask = getPackageFile("test-data/tp53_rmsk.bed")),
                                   testname="test.callVariantsGATK")
    
    bam.file <- getPackageFile("test-data/variant_calling/tp53_test.bam")

    HTSeqGenie:::callVariantsGATK(bam.file)
    checkTrue(file.exists(file.path(save.dir, "results", "test_pe.variants.vcf.gz")),
              "callVariantsGATK writes vcf.gz file")
    checkTrue(file.exists(file.path(save.dir, "results", "test_pe.variants.vcf.gz.tbi")),
              "callVariantsGATK writes vcf index file")
}
  else{
    DEACTIVATED("callVariantsGATK() tests need gatk.path option set")
  }
}

test.excludeVariantsByRegion <- function(){

  gr <- GRanges(seqnames=Rle('1'),
                 IRanges(1:2, width=1), strand=c('*','+'))
  
  variants = VariantAnnotation::VRanges(
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
  
  observed = HTSeqGenie:::excludeVariantsByRegions(variants, gr)
  checkEquals(length(observed), 0,
              " variant overlapping repeat on * strand is filtered out")

  
  ranges(variants) <- IRanges(start=2, end=2)
  observed = HTSeqGenie:::excludeVariantsByRegions(variants, gr)
  checkEquals(length(observed), 0,
              " variant overlapping repeat on + strand is filtered out")
  

  ranges(variants) <- IRanges(start=3, end=3)
  observed = HTSeqGenie:::excludeVariantsByRegions(variants, gr)
  checkEquals(length(observed), 1,
              " variant outside of mask is not filtered")
}

test.realignIndelsGATK <- function() {
  if (checkGATKJar()){
    config.filename <- "test-data/test_config.txt"

    buildTP53FastaGenome()
    
    save.dir <- setupTestFramework(config.filename=config.filename,
                                   config.update=list(
                                     num_cores=1,
                                     path.gatk=getOption('gatk.path'),
                                     path.gatk_genomes=tempdir() ),
                                   testname="test.realignIndelsGATK")
    
    bam.file <- getPackageFile("test-data/Realigner/jiggling_indels.bam")
    
    observed.file <- HTSeqGenie:::realignIndelsGATK(bam.file)
    checkTrue(file.exists(observed.file),
              "realignIndelsGATK() writes realigned bam file")
    checkEquals(unique(cigar(readGAlignments(observed.file)[3:4])),
                "23M1D77M", "... and indels in second read where shifted")    
  }
  else{
    DEACTIVATED("realignIndelsGATK() tests need gatk.path option set")
  }
}

test.realignIndelsGATK.parallel <- function() {
  if (checkGATKJar()){
    config.filename <- "test-data/test_config.txt"

    HTSeqGenie:::buildAnyFastaGenome(c("KRAS", "TP53"))
    ## genome name changes with each version bump of the underlying txdb
    ## this mess is hidden by using this function
    genome.name <- gmapR:::geneGenomeName('KRAS_TP53')
    
    save.dir <- setupTestFramework(config.filename=config.filename,
                                   config.update=list(
                                     alignReads.genome=genome.name,
                                     num_cores=2,
                                     path.gatk=getOption('gatk.path'),
                                     path.gatk_genomes=tempdir() ),
                                   testname="test.realignIndelsGATK.parallel")
    
    ## we run into locking issues if two GATK process try to create the dict simultanously
    file.copy(getPackageFile(paste0("test-data/Realigner/KRAS_TP53_demo.dict")),
              file.path(tempdir(), paste0(genome.name, ".dict")))
    bam.file <- getPackageFile("test-data/Realigner/jiggling_indels_tp53_kras.bam")

    observed.file <- HTSeqGenie:::realignIndelsGATK(bam.file)
    checkTrue(file.exists(observed.file),
              "realignIndelsGATK() writes realigned bam file")
    checkEquals(unique(cigar(readGAlignments(observed.file)[3:4])),
                "23M1D77M", "... and indels in second read where shifted")    
  }
  else{
    DEACTIVATED("realignIndelsGATK() tests need gatk.path option set")
  }
}

test.realignIndels <- function() {
  if (checkGATKJar()){
    config.filename <- "test-data/test_config.txt"

    buildTP53FastaGenome()

    save.dir <- setupTestFramework(config.filename=config.filename,
                                   config.update=list(num_cores=1,
                                     path.gatk=getOption('gatk.path'),
                                     path.gatk_genomes=tempdir()),
                                   testname="test.realignIndels")
    
    bam.file <- getPackageFile("test-data/Realigner/jiggling_indels.bam")
    analyzed.bam <- file.path(save.dir, "bams", "test_pe.analyzed.bam")
    file.copy(bam.file, analyzed.bam)
    indexBam(analyzed.bam)

    #test starts here
    realignIndels()
    
    # function overwrites existing files, so check for change in size as proof of action
    checkTrue(file.exists(analyzed.bam),
              "callVariantsGATK writes bam")
    checkTrue(file.info(bam.file)$size != file.info(analyzed.bam)$size,
              " and the size has changed")
    checkTrue(file.exists(paste0(analyzed.bam, '.bai')),
              "callVariantsGATK writes index file")
  }
  else{
    DEACTIVATED("test.realignIndels() tests need gatk.path option set")
  }
}

test_zipUp <- function(){

  testfile <- file.path(HTSeqGenie:::createTmpDir(),
                        "test.vcf")
  file.copy(system.file("extdata", "ex2.vcf", package="VariantAnnotation"),
            testfile)

 
  observed <- HTSeqGenie:::.zipUp(testfile)
  checkTrue(file.exists(observed), ".zipUp creates file")
  checkTrue(grepl('gz$', observed), "...and the file has gz suffix")

  observed2 <- HTSeqGenie:::.zipUp(observed)
  checkEquals(observed, observed2, "zipUp does not modify already zipped filed")
}
  
