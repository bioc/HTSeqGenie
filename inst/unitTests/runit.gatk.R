test.checkGATKJar <- function(){

  checkTrue(!checkGATKJar(path="foo"),
            "checkGATKJar() returns false as expected")

  path <- getOption("gatk.path")
  if(!is.null(path)){
    ## Note this tests relies on having the gatk.path
    ## correctly set in the options.
    ## In general this will no be the tests for automated tests,
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
    ##TODO: fix tp53 demo fasta file
    ## get-genome -D .local/share/gmap/  -d TP53 <- demo -A > tp53 <- genome.fasta
    
##    tp53.fasta <- system.file("extdata/hg19.p53.fasta", package="gmapR")
    save.dir <- setupTestFramework(config.filename=config.filename,
                                   config.update=list(num_cores=1,
                                     analyzeVariants.method="GATK",
                                     path.gatk=getOption('gatk.path'),
                                     path.gatk_genomes=getPackageFile("test-data/variant_calling/")),  
                                   testname="test.callVariantsGATK")
    
    bam.file <- getPackageFile("test-data/variant_calling/tp53_test.bam")

    callVariantsGATK(bam.file)
    checkTrue(file.exists(file.path(save.dir, "results", "test_pe_variants.vcf.gz")),
              "callVariantsGATK writes vcf file")
    checkTrue(file.exists(file.path(save.dir, "results", "test_pe_variants.vcf.gz.tbi")),
              "callVariantsGATK writes vcf index file")
}
  else{
    DEACTIVATED("callVariantsGATK() tests need gatk.path option set")
  }
}
  

  
