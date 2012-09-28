test.calculateTargetLength <- function(){
  config.filename <- "test-data/test_config.txt"  
  save.dir <- setupTestFramework(config.filename=config.filename, testname="test.calcTargetLengths")

  bam.file <- getPackageFile("test-data/bam_files_unit_test_PE/test_pe.gsnap.merged.concordant_uniq.bam")
  observed.table = calculateTargetLengths(bam.file, save.dir)

  ## repeat 2 times, to simulate paired reads
  isizes.from.bam <- rep(c(822,
                       15949,
                       188,
                       1016,
                       250,
                       233,
                       190,
                       9927,
                       25527,
                       234,
                       250,
                       338,
                       11181,
                       220), 2)

  expected.summary <- summary(abs(isizes.from.bam))

  checkTrue(all(observed.table['Full_Data_Set'] == expected.summary),
             'calculateTargetLength() computes correct summary')

  ## we trust summary() to do its jobs, so only test that we get something back that looks like summary output
  checkEquals(names(observed.table), c('Trimmed_Data_Set', 'Full_Data_Set'),
            "calculateTargetLengths() returns summary table")

  checkTrue(file.exists(file.path(save.dir,'reports','images','TargetLengths.pdf')),
            'calculateTargetLengths() writes pdf')
  checkTrue(file.exists(file.path(save.dir,'reports','images','TargetLengths.png')),
            'calculateTargetLengths() writes png')
}
