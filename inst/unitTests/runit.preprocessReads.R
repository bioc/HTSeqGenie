test.preprocessReads <- function() {
  ## setup test framework
  setupTestFramework(config.filename="test-data/test_config.txt", testname="test.preprocessReads")

  ## test preprocessReads
  summary_preprocess <- preprocessReads()

  ## checks
  checkEqualsNumeric(summary_preprocess$total_reads, 26,
                     "preprocessReads() does not return the correct number of reads")
  checkEquals(summary_preprocess$highqual_reads, 25,
              "preprocessReads() does not return the correct number of highqual reads")
  checkEqualsNumeric(summary_preprocess$adapter_contam, 1,
                     "preprocessReads() does not return the correct number of adapter contam")
  checkEqualsNumeric(summary_preprocess$read_length, 75,
                     "preprocessReads() does not return the correct read length")
  ## Since rRNA contam check is currently disbaled this will be now 0 instead of 1
  checkEquals(summary_preprocess$rRNA_contam_reads, 0,
              "preprocessReads() does not return the correct number of rRNA contams")
}

test.preprocessReads_single_end <- function() {
  ## setup test framework
  setupTestFramework(config.filename="test-data/test_config_single_end.txt",
                     config.update =list(
                       detectAdapterContam.force_paired_end_adapter=TRUE),
                     testname="test.preprocessReads_single_end")

  ## test preprocessReads
  summary_preprocess <- preprocessReads()

  ## checks
  checkEqualsNumeric(summary_preprocess$total_reads, 26,
                     "preprocessReads() does not return the correct number of reads")
  checkEquals(summary_preprocess$highqual_reads, 25,
              "preprocessReads() does not return the correct number of highqual reads")
  checkEqualsNumeric(summary_preprocess$adapter_contam, 1,
                     "preprocessReads() does not return the correct number of adapter contam")
  checkEqualsNumeric(summary_preprocess$read_length, 75,
                     "preprocessReads() does not return the correct read length")
  checkEquals(summary_preprocess$rRNA_contam_reads, 0,
              "preprocessReads() does not return the correct number of rRNA contams")
}

test.preprocessReads.minichunks <- function() {
  ## setup test framework
  config.update <- list(chunk_size=10, detectRRNA.do=FALSE)
  setupTestFramework(config.filename="test-data/test_config.txt", config.update=config.update,
                     testname="test.preprocessReads.minichunks")
  
  ## test preprocessReads
  summary_preprocess <- preprocessReads()

  ## checks
  checkEqualsNumeric(summary_preprocess$total_reads, 26,
                     "preprocessReads() does not return the correct number of reads")
  checkEquals(summary_preprocess$highqual_reads, 25,
              "preprocessReads() does not return the correct number of highqual reads")
  checkEqualsNumeric(summary_preprocess$adapter_contam, 1,
                     "preprocessReads() does not return the correct number of adapter contam")
  checkEqualsNumeric(summary_preprocess$read_length, 75,
                     "preprocessReads() does not return the correct read length")
  checkEquals(summary_preprocess$rRNA_contam_reads, 0,
              "preprocessReads() does not return 0 rRNA contams (detectRRNA.do is disabled)")
}
