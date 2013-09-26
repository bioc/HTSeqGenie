## full tophalf test with 3 chunks in parallel
test.alignReads <- function() {
  ## setup test framework
  setupTestFramework(config.filename="test-data/test_config.txt",
                     config.update=list(chunk_size=10),
                     testname="test.alignReads")
  preprocessReads()
  alignReads()

  ## correct number of chunks
  chunk_dirs <- getChunkDirs()
  checkEqualsNumeric(length(chunk_dirs), 3,
                     "alignReads() produces the expected number of chunks")

  ## peek into at least one dir to see if it actually computed something
  bam_files = dir(file.path(chunk_dirs[1], "bams"), pattern="*.bam$")
  checkEqualsNumeric(length(bam_files), 19,
                     "alignReads() produces 19 bam files")

  unlink(getConfig("save_dir"), recursive=TRUE)
}

test.alignReadsOneSingleEnd <- function() {
  ## setup test framework
  setupTestFramework(config.filename="test-data/test_config_single_end.txt",
                     config.update=list(num_cores=1, prepend_str="test.alignReads",
                       alignReads.max_mismatches="0"),
                     testname="test.alignReadsOneSingleEnd")

  HTSeqGenie:::alignReadsChunk(getPackageFile("test-data/unit_tests_1.fastq"))
  
  bam_dir <- file.path(getConfig("save_dir"), "bams")
  
  mapping_types= c('unpaired_uniq', 'unpaired_transloc',
                   'unpaired_mult', 'nomapping')

  base <- paste(bam_dir,"test.alignReads", sep='/')
  expected_files <- paste(base, mapping_types, "bam", sep=".")
  checkTrue(all(file.exists(expected_files)),
            "alignReads() creates all expected files for single ends")

  expected_files <- paste(base, mapping_types, "bam.bai", sep=".")
  checkTrue(all(file.exists(expected_files)),
            "alignReads() creates all expected index files for single ends")

  unlink(getConfig("save_dir"), recursive=TRUE)
}

test.alignReads.sparsechunks <- function() {
  ## setup framework, to generate sparse chunks (3 reads out of 26, using chunk_size=8)
  save_dir <- setupTestFramework(config.filename="test-data/test_config.txt"  ,
                                 config.update=list(chunk_size=8, subsample_nbreads=2),
                                 testname="test.alignReads.sparsechunks")
  ## run pipeline
  preprocessReads()

  ## check results
  summary_preprocess <- getNumericVectorDataFromFile(save_dir, "summary_preprocess")
  checkEqualsNumeric(summary_preprocess["total_reads"], 2)
  checkEqualsNumeric(sum(is.na(summary_preprocess)), 0, "no NAs were detected in summary_preprocess")
}

