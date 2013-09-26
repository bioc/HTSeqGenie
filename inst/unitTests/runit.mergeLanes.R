test.mergeLanes <- function() {
  ## set framework
  setupTestFramework(config.filename="test-data/test_config.txt",
                     config.update=list(remove_chunkdir=FALSE, chunk_size=10),
                     testname="test.mergeLanes")
  config_filename <- file.path(tempdir(), "config.txt")
  writeConfig(config_filename)
  
  ## run pipeline
  runPipelineConfig(config_filename)
  maindir <- getConfig("save_dir")
  gcounts1 <- getTabDataFromFile(maindir, "counts_gene")$count
  coverage1 <- get(load(getObjectFilename(file.path(maindir, "results"), "coverage.RData$")))
  variants1 <- get(load(getObjectFilename(file.path(maindir, "results"), "filtered_variants.RData$")))
  
  ## run mergeLanes on chunks
  indirs <- getChunkDirs()
  mergedir <- file.path(maindir, "merged")
  prepend_str <- "merged"
  ## preMergeChecks.do is disabled since variants are not computed on chunks, and ID verify cannot be performed
  mergeLanes(indirs, mergedir, prepend_str, num_cores=1, preMergeChecks.do=FALSE, config_update=list(shortReadReport.do=FALSE))
  
  gcounts2 <- getTabDataFromFile(mergedir, "counts_gene")$count
  coverage2 <- get(load(getObjectFilename(file.path(mergedir, "results"), "coverage.RData$")))
  variants2 <- get(load(getObjectFilename(file.path(mergedir, "results"), "filtered_variants.RData$")))

  ## compare results
  checkEquals(gcounts1, gcounts2, "gene counts are the same after merging")
  checkTrue(hashCoverage(coverage1)==hashCoverage(coverage2), "coverage is the same after merging")
  checkEquals(length(variants1), length(variants2), "varisnts are the same after merging")
  
  ## compare summaries
  summaries <- paste("summary", c('preprocess', 'alignment', 'counts'),
                      sep="_")

  x <- lapply(summaries, function(summary){
    checkTrue(all(parseSummaries(maindir, summary) ==
                  parseSummaries(mergedir, summary)),
              paste(summary, "is identical", sep=" "))
  })
}
