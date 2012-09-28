createSummaryCounts <- function(counts, save_dir, prepend_str) {
  ## summarize nbreads per featurename
  summary_counts <- list()
  for (feature_name in names(counts$nbreads)) { 
    z <- list(counts$nbreads[[feature_name]])
    names(z) <- paste(feature_name, "nbreads", sep="_")
    summary_counts  <- c(summary_counts , z)
  }

  ## summarize nbfeatures per featurename
  for (feature_name in names(counts$nbreads.byfeature)) {
    z <- sum(counts$nbreads.byfeature[[feature_name]]>0)
    names(z) <- paste(feature_name, "nbfeatures", sep="_")
    summary_counts  <- c(summary_counts , z)
  }

  ## write file
  df <- data.frame(name=names(summary_counts), value=unlist(summary_counts)) 
  saveWithID(df, "summary_counts", id=prepend_str, save_dir=file.path(save_dir, "results"), format="tab")  
}

mergeSummaryCounts <- function(dirs_to_merge, merged_dir, prepend_str, merge_nbfeatures) {
  ## load summary tables
  lsummary_counts <- lapply(dirs_to_merge, function(dir) {
    df <- getTabDataFromFile(dir, "summary_counts")
    summary_counts <- as.numeric(df$value)
    names(summary_counts) <- df$name
    summary_counts
  })

  ## merge nbreads by summing
  summary_counts <- do.call(cbind, lsummary_counts)
  summary_counts <- apply(summary_counts, 1, sum)

  ## use merge_nbfeatures to merge nbfeatures
  summary_counts <- summary_counts[-grep("nbfeatures", names(summary_counts))]
  summary_counts <- c(summary_counts, merge_nbfeatures)
  
  ## save merged summary table
  df <- data.frame(name=names(summary_counts), value=unlist(summary_counts))
  saveWithID(df, "summary_counts", id=prepend_str, save_dir=file.path(merged_dir, "results"), format="tab")
}
