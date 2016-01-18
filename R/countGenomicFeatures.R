##' Count overlaps with genomic features
##'
##' @title  Count overlaps with genomic features
##' @return Nothing
##' @author Gegoire Pau
##' @keywords internal
##' @export
##' @importFrom GenomicRanges seqnames
##' @importFrom IRanges split
countGenomicFeatures <- function() {
  safeExecute({
    loginfo("countGenomicFeatures.R/countGenomicFeatures: starting...")
    
    ## rebuild lbams from chunk files
    chunkdirs <- getChunkDirs()

    ## load "genomic_features": pre-calculated granges of exons, genes, cds, etc to count overlaps
    filename <- file.path(getConfig("path.genomic_features"), getConfig("countGenomicFeatures.gfeatures"))
    genomic_features <- get(load(filename))
    
    ## schedule each chunk
    listIterator.init(chunkdirs)
    processChunks(listIterator.next, function(chunkdir, save_dir) countGenomicFeaturesChunk(chunkdir, genomic_features))
    
    ## merge counts and summaries
    save_dir <- getConfig("save_dir")
    prepend_str <- getConfig("prepend_str")
    merge_nbfeatures <- mergeCounts(chunkdirs, save_dir, prepend_str)
    mergeSummaryCounts(chunkdirs, save_dir, prepend_str, merge_nbfeatures)   
    writeGenomicFeaturesReport()
    loginfo("countGenomicFeatures.R/countGenomicFeatures: done...")
  }, memtracer=getConfig.logical("debug.tracemem"))
}

##' Count reads by genomic Feature
##'
##' given a BAM-file output from gsnap (with the MD tag), count
##' hits to exons, genes, ncRNAs, etc.
##' and quantify miRNA/ncRNA contaminatino 
##' @title Count reads by genomic Feature
##' @param save_dir PAth to a pipeline run's save dir
##' @param genomic_features A list of genomic features to tally
##' @return Nothing
##' @author Cory Barr
##' @keywords internal
##' @export
countGenomicFeaturesChunk <- function(save_dir, genomic_features) {
  ## get config parameters
  paired_ends <- getConfig.logical("paired_ends")
  prepend_str <- getConfig("prepend_str")
  loginfo("countGenomicFeatures.R/countGenomicFeaturesChunk: starting...")
  
  ## read bam file
  bams <-  dir(file.path(save_dir, "bams"), full.names=TRUE)
  analyzedbam <- grep("analyzed\\.bam$", bams, value=TRUE)
  reads <- readRNASeqEnds(analyzedbam, paired_ends, remove.strandness=TRUE)
  
  ## count genomic featuress
  counts <- countFeatures(reads, genomic_features)

  ## compute rpkm and save counts
  nbanalyzedreads <- length(reads)
  for (feature_name in names(genomic_features)) {
    widths <- computeWidth(genomic_features[[feature_name]])
    nbreads.byfeature <- counts$nbreads.byfeature[[feature_name]]
    rpkms <- rpkm(counts=nbreads.byfeature, widths=widths, nbanalyzedreads)
    df <- data.frame(name=names(nbreads.byfeature), count=nbreads.byfeature, width=widths, rpkm=rpkms)
    filename <- paste("counts", feature_name, sep="_")
    saveWithID(df, filename, id=prepend_str, save_dir=file.path(save_dir, "results"), format="tab")  
  }
 
  ## write summary counts
  createSummaryCounts(counts, save_dir, prepend_str)
  
  loginfo("countGenomicFeatures.R/countGenomicFeaturesChunk: done") 
}

##' Calculate RPKM
##' 
##' @title Calculate RPKM
##' @param counts A vector of counts
##' @param widths vector of the width of each bin the counts were perfomred on
##' @param nbreads vector containing number of reads mapped to each bin
##' @return vector of RPKMs
##' @author Gregoire Pau
##' @export
##' @keywords internal
rpkm <- function (counts, widths, nbreads)  {
  (1e9*counts) / (nbreads * as.numeric(widths))
}

computeWidth <- function(gf) {
  if (class(gf)=="GRangesList") sapply(width(reduce(gf)), sum)
  else width(gf)
}

mergeCounts <- function(dirs_to_merge, merged_dir, prepend_str) {
  loginfo("countGenomicFeatures.R/mergeCounts: starting...")

  ## enumerate counts_* files
  rdata_dirs <- file.path(dirs_to_merge, "results")
  count_files <- lapply(rdata_dirs,
                          function(x) {
                            grep("\\.counts_.*\\.tab$",
                                 dir(x, full.names=TRUE),
                                 value=TRUE)
                          })
  lens <- elementLengths(count_files)
  if (!all(lens[1] == lens))
    stop("countGenomicFeatures.R/mergeCounts: there are a differing number of files storing hits by genomic features across the experiments")
  if (all(lens==0)) stop("countGenomicFeatures.R/mergeCounts: no counts_* files were found")
  
  ## make sure the same genomic features are being summed
  count_files <- lapply(count_files, sort)
  feature_names <- lapply(count_files,
                          function(x) {
                            sub("^.*\\.(counts_[^\\.]+).*$", "\\1", basename(x))
                          })
  if(!all(sapply(feature_names, function(i) all(i == feature_names[[1]]))))
    stop("countGenomicFeatures.R/mergeCounts: there are different sets of genomic features in the experiments to merge")
    
  ## group by genomic feature type
  files_by_gf <- split(unlist(count_files), rep(feature_names[[1]], length(count_files)))

  ## get nbanalyzedreads
  df <- getTabDataFromFile(merged_dir, "summary_alignment")
  nbanalyzedreads <- as.numeric(df$value[match("analyzed", df$name)])

  ## merge counts and save
  merge_nbfeatures <- list()
  for (feature_name in names(files_by_gf)) {
    filenames <- files_by_gf[[feature_name]]
    dfm <- NULL
    for (filename in filenames) {
      z <- try({
        df <- read.table(filename, sep="\t", header=TRUE, stringsAsFactors=FALSE)
        if (is.null(dfm)) dfm <- df
        else {
          dfm$count <- dfm$count + df$count
          if (any(dfm$name!=df$name)) {
            stop("countGenomicFeatures.R/mergeCounts: names of the genomic features across samples do not agree")
          }
          if (any(dfm$width!=df$width)) {
            stop("countGenomicFeatures.R/mergeCounts: widths of the genomic features across samples do not agree")
          }
        }
      })
      if (class(z)=="try-error") stop("countGenomicFeatures.R/mergeCounts: cannot load or merge filename=", filename)
    }

    ## recompute merge_nbfeatures, this will be needed to merge summary_counts
    z <- list(sum(dfm$count>0))
    names(z) <- paste(feature_name, "nbfeatures", sep="_")
    merge_nbfeatures <- c(merge_nbfeatures, z)
    
    ## recompute rpkm
    rpkms <- rpkm(counts=dfm$count, widths=dfm$width, nbanalyzedreads)
    df <- data.frame(name=dfm$name, count=dfm$count, width=dfm$width, rpkm=rpkms, stringsAsFactors=FALSE)
    saveWithID(df, feature_name, id=prepend_str, save_dir=file.path(merged_dir, "results"), format="tab")  
  }

  loginfo("countGenomicFeatures.R/mergeCounts: done")
  merge_nbfeatures
}


##' Given GRanges, counts number of hits by gene, exon, intergenic, etc
##'
##' Given a GRanges object, this function performs an overlap against
##' a previously created set of genomic regions. These genomic
##' regions include genes, coding portions of genes (CDS),
##' exons, intergenic regions, and exon groups (which
##' contain two or more exons)
##' @title Count RNA-Seq Pipeline Genomic Features
##' @param reads GRangesList object of interval, usually where reads aligned
##' @param features A list of genome annotations as GRangesList
##' @return A list of counts by feature
##' @author Cory Barr
##' @keywords internal
countFeatures <- function(reads, features) {
  loginfo("countGenomicFeatures.R/countFeatures: starting...")
  
  ## compute counts
  res <- lapply(features, function(features) {
    overlaps <- suppressWarnings(findOverlaps(features, reads))
    nbreads.byfeature <- tabulate(queryHits(overlaps), nbins=queryLength(overlaps))
    nbreads <- length(unique(subjectHits(overlaps)))
    names(nbreads.byfeature ) <- names(features)
    list(nbreads.byfeature=nbreads.byfeature, nbreads=nbreads)
  })

  ## rearrange output
  counts <- list()
  counts$nbreads.byfeature <- lapply(names(features), function(fname) res[[fname]]$nbreads.byfeature)
  counts$nbreads <- lapply(names(features), function(fname) res[[fname]]$nbreads)
  names(counts$nbreads.byfeature) <- names(features)
  names(counts$nbreads) <- names(features)
    
  loginfo("countGenomicFeatures.R/countFeatures: done")
  return(counts)
}

##' Compute statistics on count features
##'
##' @title Compute statistics on count features
##' @param save_dir A character string containing a NGS analysis directory
##' @param feature A character string containing a features name. Default is "counts_gene".
##' @return A numeric vector containing statistics about features.
##' @author Gregoire Pau
##' @export
##' @keywords internal
statCountFeatures <- function(save_dir, feature="counts_gene") {
  ## load data
  z <- getTabDataFromFile(save_dir, feature)$count

  ## compute some quantiles
  q <- quantile(z, c(0, 0.25, 0.5, 0.75, 1))
  names(q) <- c("min", "q025", "median", "q075", "max")
  
  return(c(nbfeatures=c(atleast.1=sum(z>=1), atleast.10=sum(z>=10), atleast.100=sum(z>=100)),
           nbreads=c(mean=mean(z), q)))
}
