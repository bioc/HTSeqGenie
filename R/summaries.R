##' Parse a summary from a list of save_dirs
##'
##' This function allows to parse a given summary from
##' a list of pipeline results save dirs
##' @title parse summary files from save dirs
##' @param save.dirs list of result dirs
##' @param summary.name name of summary file e.g. summary_counts
##' @return data frame with summaries
##' @author Jens Reeder
##' @export
##' @keywords internal
parseSummaries <- function(save.dirs, summary.name) {
  if (length(save.dirs) > 1) {
    res <- lapply(save.dirs, parseSummaries, summary.name)
    do.call(rbind, res)
  } else {
    save.dir <- save.dirs
    z <- try(getNumericVectorDataFromFile(save.dir, summary.name), silent=TRUE)
    if (class(z)=="try-error") {
      cat("parseSummaries: cannot parse summary=", summary.name, " in save.dir=", save.dir, "\n")
      z <- NULL
    }
    return(t(as.matrix(z)))
  }
}

.getAllSummaries <-function(dirs){

  ## get all summary stats into a list of data frames
  summary.types <- c("summary_preprocess", "summary_alignment",
                     "summary_variants", "summary_counts",
                     "summary_junctions", "summary_analyzed_bamstats",
                     "summary_target_lengths")

  dat.summary <- list()
  for (type in summary.types) {
    dat.summary[[type]] <- parseSummaries(dirs, type)
  }

  ## check wether all dirs have the same summaries
  if (length(unique(sapply(dat.summary, nrow))) != 1 )
    stop("inconsistent number of lanes from different summary files")

  return(dat.summary)
}

.getEffectiveReadLengths <- function(save_dirs) {  

  ## Could use processed_min_read_length processed_max_read_length now
  summary <- parseSummaries(save_dirs, "summary_preprocess")
  lengths <- subset(summary, select = c('read_length') )
  
  config_files <- file.path(save_dirs, "logs", "config.txt")
  configs <- lapply(config_files, parseDCF)
  
  .getEffectiveReadLength <- function (length, config){
    if( config$trimReads.do){
      return( min(length, as.integer(config$trimReads.length)) )
    }
    else return(length)
  }
  lengths <- mapply(.getEffectiveReadLength, lengths, configs)
  
  return(table(lengths))
}
 
.plotPreprocessSummary <- function(save_dirs, image_dir) {  
  summary <- parseSummaries(save_dirs, "summary_preprocess")

  filename <- file.path(image_dir, "summary_preprocess")  
  z <- summary[, grep(colnames(summary), pattern="read_length", value=TRUE, invert=TRUE)]
  
  return(.heatmap(z, basename(save_dirs), filename))
}

.plotAlignmentSummary <- function(save_dirs, image_dir) {

  summary <- parseSummaries(save_dirs, "summary_alignment")  
  filename <- file.path(image_dir, "summary_alignment")

  ## Maybe normalize to number of used reads per sample here?
##  z <- summary/summary[,1]
  z <- summary
  return(.heatmap(z, basename(save_dirs), filename))
}

.plotBamStatsSummary <- function(save_dirs, image_dir) {  
  summary <- parseSummaries(save_dirs, "summary_analyzed_bamstats")

  filename <- file.path(image_dir, "summary_analyzed_bamstats")
  keep <- c("nm.bin.0",
            "cigar.withS.nbrecords", "cigar.withN.nbrecords", "cigar.withID.nbrecords",
            "chromosome.bin.1", "chromosome.bin.Y", "chromosome.bin.MT")
 
  z <- round(100*summary[,keep]/summary[,'mapped.nbrecords'],2)
  
  dotboxplot(as.data.frame(z), filename=filename, ylab="Percentage of records")
  return(filename)
}

.plotVariantsSummary <- function(save_dirs, image_dir) {

  summary <- parseSummaries(save_dirs, "summary_variants")  
  filename <- file.path(image_dir, "summary_variants")

  return(.heatmap(summary, basename(save_dirs), filename))
}

.plotCountGFSummary <- function(save_dirs, image_dir) {
  summary <- parseSummaries(save_dirs, "summary_counts")
  nb_reads <- parseSummaries(save_dirs, "summary_preprocess")[,'processed_reads']

  filename <- file.path(image_dir, "summary_counts")
  counts_by_nbreads <- c('gene_nbreads',
                         'gene_exonic_nbreads',
                         'gene_coding_nbreads',
                         'transcript_nbreads',
                         'exon_nbreads',
                         'dexseq_nbreads',
                         'intergenic_nbreads',
                         'ncRNA_nbreads',
                         'ncRNA_nongenic_nbreads')
  z <- summary[,counts_by_nbreads]/nb_reads
 
  return(.heatmap(z, basename(save_dirs), filename))
}

.plotCountGFbyFeatureSummary <- function(save_dirs, image_dir) {
  summary <- parseSummaries(save_dirs, "summary_counts")
  nb_reads <- parseSummaries(save_dirs, "summary_preprocess")[,'processed_reads']

  filename <- file.path(image_dir, "summary_counts_by_feat")
  ## TODO require total number per feature
  counts_by_nbfeat <- c('counts_dexseq_nbfeatures',
                        'counts_exon_nbfeatures',
                        'counts_gene_nbfeatures',
                        'counts_gene_coding_nbfeatures',
                        'counts_gene_exonic_nbfeatures',
                        'counts_intergenic_nbfeatures',
                        'counts_ncRNA_nbfeatures',
                        'counts_ncRNA_nongenic_nbfeatures',
                        'counts_transcript_nbfeatures')  
  z <- summary[,counts_by_nbfeat]
 
  return(.heatmap(z, basename(save_dirs), filename))
}

.heatmap <- function (data, labels, filename){

  ## hack to get the labels fit into the plot:
  label.size <- 5 + max(unlist(lapply(labels,nchar)))/2
  .plot <- function(){    
    col <- colorRampPalette(c("white", "blue"))(255)
    ## Manually choose margin values that worked for my test example.
    ## Ideally this is done by R. Needs fix.
    heatmap(data, Colv=NA, scale='none', col=col,labRow=labels,
            cexCol=1, cexRow=1, margins = c(10,label.size))
    dev.off()
  }
  png(filename=paste(filename,"png", sep="."), width=900, height=750)
  .plot()

  pdf(file=paste(filename,"pdf",sep="."))
  .plot()
  return(filename)
}


.plotReadMappingsPerChr <- function(dirs, image_dir){
  ## read summaries of analyzed bamstats
  summary_analyzed_bamstats <- lapply(dirs, function(dir) {
    res <- getNumericVectorDataFromFile(dir, object_name="summary_analyzed_bamstats")
    res <- res[grepl("^chromosome.bin", names(res))]
    names(res) <- gsub("^chromosome.bin\\.", "", names(res))
    res
  })
  
  data <- as.data.frame(do.call(rbind, summary_analyzed_bamstats))
  data <- data[,grep("GL", colnames(data), invert=TRUE)] ## get rid of human contigs not assembled into chr
  data <- data[,grep("NT", colnames(data), invert=TRUE)]  ## and the same for mouse
  data <- data/1e6
  
  filename <- file.path(image_dir,"chr_mappings")
  dotboxplot(data, filename = filename,
             ylab="Number of mapped reads [million]")

  return(filename)
}
  
##' Write html Summary for list of runs
##'
##' @title Write HTML summary
##' @param dirs List of pipeline result dirs 
##' @param outdir Path to output directory. Does not create dir.
##' @param cutoffs list, cutoffs for each plotting/QA function
##' @return Nothing, but writes file
##' @author Jens Reeder
##' @export
##' @keywords internal
writeSummary <- function(dirs, outdir="Summary", cutoffs) {

  if(!file.exists(outdir)) stop("Output dir does not exist", outdir)
  image_dir <- file.path('images')  
  makeDir(image_dir, overwrite='overwrite')
  
  ## create plots
  prep.sum      <- .plotPreprocessSummary(dirs, image_dir)
  align.sum     <- .plotAlignmentSummary(dirs, image_dir)
  bamstats.sum  <- .plotBamStatsSummary(dirs, image_dir)
  variants.sum  <- .plotVariantsSummary(dirs, image_dir)
  map.per.chr   <- .plotReadMappingsPerChr(dirs, image_dir)
  qa            <- .createQAReport(dirs,cutoffs=cutoffs, image_dir)

  qa$calls <- .finalCallFromIndividualQACalls(qa$calls,cutoffs)
  write.csv(file=file.path(outdir,"qa.calls.csv"),qa$calls)
  write.csv(file=file.path(outdir,"qa.data.csv"),qa$qc.data)

  ## shortcut to write pdf/png images
  .hwriteImage <- function(basename){
    hwriteImage(paste(basename,'png', sep="."), p, br=TRUE,
                link=paste( basename, 'pdf', sep="."))
  }

  ## create html page
  outfile <- file.path("summary.html")
  p <- openPage(outfile)
  hwrite("Summary of NGS pipeline run summaries", p, heading=1)
    
  hwrite("This page summarize per-run statistics and should be used as a first overview before analyzing the results.",p,br=TRUE)
  hwrite('',p, br=TRUE)
  hwrite("This page is structured into three sections:",p,br=TRUE)
  hwrite("QA", link="#qa",p, br=TRUE)
  hwrite("Results", link="#results",p,br=TRUE)
  hwrite("Heatmaps of raw summaries", link="#heatmaps",p,br=TRUE)
  hwrite('',p, br=TRUE)
  
  hwrite('', p, name='qa')
  hwrite('Quality control', p, heading=3)
  if(length(qa$fails)>0){
    hwrite('Runs flagged for QA issues', p, heading=4)
    hwrite("These runs fail one or more of our QA criteria. Please inspect before using.",p,br=TRUE)
    hwrite(qa$fails, row.names=FALSE, p)
  }

  hwrite('Read lengths', p, heading=4)
  hwrite("Ideally all runs share the same read length. If not, it is up to the analyst to make sure the data can be combined for downstream analysis.",p,br=TRUE)
  hwrite(as.data.frame(.getEffectiveReadLengths(dirs)), p, br=TRUE)

  hwrite('QA plots', p, heading=4)
  Map(.hwriteImage, qa$files)
  
  hwrite('Reads mapped per chromosome', p, heading=4)
  hwrite('This plot shows the number of reads mapped to each chromosome. Note that forward and reverse reads are counted independently and some mappings span two chromosomes.', p, br=TRUE)
  .hwriteImage(map.per.chr)

  ##Nest
  hwrite('', name="results",p)
  hwrite('Individual pipeline run results',p, heading=3)  

  df <- data.frame("pipeline result dir"=basename(dirs),
                   "Alignment.reports"="Alignments",
                   "Feature.counting.reports" ="Feature counting",
                   "Variant.calling.reports"="Variants"
                   )

  #html sits in Summary/, but paths are relative to ./
  linksAln <- file.path(dirs, "reports/reportAlignment.html")
  linksFC  <- file.path(dirs, "reports/reportFeatureCounts.html")
  linksVar <- file.path(dirs, "reports/reportVariants.html")  
  
  hwrite(df, col.link=list(
               "Alignment.reports"=linksAln,
               "Feature.counting.reports"=linksFC,
               "Variant.calling.reports"=linksVar), p)
  hwrite('',p, br=TRUE)

  hwrite("", p, name="heatmaps")
  hwrite("Heatmap of raw summaries", p, heading=3, br=TRUE)
  
  hwrite(hwrite("Summary of preprocessing and read filtering",  heading=4), p, name="heatmaps", br=TRUE)
  hwrite("This heatmap shows the raw number of reads without any normalization or scaling.", p, br=TRUE)
  .hwriteImage(prep.sum)
  hwrite('',p, br=TRUE)

  hwrite("Summary of Gsnap alignments",p, heading=4)
  hwrite("This heatmap shows the number of mapped reads that fall into each gsnap category.", p, br=TRUE)
  .hwriteImage(align.sum)
  hwrite('',p, br=TRUE)

  hwrite("Summary of analyzed.bam statistics",p, heading=4)
  hwrite("This heatmap shows the a selected subset of the analyzed.bam statistics.", p, br=TRUE)
  .hwriteImage(bamstats.sum)
  hwrite('',p, br=TRUE)

  countGF.sum       <- .plotCountGFSummary(dirs, image_dir)
  countGFbyFeat.sum <- .plotCountGFbyFeatureSummary(dirs, image_dir)
  
  hwrite("Summary of feature countings", p, heading=4)
  hwrite("This heatmap shows the number of mapped reads per genomic feature category.", p, br=TRUE)  
  hwrite("Each row/sample is normalized by the number of analyzed reads of that particular sample.", p, br=TRUE)   
  .hwriteImage(countGF.sum)
  hwrite('', p, br=TRUE)
    
  hwrite("This heatmap shows the number of genomic features hit by at least one read per feature category.", p, br=TRUE)  
  hwrite("Ideally, this would be normalized by the number of features in each category, but that's not yet implemented.", p, br=TRUE)   
  .hwriteImage(countGFbyFeat.sum)
  hwrite('', p, br=TRUE)
  
  hwrite("Summary of variant calling", p, heading=4)
  hwrite("This heatmap shows some statistics about the called variants.", p, br=TRUE)  
  .hwriteImage(variants.sum)
  hwrite('',p, br=TRUE)

  close(p)  
}

.createQAReport <- function(dirs,cutoffs,image_dir) {
  plot.names = c("good_reads","perc_good_reads","bad_reads","perc_bad_reads",
    "detected_genes","junction_reads","in_out_min_max_readlen", "target_quartiles_pct")
  file_list <- file.path(image_dir,plot.names)

  ### Switching to strategy where plot functions return one or more columns of PASS/FAIL
  ###  and these are collected in one data.frame
  dat.summary <- .getAllSummaries(dirs)

  # DF to hold result of each QA test
  qc.data = data.frame(Technology=rep(cutoffs$seqType,length(dirs)),row.names = basename(dirs))

  ## get all reads count related numbers into one matrix
  dat.reads <-
    cbind(dat.summary$summary_preprocess,
          dat.summary$summary_alignment,
          dat.summary$summary_counts[,grepl("nbreads",
                                            colnames(dat.summary$summary_counts))])

  ## combine by mapping category
  dat.reads <-
    cbind(dat.reads,
          mult_mapping=rowSums(dat.reads[,grepl("_mult",colnames(dat.reads))]),
          unpaired_mapping=rowSums(dat.reads[,grepl("unpaired_",colnames(dat.reads))]),
          half_mapping=rowSums(dat.reads[,grepl("halfmapping_",colnames(dat.reads))]))

  ## Provide normalization by total, processed and analyzed number of reads
  perc.total <- dat.reads / dat.reads[,"total_reads"] * 100
  colnames(perc.total) <- paste(colnames(perc.total),"/total",sep="")

  perc.processed <- dat.reads / dat.reads[,"processed_reads"] * 100
  colnames(perc.processed) <- paste(colnames(perc.processed),"/processed",sep="")
  
  perc.analyzed <- dat.reads / dat.reads[,"analyzed"] * 100
  colnames(perc.analyzed) <- paste(colnames(perc.analyzed),"/analyzed",sep="")

  target.length.bin.cols <- colnames(dat.summary$summary_analyzed_bamstats)[ grepl("^target_length.bin",colnames(dat.summary$summary_analyzed_bamstats)) ]
  target.length.bin.591 = (dat.summary$summary_analyzed_bamstats[,"target_length.bin.591"] / rowSums( dat.summary$summary_analyzed_bamstats[ , target.length.bin.cols ])) * 100

  column.nm <- colnames(dat.summary$summary_analyzed_bamstats)[ grepl("^nm.bin",colnames(dat.summary$summary_analyzed_bamstats)) ]
  total.nm = rowSums(dat.summary$summary_analyzed_bamstats[,column.nm])
  ok.nm = rowSums(dat.summary$summary_analyzed_bamstats[,c("nm.bin.0","nm.bin.1")])
  nm2plus = ((total.nm - ok.nm) / total.nm) * 100

  perc.reads <- cbind(perc.total, perc.processed, perc.analyzed,target.length.bin.591, nm2plus, check.names=FALSE)

  ## doing the plots
  ## number/perc of good reads

  qc.data <- cbind(qc.data,.plot.reads.good(dat.reads, image_dir,cutoffs))
  qc.data <- cbind(qc.data,.plot.perc.reads.good(perc.reads, image_dir,cutoffs))

  ## number/perc of bad reads
  qc.data <- cbind(qc.data,.plot.reads.bad(dat.reads, image_dir))
  qc.data <- cbind(qc.data,.plot.perc.reads.bad(perc.reads, image_dir))

  qc.data$genesDetected <- .plotGenesDetected(dat.summary$summary_counts, image_dir)
  qc.data[,"junctionReads"] <- .plotJunctionReads(dat.summary$summary_junctions,
                              dat.summary$summary_analyzed_bamstats, image_dir,cutoffs)

  qc.data[,"inOutLenEqual"]  <- .plot.in.out.min.max.readlength(dat.summary$summary_preprocess, image_dir, cutoffs)

  qc.data <- cbind(qc.data,.plot.quartile.target.length(dat.summary$summary_target_lengths, image_dir))

  qc.data[,"withS_pct"] = round(100 * (dat.summary$summary_analyzed_bamstats[,"cigar.withS.nbrecords"] / dat.summary$summary_analyzed_bamstats[,"mapped.nbrecords"]),2) # Plotted elsewhere
  qc.data[,"MT_pct"]    = round(100 * dat.summary$summary_analyzed_bamstats[,"chromosome.bin.MT"] / dat.summary$summary_analyzed_bamstats[,"mapped.nbrecords"],2) # Plotted elsewhere

  calls <- .reportFailedRuns(dirs, qc.data, cutoffs)
  
  return(list(files=file_list, qc.data=qc.data, calls=calls))
}

# Given calls from .createQAReport, decide which samples are OK
.finalCallFromIndividualQACalls <- function(calls,cutoffs) {

  # Eventually be smart by weighting certain QA tests higher
  if ( cutoffs$seqType %in% "RNASeq" ) {
    fail.cols = c("inOutLenEqual","analyzed","geneReadsPerc","junctionReads")
  } else {
    fail.cols = c("inOutLenEqual","analyzed","geneReadsPerc","junctionReads","cigarWithS","mtReads","targetBin591")
  }

  warn.cols = setdiff(colnames(calls),c("Technology",fail.cols))
  fail.calls = calls[, fail.cols]
  warn.calls = calls[, warn.cols]
  failed.tests = apply( fail.calls, 1, function(x) { paste(colnames(fail.calls)[ x ], collapse=",") } )
  warning.tests = apply( warn.calls, 1, function(x) { paste(colnames(warn.calls)[ x ], collapse=",") } )
  calls$Final = rowSums(fail.calls) > 0
  calls$WarningTests = warning.tests
  calls$FailedTests = failed.tests
  calls$ManualFailedTests = failed.tests
  return(calls)
}

.plot.in.out.min.max.readlength <- function(preprocess, image_dir,cutoffs) {
  filename <- file.path(image_dir,"in_out_min_max_readlen")
  plotInOut <- function(preprocess) {
    par(mfrow=c(2,2))
    barplot(prop.table(table(preprocess[,"input_min_read_length"]))*100,main="Input Min Read Length",ylab="Percentage with Length")
    barplot(prop.table(table(preprocess[,"input_max_read_length"]))*100,main="Input Max Read Length",ylab="Percentage with Length")
    barplot(prop.table(table(preprocess[,"processed_min_read_length"]))*100,main="Output Min Read Length",ylab="Percentage with Length")
    barplot(prop.table(table(preprocess[,"processed_max_read_length"]))*100,main="Output Max Read Length",ylab="Percentage with Length")
    return(invisible(NULL))
  }
  png(filename=paste(filename,".png",sep=""))
  plotInOut(preprocess)
  dev.off()
  pdf(file=paste(filename,".pdf",sep=""))
  plotInOut(preprocess)
  dev.off()
  max.len.mode = median(preprocess[,"processed_max_read_length"])
  min.len.mode = median(preprocess[,"processed_min_read_length"])
  qc = (preprocess[,"processed_min_read_length"] == preprocess[,"processed_max_read_length"]) & (preprocess[,"processed_min_read_length"] == min.len.mode)
  return(qc)
}

.plot.quartile.target.length <- function(dat.reads, image_dir) {
  filename <- file.path(image_dir,"target_quartiles_pct")
  qc = as.data.frame(dat.reads[,c(2,3,5)])
  colnames(qc) = c("TargetLength25thPercentile","TargetLengthMedian","TargetLength75thPercentile")
  dotboxplot(qc,filename=filename,ylab="Target Length")
  return(qc)
}

.plot.reads.good <- function(dat.reads, image_dir, cutoffs){
  filename <- file.path(image_dir,"good_reads")
  column.reads.good <- c("total_reads","processed_reads","analyzed",
                         "concordant_uniq",
                         "gene_nbreads","gene_exonic_nbreads")
  qc = data.frame( dat.reads[,column.reads.good], check.names=FALSE )
  dotboxplot(qc/1e6,
             filename=filename, cutoff=10,
             ylab="Number of reads [million]"
  )
  return(qc)
}

.plot.perc.reads.good <- function(perc.reads, image_dir, cutoffs){
  filename <- file.path(image_dir,"perc_good_reads")
  column.preads.good <- c("processed_reads/total","highqual_reads/total","analyzed/total",
                          "concordant_uniq/processed",
                          "gene_nbreads/analyzed","gene_exonic_nbreads/analyzed")
  
  qc = data.frame(perc.reads[,column.preads.good],check.names=FALSE)
  dotboxplot(qc, filename=filename, cutoff=cutoffs[["geneReadsPerc"]], ylab="Percentage of reads",ylim=c(0,100))
  return(qc)
}

.plot.reads.bad <- function(dat.reads, image_dir){
  filename <- file.path(image_dir,"bad_reads")
  column.reads.bad <- c("adapter_contam","rRNA_contam_reads",
                        "nomapping","mult_mapping",
                        "unpaired_mapping","half_mapping")
  qc = data.frame(dat.reads[,column.reads.bad],check.names=FALSE)
  dotboxplot(qc/1e6,
             filename=filename, ylab="Number of reads [million]")
  return(qc)
}

.plot.perc.reads.bad <- function(perc.reads, image_dir){
  filename <- file.path(image_dir,"perc_bad_reads")
  column.preads.bad <- c("adapter_contam/total","rRNA_contam_reads/total",
                         "nomapping/processed","mult_mapping/processed",
                         "unpaired_mapping/processed","half_mapping/processed", "target.length.bin.591","nm2plus","nomapping/total")
  qc = data.frame(perc.reads[,column.preads.bad],check.names=FALSE)
  dotboxplot(qc,
             filename=filename, ylab="Percentage of reads",ylim=c(0,100))
  return(qc)
}

.plot.nm2.pct <- function(dat.reads, image_dir) {
  filename <- file.path(image_dir,"nm2_reads")
  column.nm <- colnames(dat.reads)[ grepl("^nm.bin",colnames(dat.reads)) ]
  total.nm = rowSums(dat.reads[,column.nm])
  ok.nm = rowSums(dat.reads[,c("nm.bin.0","nm.bin.1")])
  nm2plus.pct = ((total.nm - ok.nm) / total.nm) * 100
  dotboxplot(nm2plus.pct,filename=filename,ylab="Percent of Reads With 2+ Mismatches",ylim=c(0,100))
  return(filename)
}
  
.plotGenesDetected <- function(summary_counts, image_dir){
   filename <- file.path(image_dir,"detected_genes")
   qc = summary_counts[,"counts_gene_nbfeatures"]/1000
   dotboxplot(qc,
             filename=filename, ylab="Number of detected genes [thousand]")
  return(qc)
}

.plotJunctionReads <- function(summary_junctions, summary_bamstats, image_dir, cutoffs){
  filename <- file.path(image_dir,"junction_reads")
  data <- summary_junctions[,"count_junction_reads"]
  totals <- summary_bamstats[,"mapped.nbrecords"]
  x <- 100 *data.frame(data/totals)
  dotboxplot(x,filename=filename, ylab="Percent of mapped reads that are Junction Reads",ylim=c(0,100))
  qc = x[,1,drop=FALSE]
  return(qc)
}
  

## output the list of failed lanes and their failure mode
.reportFailedRuns <- function(dirs, qc.data, cutoffs) {

  calls = data.frame(Technology=rep(cutoffs$seqType,length(dirs)),row.names = basename(dirs))

  calls[,"goodReads"] = qc.data [,"total_reads"] < cutoffs$totalReads
  calls[,"analyzed"] = qc.data[,"analyzed"] < cutoffs$processedReadsAnalyzed
  calls[,"concordantUniq"] = qc.data[,"concordant_uniq"] < cutoffs$concordantUniq
  calls[,"analyzedPerc"] = qc.data[,"analyzed/total"] < cutoffs[["analyzedPerc"]]
  calls[,"concordantUniqPerc"] = qc.data[,"concordant_uniq/processed"] < cutoffs[["concordantUniqPerc"]]
  calls[,"geneReadsPerc"] = qc.data[,"gene_nbreads/analyzed"] < cutoffs[["geneReadsPerc"]]
  calls[,"geneExonicReadsPerc"] = qc.data[,"gene_exonic_nbreads/analyzed"] < cutoffs[["exonicPerc"]]
  if (cutoffs$seqType %in% "RNASeq") {
    calls[,"junctionReads"]= qc.data[,"junctionReads"] < cutoffs$junctionReads
  } else {
    calls[,"junctionReads"]= qc.data[,"junctionReads"] > cutoffs$junctionReads
    calls[,"targetBin591"]= qc.data[,"target.length.bin.591"] > cutoffs$targetBin591
  }
  calls[,"inOutLenEqual"] = ! qc.data[,"inOutLenEqual"]
  calls[,"cigarWithS"] = qc.data[,"withS_pct"] > cutoffs$cigarWithS
  calls[,"mtReads"] = qc.data[,"MT_pct"] > cutoffs$mtReads
  calls[,"nonMappingPct"] = qc.data[,"nomapping/processed"] > cutoffs$nonMappingPct
  calls[,"highQualPct"] = qc.data[,"highqual_reads/total"] < cutoffs$highQualPct
  return(calls)
}

dotboxplot <- function(dat, filename, cutoff=NULL, ylim=NULL, dot.col="firebrick",
                       dot.pch=19, las=2, cex=0.4, ...) {

  if ( is.null(ylim) ) {
    if ( class(dat) == "formula" ) {
      ylim <- range(eval(dat[[2]]))
    } else {
      ylim <- range(unlist(dat))
    }
  }

  .plot <- function(){
    ## Add a big margin at below x for long labels
    boxplot(dat,outline=FALSE ,ylim=ylim,
            par=par(mar=c(15, 4, 4, 2) +  0.1), las=las, cex=cex, ...)
    stripchart(dat,vertical=TRUE,method="jitter",pch=dot.pch,add=TRUE, col=dot.col, ...)
 
    for (i in seq(0,100,10)) {
      abline(h=i,lty=3,col='grey')
    }
    if(! is.null(cutoff)){
      abline(h=cutoff,lty=2,col='blue')
    }
    dev.off()
  }
 
  pdf(file=paste(filename, "pdf", sep="."))
  .plot()
  
  png(filename=paste(filename, "png", sep="."),width=600, height=800)
  .plot()
}
