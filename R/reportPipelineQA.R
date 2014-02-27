##' Generates a summary HTML for the preprocess and align step
##'
##' @title Generate Pipeline Report
##' @return Name of created HTML file
##' @author Melanie Huntley, Cory Barr, Jens Reeder
##' @export
##' @importFrom hwriter hwrite
##' @importFrom Cairo CairoPDF
##' @keywords internal
writePreprocessAlignReport <- function() {
  loginfo("reportPipelineQA.R/writePreprocessAlignReport: creating joint preprocess-alignment report...")
  pdf <- CairoPDF ## to avoid random crashes caused by parallel pdf writing using the pdf() device!

  ## setup output dirs
  dirPath <- getConfig("save_dir")
  report_dir <- file.path(dirPath, 'reports')
  image_dir <- file.path(report_dir, "images")
  paired_ends <- getConfig.logical("paired_ends")
  
  ## load summary data
  summary.preprocessTable <- getNumericVectorDataFromFile(dirPath, "summary_preprocess")
  summary.alignTable <- getNumericVectorDataFromFile(dirPath, "summary_alignment")

  reads_after_preprocess <- summary.preprocessTable[["processed_reads"]]
  ## compute filtering and align stats
  readFiltering <- .makeReadFiltering(summary.preprocessTable, summary.alignTable)
  readFilteringTable <- data.frame(Read_counts=readFiltering,
                                   Percent_of_total=sprintf("%.2f", readFiltering/readFiltering[1]*100))
  row.names(readFilteringTable) <- names(readFiltering) 

  gsnap_mapped_reads_files <- names(summary.alignTable)
  readMappings <- summary.alignTable[gsnap_mapped_reads_files]
  ReadMappingsTable <- .makePreprocessAlignPlots(readFiltering, readMappings, reads_after_preprocess,
                                                image_dir, summary.alignTable)

  sanity_check <- sanityCheck(summary.alignTable, readFiltering, paired_ends)
  if(!all(sanity_check$Status == "PASS")){
    logwarn("**********************")
    logwarn("Sanity checks failed!!")
    logwarn("**********************")
  }

  if(paired_ends) {
    concord_uniq <- getBams(dirPath)$"concordant_uniq"
    targetLengthTable <- calculateTargetLengths(concord_uniq, dirPath)
  }

  outfile=file.path(report_dir, "reportAlignment.html")
  writePreprocessAlignHTML(outfile, dirPath, sanity_check,
                           readFilteringTable,
                           ReadMappingsTable,
                           targetLengthTable)
  loginfo("reportPipelineQA.R/writePreprocessAlignReport: done")
  return(outfile)
}

##' Generates a summary HTML for the Genomic Feature counting step
##'
##' @title Generate pipeline report
##' @return Name of created HTML file
##' @author Melanie Huntley, Cory Barr, Jens Reeder
##' @export
##' @keywords internal
writeGenomicFeaturesReport <- function(){
  loginfo("reportPipelineQA.R/writeGenomicFeaturesReport: creating report of genomic features counts...")
  pdf <- CairoPDF ## to avoid random crashes caused by parallel pdf writing using the pdf() device!

  ## setup output dirs
  dirPath <- getConfig("save_dir")
  report_dir <- file.path(dirPath, 'reports')
  image_dir <- file.path(report_dir, "images")

  ## load summary data
  summary.alignTable <- getNumericVectorDataFromFile(dirPath, "summary_alignment")
  summary.countsTable <- getNumericVectorDataFromFile(dirPath, "summary_counts")
  
  GenomicFeaturesTable <- .calcHitsByGenomicFeature(dirPath, image_dir, summary.countsTable,
                                                   summary.alignTable[['analyzed']])
  GenomicFeaturesDetectedTable <- .calcGenomicFeaturesDetected(dirPath, image_dir)
  ExonsCoveredTable <- .calcExonsCovered(dirPath)

  .generateCountFeaturesPlots(dirPath, image_dir)

  outfile=file.path(report_dir, "reportFeatureCounts.html")
  writeFeatureCountsHTML(outfile, dirPath,
                         ExonsCoveredTable, GenomicFeaturesTable,
                         GenomicFeaturesDetectedTable)
  loginfo("reportPipelineQA.R/writeGenomicFeaturesReport: done")
  return(outfile)
}

.makeReadFiltering <- function(summary.preprocessTable,
                              summary.alignTable){
  reads_after_preprocess<- summary.preprocessTable[["processed_reads"]]
  mapped_reads <- reads_after_preprocess - summary.alignTable[["nomapping"]]
  
  unique_mapping   <- sum(summary.alignTable[grep("uniq", names(summary.alignTable))])
  multi_mapping    <- sum(summary.alignTable[grep("mult", names(summary.alignTable))])
  transloc_mapping <- sum(summary.alignTable[grep("transloc", names(summary.alignTable))])
  readFiltering <- c("Total"           = summary.preprocessTable[["total_reads"]],
                     "High-Qual"       = summary.preprocessTable[["highqual_reads"]],
                     "Adapter contam"  = summary.preprocessTable[["adapter_contam"]],
                     "rRNA reads"      = summary.preprocessTable[["rRNA_contam_reads"]],
                     "Post-preprocess" = reads_after_preprocess,
                     "Mapped"          = mapped_reads,
                     "Uniquely"        = unique_mapping,
                     "Multiple"        = multi_mapping,
                     "Translocation"   = transloc_mapping,
                     "Unmapped"        = summary.alignTable[["nomapping"]],
                     "Analyzed"        = summary.alignTable[["analyzed"]]
                     )
  return(readFiltering)
}

.makePreprocessAlignPlots <- function(readFiltering, readMappings, total_gsnap_processed,
                                     image_dir, summaryTable){
#####################################
  ## Plot barchart of reads and
  ## filtering steps
#####################################

  relativeBarPlot(data=readFiltering,
                  total=readFiltering[1], ## total_reads
                  labels=names(data),
                  title="Summary of read filtering and mapping",
                  ylab = "Percent of reads",
                  filename=file.path(image_dir, "ReadSummaryFilteringMapping"))
  
  ####################################
  ## Plot mapping types
  ####################################
  ## Make names look nicer in plot
  names(readMappings) <- gsub("_", " ", names(readMappings))
  readMappings <- readMappings[grep("analyzed", names(readMappings), invert=TRUE)]
  labels <- sub(" ", "\n", names(readMappings))

  ymax <- round(max(readMappings/total_gsnap_processed*100))

  relativeBarPlot(data=readMappings,
           total=total_gsnap_processed,
           title="Summary of read mappings",
           ylab="Percent of mapped reads", 
           labels=labels,
           cex.names=0.8,
           ymax=ymax,
           filename=file.path(image_dir, "ReadSummaryMappings"))
  
  ReadMappingsTable <- data.frame(Read_counts=readMappings,
                                  Percent_of_total=sprintf("%.2f", (readMappings/total_gsnap_processed)*100))
  row.names(ReadMappingsTable) <- names(readMappings)
  return(ReadMappingsTable)
}

.calcHitsByGenomicFeature <- function(dirPath, image_dir, summary.Table, nbreads){
  ####################################
  ## Reads by chromosome
  ####################################
  ## genic reads by chromosome
  hits_by_gene <- getTabDataFromFile(dirPath, "counts_gene")$count
  
  ####################################
  ## types of genomic features mapped to
  ####################################
  hits_by_gene_exonic <- getTabDataFromFile(dirPath, "counts_gene_exonic")$count

  hits_by_gene_coding <- getTabDataFromFile(dirPath, "counts_gene_coding")$count

  intron_reads <- summary.Table[["gene_nbreads"]] - summary.Table[["exon_nbreads"]]

  hits_by_genomicFeature <- c("Genes" = summary.Table[["gene_nbreads"]],
                              "Exons" = summary.Table[["exon_nbreads"]],
                              "CDS" = summary.Table[["gene_coding_nbreads"]],
                              "Introns" = intron_reads,
                              "ncRNA" = summary.Table["ncRNA_nbreads"],
                              "ncRNA_nongenic" = summary.Table["ncRNA_nongenic_nbreads"],
                              "Intergenic" = summary.Table[["intergenic_nbreads"]])
  
  ymax <- floor(max((hits_by_genomicFeature/nbreads)*100)+0.5)
  
  relativeBarPlot(data=hits_by_genomicFeature,
                  total=nbreads,
                  labels=names(hits_by_genomicFeature),
                  ymax=ymax,
                  ylab="Percent of mapped reads",
                  title="Summary of genomic features covered by read mappings",
                  filename = file.path(image_dir, "ReadSummaryGenomicFeaturesMapping"))
  
  GenomicFeaturesTable <- data.frame(Read_counts=hits_by_genomicFeature,
                                     Percent_of_total=sprintf("%.2f", hits_by_genomicFeature / nbreads * 100))
  row.names(GenomicFeaturesTable) <- names(hits_by_genomicFeature)

  return(GenomicFeaturesTable)
}

.calcGenomicFeaturesDetected <- function(dirPath, image_dir){

  ## genic reads by chromosome
  hits_by_gene <- getTabDataFromFile(dirPath, "counts_gene")$count 
  hits_by_gene_exonic <- getTabDataFromFile(dirPath, "counts_gene_exonic")$count
  hits_by_gene_coding <- getTabDataFromFile(dirPath, "counts_gene_coding")$count
  hits_by_exon <- getTabDataFromFile(dirPath, "counts_exon")$count
  hits_by_transcript <- getTabDataFromFile(dirPath, "counts_transcript")$count
  hits_by_ncRNA <- tryCatch(getTabDataFromFile(dirPath, "counts_ncRNA")$count,
                            error=function(e) numeric(0))
  hits_by_ncRNA_nongenic <- tryCatch(getTabDataFromFile(dirPath, "counts_ncRNA_nongenic")$count,
                                     error=function(e) numeric(0))
  hits_by_intergenic <- getTabDataFromFile(dirPath, "counts_intergenic")$count

  ####################################
  ## number of genomic features
  ## detected
  ####################################
  genes_detected <- sum(hits_by_gene > 0)
  total_genes <- length(hits_by_gene)
  
  genes_exonic_detected <- sum(hits_by_gene_exonic > 0)
  total_genes_exonic <- length(hits_by_gene_exonic)
  
  genes_coding_detected <- sum(hits_by_gene_coding > 0)
  total_genes_coding <- length(hits_by_gene_coding)

  exons_detected <- sum(hits_by_exon > 0)
  total_exons <- length(hits_by_exon)

  transcripts_detected <- sum(hits_by_transcript > 0)
  total_transcripts <- length(hits_by_transcript)

  ncRNA_detected <- sum(hits_by_ncRNA > 0)
  total_ncRNA <- length(hits_by_ncRNA)

  ncRNA_nongenic_detected <- sum(hits_by_ncRNA_nongenic > 0)
  total_ncRNA_nongenic <- length(hits_by_ncRNA_nongenic)
  
  intergenic_detected <- sum(hits_by_intergenic > 0)
  total_intergenic <- length(hits_by_intergenic)

  genomicFeatures_detected <- c("Genes" = genes_detected,
                                "Genes (exonic)" = genes_exonic_detected,
                                "CDS" = genes_coding_detected,
                                "Exons" = exons_detected,
                                "Transcripts" = transcripts_detected,
                                "ncRNA" = ncRNA_detected,
                                "ncRNA_nongenic" = ncRNA_nongenic_detected,
                                "Intergenic" = intergenic_detected)
 
  total <- c(total_genes, total_genes_exonic, total_genes_coding, total_exons,
             total_transcripts, total_ncRNA, total_ncRNA_nongenic, total_intergenic)
  ymax <- floor(max((genomicFeatures_detected/total)*100)+0.5)                     

  relativeBarPlot(data = genomicFeatures_detected,
                  total = total,
                  labels=names(genomicFeatures_detected),
                  ymax=ymax,
                  ylab="Percent of genomic feature category detected",
                  title="Summary of genomic features detected",
                  filename=file.path(image_dir, "GenomicFeaturesDetected"))
  
  GenomicFeaturesDetectedTable <- data.frame(Features_detected=genomicFeatures_detected,
                                             Percent=sprintf("%.2f",(genomicFeatures_detected/total)*100))
  row.names(GenomicFeaturesDetectedTable) <- names(genomicFeatures_detected)
  return(GenomicFeaturesDetectedTable)
}

.calcExonsCovered <- function(dirPath, min.coverage=10){  
  ### Compute exons with coverage at least min.coverage
  hits_by_exon <- getTabDataFromFile(dirPath, "counts_exon")$count
  total_exons <- length(hits_by_exon)
  
  exons_covered <- sum(hits_by_exon > min.coverage)
  ExonsCoveredTable <- data.frame(ExonsCovered=exons_covered,
                                  PercentOfExonsCovered=sprintf("%.2f",(exons_covered/total_exons)*100))
  row.names(ExonsCoveredTable) <- "At_least_10_reads"
  return(ExonsCoveredTable)
}

.generateCountFeaturesPlots <- function(dirPath, image_dir){
  ## Since we have to redo almost the same plot multiple times, we
  ## define this local helper function that will use variables defined in this function scope
  .plotDFhelper <- function(counts, limit, filename){
    data <- counts[counts$count <= limit,]
    fp <- file.path(image_dir, filename)
    # there is nothing to plot when all data is filtered
    if(length(data$count) > 0){
      plotDF(data,
             xlab=xlab,
             ylab=ylab,
             filename=fp)
      invisible(fp)
    }
  }
  
  hits_by_gene_exonic <- getTabDataFromFile(dirPath, "counts_gene_exonic")$count

  ## genes with AT LEAST plots
  df <- how_many(hits_by_gene_exonic, decreasing=TRUE)
  ylab <- "Percent of genes (exonic) with at least x reads"
  xlab <- "Reads mapped to a gene"  

  .plotDFhelper(df, limit=max(df$count), filename="DistrOfReadsToGenesExonic")
  .plotDFhelper(df, limit=2000, filename="DistrOfReadsToGenesExonic2000")
  .plotDFhelper(df, limit=200, filename="DistrOfReadsToGenesExonic200")
  
  ## RPKMs
  rpkm_by_gene_exonic <- getTabDataFromFile(dirPath, "gene_exonic")$rpkm
  df <- how_many(rpkm_by_gene_exonic, decreasing=TRUE)
  ## if there are no mappings, rpkm is undefined
  if(all(is.finite(df$count))){
    xlab <- "RPKM"
    ylab="Percent of genes (exonic) with an RPKM of at least x"
    .plotDFhelper(df, limit=max(df$count), filename="DistrOfRPKMsToGenesExonic")    
    .plotDFhelper(df, limit=200, filename="DistrOfRPKMsToGenesExonic200")
    .plotDFhelper(df, limit=50, filename="DistrOfRPKMsToGenesExonic50")
  }
}

##' writePreprocessAlignHTML
##'
##' @title writePreprocessAlignHTML
##' @param outfile a path
##' @param dirPath a path
##' @param sanity_check a logical
##' @param readFilteringTable a table
##' @param ReadMappingsTable a table
##' @param targetLengthTable a table
##' @return Nothing
##' @author Gregoire Pau
##' @keywords internal
##' @importFrom hwriter openPage hwrite hwriteImage
writePreprocessAlignHTML <- function(outfile, dirPath, sanity_check,
                                      readFilteringTable, ReadMappingsTable,
                                      targetLengthTable){
  p <- openPage(outfile)
  hwrite('Summary and QA of Pipeline Output', p, heading=1)
  hwrite(paste('Pipeline output directory:',dirPath),p,br=TRUE)

  hwrite('Config file: ', p)
  hwrite('config.txt', p, link='../logs/config.txt', br=TRUE)
  
  hwrite('Valid input and output files', p, heading=2)
  hwrite('This section checks that the number of reads  sent to gSNAP for alignment are
          accounted for in the gSNAP output files. Additionally we check that there are more unique mappings than multi mappings.',p,br=TRUE)
  hwrite('',p,br=TRUE)
  hwrite(sanity_check,p)
  hwrite('',p,br=TRUE)

  hwrite('ShortRead QA summary reports for filtered reads',p,heading=2)
  hwrite('Full report for paired end 1', p,
         link=file.path("shortReadReport_1", "index.html"), br=TRUE)
  if(getConfig.logical("paired_ends")) {
    hwrite('Full report for paired end 2',p,
           link=paste("shortReadReport_2", "index.html",sep="/"),br=TRUE)
  }

  hwrite('Read counts and filtering',p,heading=2)
  hwrite('Here we tabulate the number of reads at each step of preprocessing and alignment. We begin with the total
number of reads reported by the sequencing machine software (this often includes a filtering step by the native
software). We then apply our own quality filtering step and remove reads that fail our quality filtering threshold.
In order to pass our quality filtering threshold a read must have more than 70% of the nucleotides with a quality
score of 23 or higher.
Reads are then mapped to a curated database of ribosomal RNA sequences to detect and filter ribosomal reads from
the dataset. We then estimate adapter contamination by aligning the reads to a curated set of adapter sequences.
Adapter contamined reads are not removed from the dataset, however, as portions of those reads may still align well
to the genome. The reads are then mapped to the genome and we report the number of reads that mapped (broken down by
uniquely mapping or multply mapping) and the number of reads for which we find no mapping.',p,br=TRUE)
  hwrite('',p,br=TRUE)
  hwrite(readFilteringTable, p, style='text-align:right')
  hwrite('',p, br=TRUE)
  hwriteImage(file.path('images','ReadSummaryFilteringMapping.png'), p, br=TRUE,
              link=file.path("images","ReadSummaryFilteringMapping.pdf"))
  hwrite('',p, br=TRUE)
        
  hwrite('Read mappings', p, heading=2)
  hwrite('Paired end read mappings can fall into several mapping categories depending on the number of times each
end finds a suitable mapping position in the genome, and the distance and orientation of the two ends with
respect to each other. When the two ends are antisense to each other and within 200,000bp (but not overlapping)
they are considered concordant. Otherwise they are considered unpaired. If only one end maps the paired end reads
are considered to be halfmapping. When there is a single mapping to the genome, the reads are considered to be
uniquely  mapped, whereas multiple mappings are called multi-mappers.', p, br=TRUE)
  hwrite('', p, br=TRUE)
  hwrite(ReadMappingsTable, p, style='text-align:right')
  hwrite('', p, br=TRUE)
  hwriteImage(file.path('images', 'ReadSummaryMappings.png'), p, br=TRUE, link=file.path("images", "ReadSummaryMappings.pdf"))
  hwrite('', p, br=TRUE)    

  if(getConfig.logical("paired_ends")) {  
    hwrite('Distribution of target lengths for mapped reads',p,heading=2)
    hwrite(targetLengthTable, p, style='text-align:right')
    hwrite('',p,br=TRUE)
    hwriteImage(paste('images','TargetLengths.png',sep="/"),p,br=TRUE,link=paste("images","TargetLengths.pdf",sep="/"))
    hwrite('',p,br=TRUE)
  }

  hwrite('Record statistics on "analyzed.bam"',p,heading=2)
  bstats <- getNumericVectorDataFromFile(dirPath, object_name="summary_analyzed_bamstats")
  dat <- c(bstats["nm.bin.0"], bstats["mapq.bin.40"],
           bstats["cigar.withS.nbrecords"], bstats["cigar.withN.nbrecords"], bstats["cigar.withID.nbrecords"],
           bstats["chromosome.bin.1"], bstats["chromosome.bin.Y"], bstats["chromosome.bin.MT"])
  labels <- c("with 0 mismatch", "with MAPQ 40",
              "with soft clipping [S]", "with splicing [N]", "with indel [ID]",
              "on chromosome 1", "on chromosome Y", "on mitochondrion")
  names(dat) <- labels
  
  df <- data.frame("Number of records"=c(bstats["mapped.nbrecords"], dat))
  df$"Percentage" <- round(100*df[,1]/bstats["mapped.nbrecords"], 2)
  rownames(df) <- c("total", labels)
  hwrite(df, p, br=TRUE)
 
  filename <-  file.path(dirPath, 'reports', "images", "analyzed_bamstats")
  relativeBarPlot(dat, bstats["mapped.nbrecords"], labels,
                  title="Summary of record statistics on 'analyzed.bam'",
                  filename=filename, ylab="Percentage of records")
  hwriteImage(file.path('images', 'analyzed_bamstats.png'), p, br=TRUE, link=file.path("images", "analyzed_bamstats.pdf"))
  
  close(p)
}

##' writeFeatureCountsHTML
##'
##' @title writeFeatureCountsHTML
##' @param outfile a path
##' @param dirPath  a path
##' @param ExonsCoveredTable a table
##' @param GenomicFeaturesTable a table
##' @param GenomicFeaturesDetectedTable a table
##' @return Nothing
##' @export
##' @keywords internal
##' @author Gregoire Pau
##' @importFrom hwriter openPage hwrite hwriteImage closePage
writeFeatureCountsHTML <- function(outfile, dirPath,
                                   ExonsCoveredTable, GenomicFeaturesTable,
                                   GenomicFeaturesDetectedTable){
  p <- openPage(outfile)
  hwrite('Read mappings by genomic feature', p, heading=2)
  hwrite('Here we tabulate the number of uniquely mapped reads that map to various categories of genomic features: coding genes (including exons and introns), genic exons (including coding exons and UTRs), coding exons (excluding UTRs),introns, non-coding RNAs, non-coding RNAs having no overlap with coding genes, and intergenic regions.',p,br=TRUE)
  hwrite('', p, br=TRUE)
  hwrite(GenomicFeaturesTable, p, style='text-align:right')
  hwrite('', p, br=TRUE)
  hwriteImage(file.path('images', 'ReadSummaryGenomicFeaturesMapping.png'),
              p,
              br=TRUE,
              link=file.path("images", "ReadSummaryGenomicFeaturesMapping.pdf"))
  hwrite('', p, br=TRUE)

  hwrite('Exons covered by at least 10 reads',p,heading=2)
  hwrite('Here we tabulate the total number of exons for which at least 10 reads map uniquely. We also calculate this quantity
          as a percent of all possible exons.',p,br=TRUE)
  hwrite('',p,br=TRUE)
  hwrite(ExonsCoveredTable,p,style='text-align:right')
  hwrite('',p,br=TRUE)

  hwrite('Genomic features detected',p,heading=2)
  hwrite('Here we tabulate for each category the number of genomic features we detect with uniquely mapped reads. We calculate the percent, for each category, of features detected out of all possible features (ie. for genes, percent = (number of genes detected by reads)/(total number of genes annotated for the genome)).',p,br=TRUE)
  hwrite('',p,br=TRUE)
  hwrite(GenomicFeaturesDetectedTable,p,style='text-align:right')
  hwrite('',p,br=TRUE)
  hwriteImage(paste('images','GenomicFeaturesDetected.png',sep="/"),p,br=TRUE,link=paste("images","GenomicFeaturesDetected.pdf",sep="/"))
  hwrite('',p,br=TRUE)
  
  
  hwrite('Distribution of reads by genes (exonic)',p,heading=2)
  hwrite('For various coverages (measured by the number of reads mapping uniquely to the exons of a gene) we calculate the proportion of genes that obtain AT LEAST said coverage.',p,br=TRUE)
  hwrite('',p,br=TRUE)
  hwriteImage(paste('images','DistrOfReadsToGenesExonic.png',sep="/"),p,br=TRUE,link=paste("images","DistrOfReadsToGenesExonic.pdf",sep="/"))
  hwrite('',p,br=TRUE)
  hwrite('Here we restrict the x-axis to 2000 reads mapping to the exons of a gene.',p,br=TRUE)
  hwrite('',p,br=TRUE)
  hwriteImage(paste('images','DistrOfReadsToGenesExonic2000.png',sep="/"),p,br=TRUE,link=paste("images","DistrOfReadsToGenesExonic2000.pdf",sep="/"))
  hwrite('',p,br=TRUE)
  hwrite('Here we restrict the x-axis to 200 reads mapping to the exons of a gene.',p,br=TRUE)
  hwrite('',p,br=TRUE)
  hwriteImage(paste('images','DistrOfReadsToGenesExonic200.png',sep="/"),p,br=TRUE,link=paste("images","DistrOfReadsToGenesExonic200.pdf",sep="/"))
  
  hwrite('Distribution of RPKMs for genes (exonic)',p,heading=2)
  hwrite('For various RPKMs (the number of Reads mapping uniquely to the exons of a gene Per Kilobase
          of the mappable length of the gene (exonic) per Million mapped reads) we
          calculate the proportion of genes that obtain AT LEAST said RPKM.',p,br=TRUE)
  hwrite('',p,br=TRUE)
  hwriteImage(paste('images','DistrOfRPKMsToGenesExonic.png',sep="/"),p,br=TRUE,link=paste("images","DistrOfRPKMsToGenesExonic.pdf",sep="/"))
  hwrite('',p,br=TRUE)
  hwrite('Here we restrict the x-axis to RPKMs of 200.',p,br=TRUE)
  hwrite('',p,br=TRUE)
  hwriteImage(paste('images','DistrOfRPKMsToGenesExonic200.png',sep="/"),p,br=TRUE,link=paste("images","DistrOfRPKMsToGenesExonic200.pdf",sep="/"))
  hwrite('',p,br=TRUE)
  hwrite('Here we restrict the x-axis to RPKMs of 50.',p,br=TRUE)
  hwrite('',p,br=TRUE)
  hwriteImage(paste('images','DistrOfRPKMsToGenesExonic50.png',sep="/"),p,br=TRUE,link=paste("images","DistrOfRPKMsToGenesExonic50.pdf",sep="/"))
  
  closePage(p, splash=FALSE)
}

sanityCheck <- function(summary.alignTable, readFiltering, paired_ends){

  gsnap_processed <- sum(summary.alignTable[names(summary.alignTable) != "analyzed"])
  
  sanity_check <- c("All reads passed into gsnap are present in gsnap output" = "FAIL",
                    "There are more unique than multiple mappers " =            "FAIL")
  if(gsnap_processed == readFiltering[["Post-preprocess"]]){
    sanity_check[1] <- "PASS"
  }
  if(readFiltering[["Uniquely"]] > readFiltering[["Multiple"]]) {
    sanity_check[2] <- "PASS"
  }
  sanity_check <- data.frame(Status=sanity_check)
  return(sanity_check)
}

##' Make relative bar plots
##' 
##' @param data vector of raw, absolute counts
##' @param total number to normalize by, can be vector of same length as data
##' @param labels x-axes labels, category labels for data
##' @param title Title of the plot 
##' @param filename plots will be saved under [filename].png and [filename].pdf 
##' @param ylab label of y axis
##' @param cex.names scaling param of lables, passed to plot
##' @param ymax extent of y-axis
##' @return Nothing, creates two files instead
##' @export
##' @keywords internal
relativeBarPlot <- function(data, total, labels, title, filename, ylab="Percent", cex.names=0.9, ymax=100) {
  rel_data <- (data/total)*100
  .plot_rbp <-function() {
    rbp <- barplot(rel_data, names.arg=labels,
                   col="deepskyblue4", border=FALSE,
                   cex.names=cex.names,
                   ylim=c(0, ymax+5),
                   ylab=ylab,
                   main=title)
    ## print absolute number on top of bar
    text(rbp, rel_data+3, formatC(data, format="d", big.mark=','), col="black")
    
    percents <- c(sprintf("%1.0f", rel_data))
    percents <- sapply(percents, function(x) { paste(x, "%", sep="") })
    ## print relative numbers in top part of bar 
    text(rbp, rel_data-3,
         labels=percents,col="white")
    dev.off()
  }

  ## Print once in pdf and png
  if(is.finite(ymax)) {
    pdf(width=15*1.2, height=8, file=paste(filename, "pdf", sep="."))
    .plot_rbp()
    
    png(filename=paste(filename, "png", sep="."),width = 900*1.2, height = 480, units = "px")
    .plot_rbp()
  }
}

##' Make continuous plots of distribution function
##' 
##' @param df distribution function, given as absolute count and percent
##' @param filename plots will be saved under [filename].png and [filename].pdf 
##' @param xlab label of x axis
##' @param ylab label of y axis
##' @return Nothing, creates two files instead
##' @author Jens Reeder
##' @export
##' @keywords internal
plotDF <- function(df, ylab, xlab, filename) {

  df$percent <- (df$num/max(df$num))*100
  
  .plot <- function(){
    plot(df$count, df$percent,
         col="deepskyblue4",
         ylab=ylab, xlab=xlab,
         pch=16, cex=0.8,
         ylim=c(0,max(df$percent)),
         )
    
    lines(df$count, df$percent, col="deepskyblue4", type="l")
    dev.off()
  }
  
  pdf(width=15*1.2, height=8, file=paste(filename, "pdf", sep="."))
  .plot()
  
  png(width=900*1.2, height=480, filename=paste(filename, "png", sep="."))
  .plot()
}

how_many <- function(counts, decreasing=TRUE) {
  if(length(counts) == 0)  {
    df <- data.frame(num = 1, count = 1)
  }  else  {
    df <- data.frame(num=seq(1:length(counts)),
                     count=sort(counts,decreasing=decreasing))
    df <- df[length(counts):1,]
    df <- df[!duplicated(df$count),]
  }
  df
}

