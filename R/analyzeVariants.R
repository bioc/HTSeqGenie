##' @importFrom VariantTools VariantTallyParam callVariants variantGR2Vcf tallyVariants qaVariants
##' @importFrom VariantAnnotation writeVcf
NULL

##' Calculate and process Variants
##'
##' @title Calculate and process Variants
##' @return Nothing
##' @author Jens Reeder
##' @export
analyzeVariants <- function(){
  safeExecute({
    save.dir <- getConfig('save_dir')
    bam.file <- file.path(save.dir, 'bams',
                          paste(getConfig('prepend_str'),
                                'analyzed.bam', sep="."))
    reports.dir <- file.path(save.dir, "reports")
    
    variants <- wrap.callVariants(bam.file)
    if (length(variants$filtered.variants) > 0) {
      ## all these only make sense if we have actual variants
      writeVCF(variants$filtered.variants)
      processVariants(variants, reports.dir)
    } else {
      logwarn("No variants found. Check your inputs and log messages for errors")
    }
  }, memtracer=getConfig.logical("debug.tracemem"))
}

##' Process variant callings
##'
##' For a given eset of variants, calculates
##' variant depth, confusion matrix, transition ratio
##' and the mutation matrix
##' @param variants GrangesList of variants
##' @param reports.dir Directory to store the resulting images and html
##' @return Name of generated html file
##' @export
processVariants <- function(variants, reports.dir){

  filtered <- variants$filtered.variants

  depth.plot <- plotVariantDepth(filtered)
  transition.ratio <- getTransitionRatio(x = values(filtered)$ref,
                                         y = values(filtered)$alt)
  mutation.matrix <- plotMutationMatrix(filtered)

  .writeVariantSummary(variants, transition.ratio)
  html.file <- .writeVariantReport(reports.dir, depth.plot, mutation.matrix,
                                   transition.ratio)

  return(html.file)
}

##' Call Variants in the pipeline framework
##' 
##' A wrapper around VariantTools callVariant
##' framework.
##' 
##' @title Variant calling
##' @param bam.file Aligned reads as bam file
##' @return Variants as granges
##' @author Jens Reeder
##' @export
wrap.callVariants <- function(bam.file) {  
  loginfo("analyzeVariants.R/wrap.callVariants: Calling variants...")

  prepend_str <- getConfig('prepend_str')
  genome      <- getConfig("alignReads.genome")
  genome.dir  <- getConfig("path.gsnap_genomes")
  num.cores   <- getConfig.integer("num_cores")
  bqual       <- getConfig.numeric("analyzeVariants.bqual")
  save.dir    <- getConfig('save_dir')

  max.len <- NA
  z <- try( summary <- parseSummaries(save.dir, "summary_preprocess"), silent=TRUE)
  if(class(z)=="try-error"){
    ## We could implement a fallback to read from bam file here, but for now just stop
    stop("Can't read preprocess summary, so can't figure out read length.")
  }else{
    min.len <- summary[,'processed_min_read_length']
    max.len <- summary[,'processed_max_read_length']
    if(min.len != max.len) logwarn("Inconsistent read lengths can have a negative effect on variant calling as it uses the maximum read length for one of its filters")
  }
  
  tally.param <- VariantTallyParam(GmapGenome(genome = genome,
                                              directory = path.expand(genome.dir)),
                                   readlen = max.len,
                                   high_base_quality = bqual
                                   )
  ## we use the low level interface here, so we have access to the raw variants
  raw.variants <- tallyVariants(bam.file, tally.param, mc.cores=num.cores)
  qa.variants  <- qaVariants(raw.variants)
  variants     <- callVariants(qa.variants)

  loginfo("analyzeVariants.R/wrap.callVariants: Saving GRanges of raw and filtered variants...")
  ## TODO: sclapply this
  saveWithID(raw.variants,
             "raw_variants_granges",
             prepend_str,
             save_dir=file.path(getConfig('save_dir'), "results"),
             compress=FALSE)

  saveWithID(variants,
             "filtered_variants_granges",
             prepend_str,
             save_dir=file.path(getConfig('save_dir'), "results"),
             compress=FALSE)

  loginfo("analyzeVariants.R/wrap.callVariants: ...done")
  return(list(raw.variant = raw.variants,
              filtered.variants = variants))
}

## Generate a HTML report
.writeVariantReport <- function(reports.dir, depth.by.strand.plot, mutation.matrix,
                                transition.ratio){

  out.file <- file.path(reports.dir, "reportVariants.html")
  p=openPage(out.file)
  
  hwrite('Read Depth by Strand', p, heading=2)
  hwrite('Total read depth vs. fraction of reads on the positive strand', p)
  hwrite('', p, br=TRUE)
  hwriteImage(file.path("images", paste(basename(depth.by.strand.plot), "png", sep=".")),
              p, br=TRUE,
              link=file.path("images", paste(basename(depth.by.strand.plot), "pdf", sep=".")))
  hwrite('', p, br=TRUE)
  
  hwrite('Mutation matrix', p, heading=2)
  hwrite('Dsiplays the number of times a variant changed one nucleotide to another.', p, br=TRUE)
  hwrite('', p, br=TRUE)
  hwriteImage(file.path("images", "mutation_matrix.png"),
              p, br=TRUE,
              link=file.path("images", "mutation_matrix.pdf"))
  hwrite('', p, br=TRUE)
  hwrite(as(t(mutation.matrix), "matrix"), p)
  hwrite('', p, br=TRUE)

  hwrite('Nucleotide Transition vs Transversion Ratio', p, heading=2)
  hwrite('The ratio obtained by dividing the number of nucleotide transitions by the total number of variants.')
  hwrite('"Transitions" are defined as a pyramidine changing to another pyramidine (A <-> G) or a purine to another purine (C <-> T).', p, br=TRUE)
  hwrite('', p, br=TRUE)
  hwrite(sprintf("%.4f", transition.ratio), p)
  hwrite('', p, br=TRUE)
  
  closePage(p, splash=FALSE)
  return(out.file)
}

.writeVariantSummary <- function(variants, transition.ratio){

  raw      <- variants$raw.variants
  filtered <- variants$filtered.variants

  summary_variants <- list()
  summary_variants$'raw_variants_count' <-  length(raw)
  summary_variants$'filtered_variants_count' <-  length(filtered)
  summary_variants$'transitions_transversions_ratio' <- transition.ratio
  
  df <- data.frame(name=names(summary_variants), value=unlist(summary_variants)) 
  saveWithID(df, "summary_variants", id=getConfig('prepend_str'),
             save_dir=file.path(getConfig('save_dir'), "results"), format="tab")

  return(summary_variants)
}  

##' Write variants to VCF file
##' 
##' @title writeVCF
##' @param variants.granges Genomic Variants as Granges object
##' @return VCF file name
##' @author Jens Reeder
##' @export
writeVCF <- function(variants.granges){
  loginfo("analyzeVariants.R/writeVCF: writing vcf file...")
  
  vcf.filename <- file.path(getConfig('save_dir'), "results",
                            paste(getConfig('prepend_str'),
                                  "_variants.vcf", sep=""))

  vcf <- variantGR2Vcf(variants.granges,                               
                       sample.id=getConfig('alignReads.sam_id'),
                       project="") 
  ## this is VariantAnnotation::writeVcf
  writeVcf(vcf, vcf.filename, index=TRUE)

  loginfo("analyzeVariants.R/writeVCF: ...done")
  return(vcf.filename)
}

##' Plot the variant depth
##'
##' Plot the variant depth
##' @title Plot the variant depth
##' @param variants Variants as computed by SVNsOmuC
##' @return Name of plot file
##' @author Cory Barr, Jens Reeder
##' @export
plotVariantDepth <- function(variants){
  loginfo("analyzeVariants.R/plotVariantDepth: Plotting depth by strand...")
  depth_by_strand_plot <- file.path(getConfig('save_dir'),
                                    "reports", "images", "depth_by_strand")
  if(!file.exists(dirname(depth_by_strand_plot)))
    dir.create(dirname(depth_by_strand_plot), recursive=TRUE)

  plotDepthByStrand(variants, file=depth_by_strand_plot)
  loginfo("analyzeVariants.R/plotVariantDepth: ...done")
  return(depth_by_strand_plot)
}

##' Calculate Confusion Matrix
##'
##' Calculate Confusion Matrix from Granges object containing variants
##' @title Calculate Confusion Matrix
##' @param variants Variants as GRanges
##' @return confusion matrix as table
##' @author Cory Barr, Jens Reeder
##' @export
calcConfusionMatrix <-function(variants){
  loginfo("analyzeVariants.R/calcConfusionMatrix: Calculating nucleotide confusion matrix for variants...")
  confusion_matrix <- table(factor(as.character(values(variants)$ref),
                                   levels=c('A','C','G','T')),
                            factor(as.character(values(variants)$alt),
                                   levels=c('A','C','G','T')))
  loginfo("analyzeVariants.R/calcConfusionMatrix: ...done")
  return(confusion_matrix)
}

##' Plot Mutation Matrix
##'
##' Plot Mutation Matrix
##' @title Plot Mutation Matrix
##' @param variants Variants as Granges computed by SNVsOmuC
##' @return Nothing, stores plots as files
##' @author Jens Reeder, Jeremiah Degenhardt
##' @export
plotMutationMatrix <- function(variants){

  filename <- file.path(getConfig('save_dir'), "reports", "images",
                        "mutation_matrix")
  z <- calcConfusionMatrix(variants)

  .plot.mosaic <- function(){
    mosaicplot(z, col = c("firebrick3", "dodgerblue3", "forestgreen", "goldenrod3"),
               cex.axis = 2,
               main = "Mutation Matrix",
               xlab= "Reference allele",
               ylab= "Variant allele")
    dev.off()
  }
  pdf(height=8,width=15,file=paste(filename, "pdf", sep="."))
  .plot.mosaic()

  png.file <- paste(filename, "png", sep=".")
  png(filename= png.file,width = 900, height = 480, units = "px")
  .plot.mosaic()
  return(z)
}

##' Get Nucleotide Transition Ratio
##'
##' Get Nucleotide Transition Ratio
##' @title Get Nucleotide Transition Ratio
##' @param x DNAStringSet
##' @param y DNAStringSet
##' @return numeric vector of length 1
##' @author Cory Barr
##' @export
getTransitionRatio <- function(x, y) {
  loginfo("analyzeVariants.R/calcTransitionRatio: Calculating nucleotide transition ratio...")
  pasted <- paste(as.character(x), as.character(y), sep='')
  ratio <- mean(pasted == 'AG' | pasted == 'GA' | pasted == 'CT' | pasted == 'TC')
  loginfo("analyzeVariants.R/calcTransitionRatio: ...done")
  return(ratio)
} 
