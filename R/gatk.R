##' Run a command from the GATK
##'
##' Execute the GATK jar file using
##' the method specified as arg.
##' Stops if the command executed fails.
##' 
##' @title gatk
##' @param gatk.jar.path Path to the gatk jar file
##' @param method Name of the gatk method, e.g. UnifiedGenotyper 
##' @param args additional args passed to gatk
##' @return 0 for success, stops otherwise
##' @export
##' @author Jens Reeder
gatk <- function(gatk.jar.path=getOption("gatk.path"), method, args){
  
  args <- cbind(paste("-jar", gatk.jar.path,
                      "-T", method),
                args)
  loginfo(paste("analyzeVariants/gatk: Calling gatk: java", paste(args, collapse=' ')))
  retcode <- system2('java', args, stdout=FALSE)
  if(retcode != 0){
    stop( paste("GATK command [ java", paste(args, collapse=" "), "] failed."))
  }
  retcode
} 

##' Variant calling via GATK
##'
##' Call variants via GATK using the pipeline framework.
##' Requires a GATK compatible genome with a name matching
##' the alignment genome to be installed in 'path.gatk_genome'
##' 
##' @param bam.file Path to bam.file
##' @return Path to variant file
##' @author Jens Reeder
callVariantsGATK <- function(bam.file){
  gatk.genomes <- getConfig("path.gatk_genomes")
  gatk.params  <- getConfig("gatk.params")
  jar.path     <- getConfig("path.gatk")
  genome       <- getConfig("alignReads.genome")
  num.cores    <- getConfig.integer("num_cores")  
  save.dir     <- getConfig('save_dir')
              
  vcf.file <- file.path(save.dir, "results",
                        paste(getConfig('prepend_str'),
                              "variants","vcf", sep="."))

  genome <- file.path(gatk.genomes, paste0(genome,".fa"))
  args <- paste("--num_threads", min(4,num.cores),
                ## enabling more threads does not seem to have a positive effect on runtime
                "-R", genome,
                "-I", bam.file,
                "-o", vcf.file,
                gatk.params)

  gatk(gatk.jar.path=jar.path,
       method="UnifiedGenotyper", args=args)
  if(!file.exists(vcf.file)) stop("callVariantsGATK failed to create vcf file.")

  if(getConfig.logical("gatk.filter_repeats")){
    filtered.vcf.file <- filterGATKVars(vcf.file)
  } else {
    filtered.vcf.file <- vcf.file
  }
  
  if( grepl('bgz$', filtered.vcf.file)){
    zipped.vcf <- filtered.vcf.file
  } else{
    zipped.vcf <- bgzip(filtered.vcf.file)
  }
  indexTabix(zipped.vcf, format="vcf")

  ## clean up unzipped vcf and index
  unlink(vcf.file)
  unlink(paste(vcf.file,'idx', sep='.'))

  return(zipped.vcf)
}

##' Check for the GATK jar file
##
##' @param path Path to the GATK jar file
##' @return TRUE if tool can be called, FALSE otherwise
##' @export
checkGATKJar <- function(path=getOption("gatk.path")){

  if(is.null(path)) return(FALSE)
  if (!file.exists(path))
    return(FALSE)

  retval <- system2("java", paste("-jar", path, "-h"), stderr=FALSE, stdout=FALSE)
  if (retval==0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

filterGATKVars <- function(vcf.file){  
  genome     <- getConfig("alignReads.genome")
  rmask.file <- getConfig("analyzeVariants.rep_mask")

  loginfo("gatk.R::filterGATKVars: filtering variants for low complexity regions and repeats")  

  mask <- rtracklayer::import(rmask.file, asRangedData=FALSE)
  ## Variants are always annotated on +, so we make sure the repeats are also on +
  strand(mask) <- '+'

  variants <- VariantAnnotation::readVcf(vcf.file, genome)
  variants <- excludeVariantsByRegions(variants, mask)
 
  filename <- writeVcf(variants, vcf.file, index=TRUE)
  ## returns filename in the form vcf.file.bgz
  return(filename)
}

##' Filter variants by regions
##'
##' This function can be used to filter variants in a given region,
##' e.g. low complexity and repeat regions
##' @title Filter variants by regions
##' @param variants Variants as Vranges, GRanges or VCF object
##' @param mask region to mask, given as GRanges
##' @return The filtered variants 
##' @author Jens Reeder 
excludeVariantsByRegions <- function(variants, mask) {

  if(class(variants)=="CollapsedVCF"){
    len <- length(rowData(variants))
  }else{
    len <- length(variants)
  }
  keep = rep(TRUE,len)
  ov = findOverlaps(variants, mask)
  keep[ov@queryHits]=FALSE
  variants = variants[keep,]
 
  return(variants)
}

##' Realign indels in pipeline context
##'
##' High level function call to realign
##' indels in the analyzed.bam file using GATK
##' @title realignIndels
##' @return Nothing
##' @author Jens Reeder
##' @export
##' @importFrom tools file_path_sans_ext
realignIndels <- function(){
  loginfo("gatk.R/realignIndels: starting to realign indels")

  save_dir <- getConfig("save_dir")
  analyzed.bam <- getBams(save_dir)$"analyzed"

  safeExecute({
    realigned.bam <- realignIndelsGATK(analyzed.bam)
    file.rename(realigned.bam, analyzed.bam)
    # Realigner put a file with bai instead of bam right next to the bam file.
    # In our world we call these files .bam.bai files
    file.rename(paste0(file_path_sans_ext(realigned.bam), ".bai"),
                paste0(analyzed.bam, ".bai"))
  }, memtracer=getConfig.logical("debug.tracemem"))

  loginfo("gatk.R/realignIndels: finished realigning Indels")
}

##' realign indels via GATK
##'
##' Realing indels using the GATK tools RealignerTargetCreator
##' and IndelRealigner.
##' Requires a GATK compatible genome with a name matching
##' the alignment genome to be installed in 'path.gatk_genome'
##' 
##' @param bam.file Path to bam.file
##' @return Path to realigned bam file
##' @author Jens Reeder
realignIndelsGATK <- function(bam.file){
  gatk.genomes <- getConfig("path.gatk_genomes")
  jar.path     <- getConfig("path.gatk")
  genome       <- getConfig("alignReads.genome")
  num.cores    <- getConfig.integer("num_cores")  
  save.dir     <- getConfig('save_dir')
  genome       <- file.path(gatk.genomes, paste0(genome,".fa"))
  
  interval.file <- file.path(save.dir, "bams",
                             paste(getConfig('prepend_str'),
                                   "realigner.intervals", sep="."))

  args <- paste("--num_threads", min(4,num.cores),
                "-R", genome,
                "-I", bam.file,
                "-o", interval.file)
  
  gatk(gatk.jar.path=jar.path, method="RealignerTargetCreator", args=args)
  if(!file.exists(interval.file)) stop("realignIndelsGATK failed to create interval file.")
  
  realigned.bam.file <- file.path(save.dir, "bams",
                                  paste(getConfig('prepend_str'),
                                        "realigned", "bam", sep="."))

  args <- paste("-R", genome,
                "-I", bam.file,
                "-o", realigned.bam.file,
                "--targetIntervals", interval.file)

  gatk(gatk.jar.path=jar.path, method="IndelRealigner", args=args)
  if(!file.exists(realigned.bam.file)) stop("realignIndelsGATK failed to create bam file.")
  file.remove(interval.file)

  return(realigned.bam.file)
}
