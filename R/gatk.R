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
##' @param maxheap Maximal heap space allocated for java,
##'               GATK recommends 4G heap for most of its apps
##' @return 0 for success, stops otherwise
##' @export
##' @author Jens Reeder
gatk <- function(gatk.jar.path=getOption("gatk.path"), method, args, maxheap="4g"){
  
  args <- cbind(paste("-XX:+UseSerialGC", ## keeps java from taking over a multicore server using lots of GC threads
                      paste0("-Xmx", maxheap),
                      "-jar", gatk.jar.path,
                      "-T", method),
                args)
  loginfo(paste("gatk.R/gatk: Calling gatk: java", paste(args, collapse=' ')))
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
                "--genotype_likelihoods_model", 'BOTH', ## choice of SNP, INDEL or BOTH
                gatk.params)

  gatk(gatk.jar.path=jar.path,
       method="UnifiedGenotyper", args=args)
  if(!file.exists(vcf.file)) stop("callVariantsGATK failed to create vcf file.")

  filtered.vcf.file <- filterGATKVars(vcf.file)
  zipped.vcf <- .zipUp(filtered.vcf.file)
  indexTabix(zipped.vcf, format="vcf")

  ## clean up unzipped vcf and index
  unlink(vcf.file)
  unlink(paste(vcf.file,'idx', sep='.'))

  return(zipped.vcf)
}

.zipUp <- function(filename){  
  if( grepl('gz$', filename)){
    zipped.vcf <- filename
  } else{
    zipped.vcf <- bgzip(filename, dest=paste0(filename, ".gz"))
  }
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

  if(!getConfig.logical("gatk.filter_repeats")){
    return(vcf.file)
  }
  
  genome     <- getConfig("alignReads.genome")
  rmask.file <- getConfig("analyzeVariants.rep_mask")

  loginfo("gatk.R::filterGATKVars: filtering variants for low complexity regions and repeats")  

  mask <- rtracklayer::import(rmask.file, asRangedData=FALSE)
  ## Variants are always annotated on +, so we make sure the repeats are also on +
  strand(mask) <- '+'

  variants <- VariantAnnotation::readVcf(vcf.file, genome)
  variants <- excludeVariantsByRegions(variants, mask)
 
  filename <- writeVcf(variants, vcf.file, index=FALSE)
  ## returns filename in the form xxx.vcf
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
realignIndels <- function(){
  loginfo("gatk.R/realignIndels: starting to realign indels")

  save_dir <- getConfig("save_dir")
  analyzed.bam <- getBams(save_dir)$"analyzed"

  realigned.bam <- safeExecute({
    realignIndelsGATK(analyzed.bam)
  }, memtracer=getConfig.logical("debug.tracemem"))
  
  file.rename(realigned.bam, analyzed.bam)
  file.rename(paste0(realigned.bam, ".bai"),
              paste0(analyzed.bam, ".bai"))

  loginfo("gatk.R/realignIndels: done")
}

##' Realign indels via GATK
##'
##' Realing indels using the GATK tools RealignerTargetCreator
##' and IndelRealigner.
##' Requires a GATK compatible genome with a name matching
##' the alignment genome to be installed in 'path.gatk_genome'
##'
##' Since GATKs IndelRealigner is not parallelized, we run it in
##' parallel per chromosome.
##' 
##' @param bam.file Path to bam.file
##' @return Path to realigned bam file
##' @author Jens Reeder
realignIndelsGATK <- function(bam.file){
  gatk.genomes <- getConfig("path.gatk_genomes")
  gatk.params  <- getConfig("gatk.params")
  jar.path     <- getConfig("path.gatk")
  genome       <- getConfig("alignReads.genome")
  num.cores    <- getConfig.integer("num_cores")  
  save.dir     <- getConfig('save_dir')
  genome       <- file.path(gatk.genomes, paste0(genome,".fa"))
  chunk.dir    <- getConfig('chunk_dir')
  
  if(any(grepl("-L", gatk.params))){
    ## user provided -L flag, so we can't parallelize by chr
    realigned.bam.files <- realignOne(bam.file, genome, gatk.params,
                                      jar.path, chunk.dir)
  } else {
    ## temporarily increase loglevel so gatk calls for each chr do not clutter stdout
    setLevel("WARN") 
    fun <- function(chr, ...) {
      realignOne(bam.file, genome, paste("-L", chr, gatk.params),
                 jar.path, chunk.dir)
    }
    listIterator.init( getGenomeSegments() )
    realigned.bam.files <- sclapply(listIterator.next, fun, max.parallel.jobs=num.cores)
    setLevel(getConfig("debug.level"))
  }
  
  realigned.bam.file <- file.path(save.dir, "bams",
                                  paste(getConfig('prepend_str'),
                                        "realigned", "bam", sep="."))

  bam <- catBams(realigned.bam.files, realigned.bam.file)
  return(bam)
}
  
realignOne <- function(bam.file, genome, gatk.params, jar.path, chunk.dir){

  tmp.dir <- createTmpDir(prefix="realigner", dir=chunk.dir)
  interval.file <- file.path(tmp.dir, "realigner.intervals")
  realigned.bam.file <- file.path(tmp.dir, "realigned.bam")
  
  args <- paste("-R", genome,
                "-I", bam.file,
                "-o", interval.file,
                gatk.params)
  
  gatk(gatk.jar.path=jar.path, method="RealignerTargetCreator", args=args, maxheap="2g")
  if(!file.exists(interval.file)) stop("realignIndelsGATK failed to create interval file.")
  
  args <- paste("-R", genome,
                "-I", bam.file,
                "-o", realigned.bam.file,
                "--targetIntervals", interval.file,
                gatk.params)

  gatk(gatk.jar.path=jar.path, method="IndelRealigner", args=args, maxheap="2g")
  if(!file.exists(realigned.bam.file)) stop("realignIndelsGATK failed to create bam file.")
  file.remove(interval.file)

  return(realigned.bam.file)
}

getGenomeSegments <- function() {
  return (seqnames( seqinfo( getGmapGenome())))
}

