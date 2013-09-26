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

  zipped.vcf <- bgzip(vcf.file)
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
