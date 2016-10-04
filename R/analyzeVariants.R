##' Calculate and process Variants
##'
##' @title Calculate and process Variants
##' @return Nothing
##' @author Jens Reeder
##' @export
analyzeVariants <- function(){
  loginfo("analyzeVariants/analyzeVariants: starting ...")
  
  ## get config parameters
  save_dir <- getConfig('save_dir')
  analyzeVariants.method <- getConfig("analyzeVariants.method", stop.ifempty=TRUE)
  prepend_str <- getConfig("prepend_str")
  bam.file <- getAnalyzedBamFile()
  
  safeExecute({
    if(analyzeVariants.method=='GATK'){
      vcf.filename <- callVariantsGATK(bam.file)
    }
    if(analyzeVariants.method=='VariantTools'){ 
      vcf.filename <- callVariantsVariantTools(bam.file)
    }
    ## if we don't have any variants, we don't have a vcf file
    if (!is.null(vcf.filename)) {
      vcfstat <- vcfStat(vcf.filename)
      df <- data.frame(name=names(vcfstat), value=vcfstat, stringsAsFactors=FALSE)
      saveWithID(df, "summary_variants", id=prepend_str, save_dir=file.path(save_dir, "results"), format="tab")

      if(!is.null(getConfig("analyzeVariants.vep_options"))){
        annotateVariants(vcf.filename)
      }
    }
  }, memtracer=getConfig.logical("debug.tracemem"))
  loginfo("analyzeVariants/analyzeVariants: done")
}

callVariantsVariantTools <- function(bam.file) {
  variants <- wrap.callVariants(bam.file)

  filename <- file.path(getConfig('save_dir'), "results",
                        paste(getConfig('prepend_str'),
                              ".variants.vcf", sep=""))
  vcf.filename <- writeVCF(variants$filtered.variants, filename=filename)

  if(getConfig.logical("analyzeVariants.genotype")){
    .callGenotypes(variants$raw.variants)
  } 
  return(vcf.filename)
}

##' Call Variants in the pipeline framework
##' 
##' A wrapper around VariantTools callVariant
##' framework.
##' 
##' @title Variant calling
##' @param bam.file Aligned reads as bam file
##' @return Variants as Vranges
##' @author Jens Reeder
##' @export
##' @importFrom VariantTools TallyVariantsParam VariantPostFilters postFilterVariants VariantCallingFilters
##' @importMethodsFrom VariantTools callVariants tallyVariants
##' @importFrom BiocParallel MulticoreParam
wrap.callVariants <- function(bam.file) {
  ## get config parameters
  num.cores   <- getConfig.integer("num_cores")

  ## tally variants
  loginfo("analyzeVariants.R/wrap.callVariants: Tallying variants...")  
  tally.param  <- buildTallyParam()
  tally.variants <- tallyVariants(bam.file, tally.param,
                                  BPPARAM = MulticoreParam(workers=num.cores))

  loginfo("analyzeVariants.R/wrap.callVariants: calling variants...")
  variants <- callVariants(tally.variants, calling.filters=buildCallingFilters())
  variants <- .postFilterVariants(variants)
  
  loginfo("analyzeVariants.R/wrap.callVariants: Saving GRanges of raw and filtered variants...")

  .saveVTvars(tally.variants, variants)
  
  loginfo("analyzeVariants.R/wrap.callVariants: ...done")
  return(list(raw.variants = tally.variants,
              filtered.variants = variants))
}

##' @importFrom gmapR GmapGenome
getGmapGenome <- function(){
  genome      <- getConfig("alignReads.genome")
  genome.dir  <- getConfig("path.gsnap_genomes")
  
  gmap.genome <- GmapGenome(genome = genome,
                            directory = path.expand(genome.dir))
  return(gmap.genome)
}

##' Build tally parameters
##'
##' @title Build tally parameters
##' @return a \code{VariantTallyParam} object
##' @importFrom GenomicRanges tileGenome
##' @keywords internal
##' @author Gregoire Pau
buildTallyParam <- function(){
  indels      <- getConfig.logical("analyzeVariants.indels")
  bqual       <- getConfig.numeric("analyzeVariants.bqual")
  rmask.file  <- getConfig("analyzeVariants.rep_mask")
  positions.file <- getConfig("analyzeVariants.positions")
  gmap.genome <- getGmapGenome()
  
  ## we use the low level interface here, so we have access to the raw variants
  args <- list(gmap.genome,
               high_base_quality=as.integer(bqual),
               indels=indels)

  if(!is.null(rmask.file)){
    loginfo("analyzeVariants.R/buildTallyParam: Using repeat masker track for tally filtering")
    mask <- rtracklayer::import(rmask.file, asRangedData=FALSE)
    args <- c(args, mask=mask)
  }

  if(!is.null(positions.file)){
    loginfo("analyzeVariants.R/buildTallyParam: restricting variant calls using 'which'")
    regions <- readRDS(positions.file)
    args <- c(args, which=regions)
  } else{
    ## tile width of 1e7 turns out to be a good compromise between speed and memory.
    ## Large enough to not overflow even for large WGS and fast enough for samples with
    ## just a few reads
    nbreads <- tryCatch(getNumericVectorDataFromFile(getConfig("save_dir"), "summary_alignment"),
                        error=function(...) NA)
    nbreads <- nbreads["analyzed"]
    if (!is.na(nbreads) & nbreads>100e6) {
      loginfo("analyzeVariants.R/buildTallyParam: using genome tiling for memory saving")
      tiles <- unlist(tileGenome(seqinfo(gmap.genome), tilewidth=1e7))
      args <- c(args, which=tiles)
    }
  }
  
  tally.param <- do.call(TallyVariantsParam, args)

  return(tally.param)
}

##Given a list of calling filter names, construct actual Filter objects
buildCallingFilters <- function(){

  calling.filter.names <- getConfig.vector("analyzeVariants.callingFilters")
  calling.filters <- VariantCallingFilters()
  
  return( calling.filters[ calling.filter.names ] )
}

.postFilterVariants <- function(variants){
  loginfo("analyzeVariants.R/wrap.postFilerVariants: post filtering variants...")
  dbsnp.file  <- getConfig("analyzeVariants.dbsnp")
  filter.names <- getConfig.vector("analyzeVariants.postFilters")

  if (!is.null(dbsnp.file) & "avgNborCount"%in%filter.names) {
    loginfo("analyzeVariants.R/wrap.postFilterVariants: Using dbsnp whitelist for variant filtering.")
    dbsnp <- get(load(dbsnp.file))
    all.filters <- VariantPostFilters(whitelist=dbsnp)
  } else {
    all.filters <- VariantPostFilters()
  }
  filters = all.filters[ filter.names ]
  
  vars <- postFilterVariants(variants, post.filters=filters)
  
  return(vars)
}

.saveVTvars <- function(raw.vars, filtered.vars) {
  prepend_str <- getConfig('prepend_str')
  save_dir    <- getConfig('save_dir')

  saveWithID(raw.vars,
             "raw_variants",
             prepend_str,
             save_dir=file.path(getConfig('save_dir'), "results"),
             compress=FALSE)

  saveWithID(filtered.vars,
             "filtered_variants",
             prepend_str,
             save_dir=file.path(getConfig('save_dir'), "results"),
             compress=FALSE)
}

##' @importFrom gmapR GmapGenome
##' @importFrom rtracklayer BigWigFile
##' @importMethodsFrom VariantTools callGenotypes
##' @importFrom VariantTools CallGenotypesParam
.callGenotypes <- function(variants){
  loginfo("analyzeVariants.R/.callGenotypes: calling genotypes...")  

  save.dir    <- getConfig("save_dir")
  prepend.str <- getConfig("prepend_str")
  ## Genotyping is heavy on memory, 14GB, per chunk
  ## Until we optimize this is VT, we simpley trade speed for less memory
  num.cores   <- max(1, getConfig.integer("num_cores") %/% 2)
  gmap.genome <- getGmapGenome()
  
  bwfile <- file.path(save.dir, "results", paste(prepend.str, "coverage", "bw", sep="."))
  
  cg.param <- CallGenotypesParam(genome=gmap.genome)

  genotypes <- callGenotypes(variants, BigWigFile(bwfile), param=cg.param,
                             BPPARAM = MulticoreParam(workers=num.cores))

  genotypes <- genotypes[, c("GT", "GQ", "PL", "MIN_DP")]
  geno.file <- file.path(save.dir, "results", paste0(prepend.str, ".genotypes.vcf"))  
  writeVCF(genotypes, filename=geno.file)
  
  loginfo("analyzeVariants.R/.callGenotypes: done")  
  return(geno.file)
}  

##' Write variants to VCF file
##' 
##' @title writeVCF
##' @param variants.vranges Genomic Variants as VRanges object
##' @param filename Name of vcf file to write
##' @return VCF file name
##' @author Jens Reeder
##' @export
##' @importMethodsFrom VariantAnnotation writeVcf sampleNames sampleNames<-
##' @importFrom Rsamtools bgzip indexTabix
writeVCF <- function(variants.vranges, filename){
  loginfo("analyzeVariants.R/writeVCF: writing vcf file...")

  vcf.gz <- NULL
  if (length(variants.vranges)>0) {
    ## fill in sample name if missing from vranges
    if(all(is.na(sampleNames(variants.vranges)))){
      sam_id <- getConfig("alignReads.sam_id")
      if (is.null(sam_id)) sam_id <- ""
      sampleNames(variants.vranges) <- sam_id
    }
    ## build vcf object
    ## make sure we sort, as otherwise writing the index might crash
    vcf <- as(sort(variants.vranges), "VCF")
    ## if we let writeVcf index, it makes us a .bgz file which breaks IGV
    ## thus we compress and index ourselves using the accepted .gz suffix
    vcf.filename <- writeVcf(vcf, filename, index=FALSE)
    vcf.gz <- bgzip(vcf.filename, dest=paste0(vcf.filename, ".gz"))
    indexTabix(vcf.gz, format = "vcf")
    unlink(vcf.filename)
  } 
  loginfo("analyzeVariants.R/writeVCF: ...done")
  return(vcf.gz)
}

##' Compute stats on a VCF file
##'
##' @title Compute stats on a VCF file
##' @param vcf.filename A character pointing to a VCF (or gzipped VCF) file
##' @return A numeric vector
##' @author Gregoire Pau
##' @export
##' @importFrom VariantAnnotation scanVcf
vcfStat <- function(vcf.filename) {
  ## read vcf filename
  vcf <- scanVcf(vcf.filename)[[1]]

  ## general
  nb.variants <- length(vcf$REF)
  zsnv <- width(vcf$REF)==1 & nchar(vcf$ALT)==1
  nb.snvs <- sum(zsnv)
  nb.indels <- nb.variants - nb.snvs
  
  ## snvs nuc
  nuc <- c("A", "C", "G", "T")
  ref <- factor(as.character(vcf$REF[zsnv]), nuc)
  alt <- factor(vcf$ALT[zsnv], nuc)
  tab <- as.numeric(table(ref, alt))
  ntab <- data.frame(ref=nuc, alt=rep(nuc, each=length(nuc)), stringsAsFactors=FALSE)
  names(tab) <- paste("nb.snvs.", ntab$ref, ">", ntab$alt, sep="")
  tab <- tab[ntab$ref!=ntab$alt]
  
  ## snvs freq
  df <- data.frame(nref=numeric(0), nalt=numeric(0))
  try({
    if (!is.null(vcf$GENO$AD)) {
      df <- as.data.frame(do.call(rbind, vcf$GENO$AD[zsnv]))
      colnames(df) <- c("nref", "nalt")
    }
  }, silent=TRUE)
  df$freqalt <- df$nalt/(df$nalt+df$nref)
  qq <- seq(0, 1, length=21)
  q.snvs.nalt <- quantile(df$nalt, qq, na.rm=TRUE)
  names(q.snvs.nalt) <- paste(sprintf("q%03d", round(qq*100)), ".snvs.nalt", sep="")
  q.snvs.freqalt <- quantile(df$freqalt, qq, na.rm=TRUE)
  names(q.snvs.freqalt) <- paste(sprintf("q%03d", round(qq*100)), ".snvs.freqalt", sep="")
  
  ## final
  c(nb.variants=nb.variants, nb.snvs=nb.snvs, nb.indels=nb.indels,
    tab,
    q.snvs.nalt, q.snvs.freqalt)
}

checkVEP <- function(){
  retval <- system2("variant_effect_predictor.pl", "--help",
                    stderr=FALSE, stdout=FALSE)
  if (retval==0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

##' @export
annotateVariants <- function(vcf.file){  
  loginfo("analyzeVariants/annotateVariants: starting...")
  num.cores <- getConfig.integer("num_cores")
  options   <- getConfig("analyzeVariants.vep_options")
  out.file  <- file.path(getConfig('save_dir'), "results",
                         paste(getConfig('prepend_str'),
                               ".variants_vep.vcf.bgz", sep=""))
  
  if(num.cores>1){
    options <- paste(options, "--fork", num.cores)
  }
  options <- paste(options, "--input_file", vcf.file, "--tabix", "--output_file", out.file)

  retval <- system2("variant_effect_predictor.pl", options, stderr=FALSE, stdout=FALSE)
  if (retval!=0){ 
    stop("variant_effect_predictor.pl failed." , options)
  }
  loginfo("analyzeVariants/annotateVariants: done")
  return(out.file)
}
