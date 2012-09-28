##' Parse a progress.log file
##'
##' @title Parse a progress.log file
##' @param filename A character string, pointing to a "progres.log" file
##' @param plot A logical, indicating whether to draw a barplot. Default is FALSE.
##' @param main A character string, indicating the main label of the barplot.
##' @param maxt A numeric, indicating the maximal y-range of the barplot. Default is 30.
##' @return A vector containing processing times (in h) and peakmem (in GBi)
##' @author Gregoire Pau
##' @export
##' @keywords internal
parseProgressLog <- function(filename, plot=FALSE, main="Processing times", maxt=30) {
  ## read progress.log file
  prog <- readLines(filename)
  prog <- strsplit(prog, " ")
  prog <- lapply(prog, function(z) c(z[1:2], paste(z[-(1:2)], collapse=" ")))               
  prog <- do.call(rbind, prog)
  
  ## compute durations
  df <- data.frame(date=as.POSIXct(paste(prog[,1], prog[,2])), msg=prog[,3], stringsAsFactors=FALSE)
  duration <- function(txt1, txt2) difftime(df$date[grep(txt1, df$msg)], df$date[grep(txt2, df$msg)], units="hours")
  total <- duration("Pipeline run successful", "Hi!")
  preprocess <- duration("preprocessReads.R/preprocessReads: done", "preprocessReads.R/preprocessReads: starting...")
  align <- duration("alignReads.R/alignReads: done", "alignReads.R/alignReads: starting")
  variants <- duration("analyzeVariants.R/calcConfusionMatrix: ...done", "analyzeVariants.R/wrap.callVariants: Calling variants...")
  
  ## build dat
  duration <- c(total=total, preprocess=preprocess, align=align, variants=variants)
  duration <- c(duration, other=total-sum(duration[-match("total", names(duration))]))

  ## plot
  if (plot) {
    legend <- paste(round(duration, 2), "h")
    col <- c("#ffffaa", "#ffaaaa", "#aaffff", "#aaffaa", "#aaaaaa")
    mp <- barplot(duration, col=col, ylab="Duration (h)", ylim=c(0, maxt), main=main)
    text(mp, duration + 1, legend, xpd = TRUE)
  }

  ## compute peakmem
  peakmem <- grep("wired.mem=", df$msg, value=TRUE)
  peakmem <- as.numeric(substr(peakmem, regexpr("wired.mem=", peakmem)+10, regexpr("\ GBi", peakmem)-1))
  peakmem <- max(peakmem)
  
  c(duration=duration, peakmem=peakmem)
}
