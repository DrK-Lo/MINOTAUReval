#' Creates a manhattan plot
#' 
#' Creates a manhattan plot from PLINK assoc output (or any data frame with 
#' chromosome, position, and p-value).
#' 
#' @param x A data.frame with columns "BP," "CHR," "P," and optionally, "SNP."
#' @param chr A string denoting the column name for the chromosome. Defaults to 
#'   PLINK's "CHR." Said column must be numeric. If you have X, Y, or MT 
#'   chromosomes, be sure to renumber these 23, 24, 25, etc.
#' @param bp A string denoting the column name for the chromosomal position. 
#'   Defaults to PLINK's "BP." Said column must be numeric.
#' @param stat A string denoting the column name for the statistic being plotted. Defaults to 
#'   PLINK's "P." Said column must be numeric.
#' @param snp A string denoting the column name for the SNP name (rs number). 
#'   Defaults to PLINK's "SNP." Said column should be a character.
#' @param col A character vector indicating which colors to alternate.
#' @param alpha A number for data transparency between 0 and 1, where 1 is no transparency.
#' @param chrlabs A character vector equal to the number of chromosomes
#'   specifying the chromosome labels (e.g., \code{c(1:22, "X", "Y", "MT")}).
#' @param outThresh1 Alpha value to be used to draw an outlier threshold line.
#' @param outType1 How alpha is split across the tails of the distribution. "two-tail" 
#'   indicates it is split evenly, while "lower" and "upper" represent alpha contained
#'   at the bottom or top of the distribution.
#' @param outThresh2 A second alpha value to be used to draw another outlier threshold line.
#' @param outType2 How alpha is split across the tails of the distribution. "two-tail" 
#'   indicates it is split evenly, while "lower" and "upper" represent alpha contained
#'   at the bottom or top of the distribution.
#' @param outData Vector of data to use to infer the outler threshold in 
#' \code{outThresh1} and \code{outThresh2}. Defaults to the \code{stat} column
#'   of the input data.frame. User should use this to supply the full data distribution
#'   when plotting only a subset of the full dataset to visualize proper outlier
#'   thresholds.
#' @param high_reg A data.frame of inclusive chromosomal regions to highlight with 
#'   chromosome number in column 1, region start position (lower bound) in column 2, 
#'   and region end position (upper bound) in column 3.
#' @param high_reg_col Color to use when highlighting chromosomal regions.
#' @param high_snp_pos A data.frame of SNPs positions in your dataset to highlight,
#'   with chromosome number in column 1 and position in column 2. These SNPs 
#'   should all be in your dataset.
#' @param high_snp_pos_col Color to use when highlighting SNPs from positions.
#' @param high_snp_id A character vector of SNPs in your dataset to highlight. 
#'   These SNPs should all be in your dataset.
#' @param high_snp_id_col Color to use when highlighting SNPs from IDs.
#' @param logp If TRUE, the -log10 of the statistic is plotted. This is useful when 
#'   plotting p-values, as it isn't very useful to plot raw p-values. Plotting the 
#'   raw value could be useful for other genome-wide plots, for example, peak heights, 
#'   bayes factors, test statistics, other "scores," etc.
#' @param log2 If TRUE, the log2 of the statistic is plotted.
#' @param log10 If TRUE, the log10 of the statistic is plotted.
#' #' @param annotatePval If set, SNPs below this p-value will be annotated on the plot.
#' #' @param annotateTop If TRUE, only annotates the top hit on each chromosome that is below the annotatePval threshold. 
#' @param ... Arguments passed on to other plot/points functions
#'   
#' @return A manhattan plot.
#'   
#' @keywords visualization manhattan
#'   
#' @examples
#' manhattan(gwasResults)
#'   
#' @importFrom calibrate textxy  
#'   
#' @export


flex_manhattan <- function(x, chr="CHR", bp="BP", stat="P", snp="SNP", 
                      col=c("gray10", "gray60"), alpha=NULL, chrlabs=NULL,
                      outThresh1=NULL, outType1="two-tail", 
                      outThresh2=NULL, outType2="two-tail", outData=NULL,
                      high_reg=NULL, high_reg_col="darkgreen", 
                      high_snp_pos=NULL, high_snp_pos_col="green3", 
                      high_snp_id=NULL, high_snp_id_col="green3", 
                      logp=FALSE, log2=FALSE, log10=FALSE, 
#                      annotateOutliers=NULL, outlierType=NULL, 
#                      annotateTop=TRUE,
                      ...) {
 
  # Not sure why, but package check will warn without this.
  CHR=BP=STAT=index=NULL
  
  # Check for sensible dataset
  ## Make sure you have chr, bp and p columns.
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(stat %in% names(x))) stop(paste("Column", stat, "not found!"))
  ## warn if you don't have a snp column
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  ## make sure chr, bp, and p columns are numeric.
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[stat]])) stop(paste(stat, "column should be numeric."))
  
  # Create a new data.frame with columns called CHR, BP, and P.
  d=data.frame(CHR=x[[chr]], BP=x[[bp]], STAT=x[[stat]])
  
  # If the input data frame has a SNP column, add it to the new data frame you're creating.
  if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
  
  # If more than 10,000 datapoints, add some transparancy
  if (!is.null(alpha)) {
    trans=alpha*255  
  } else if (nrow(d) < 5000) {
    trans=1*255
  } else {
    trans=0.5*255
  }
  # Set positions, ticks, and labels for plotting
  ## Sort and keep only values where is numeric.
  #d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(STAT)))
  d <- d[order(d$CHR, d$BP), ]
  #d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
  if (logp) {
    d$logp <- -log10(d$STAT)
    ylabel=expression(-log[10](italic(stat)))
  } else if (log2) {
    d$logp <- log2(d$STAT)
    ylabel=expression(log[2](italic(stat)))
  } else if (log10) {
    d$logp <- log10(d$STAT)
    ylabel=expression(log[10](italic(stat)))
  } else {
    d$logp <- d$STAT
    ylabel=stat
  }
  d$pos=NA
  
  
  # Fixes the bug where one chromosome is missing by adding a sequential index column.
  d$index=NA
  ind = 0
  for (i in unique(d$CHR)){
    ind = ind + 1
    d[d$CHR==i,]$index = ind
  }
  
  # This section sets up positions and ticks. Ticks should be placed in the
  # middle of a chromosome. The a new pos column is added that keeps a running
  # sum of the positions of each successive chromsome. For example:
  # chr bp pos
  # 1   1  1
  # 1   2  2
  # 2   1  3
  # 2   2  4
  # 3   1  5
  nchr = length(unique(d$CHR))
  if (nchr==1) { ## For a single chromosome
    ## Uncomment the next two linex to plot single chr results in Mb
    #options(scipen=999)
    #d$pos=d$BP/1e6
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
    xlabel = paste('Chromosome',unique(d$CHR),'position')
    labs = ticks
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
        d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
      }
      # Old way: assumes SNPs evenly distributed
      # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
      # New way: doesn't make that assumption
      ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
    }
    xlabel = 'Chromosome'
    #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
    labs <- unique(d$CHR)
  }
  
  # Initialize plot
  xmax = ceiling(max(d$pos, na.rm=TRUE) * 1.03)
  xmin = floor(max(d$pos, na.rm=TRUE) * -0.03)
  
  if (max(d$logp, na.rm=TRUE) < 1 && min(d$logp, na.rm=TRUE) >= 0) {
    ymax = (ceiling(max(d$logp, na.rm=TRUE)) / (1/max(d$logp, na.rm=TRUE))) * 1.03
    ymin = (floor(min(d$logp, na.rm=TRUE)) / (1/min(d$logp, na.rm=TRUE))) * 1.03
  } else {
    ymax = ceiling(max(d$logp, na.rm=TRUE) * 1.03)
    ymin = floor(min(d$logp, na.rm=TRUE) * 1.03)
  }
  
  
  # The old way to initialize the plot
  # plot(NULL, xaxt='n', bty='n', xaxs='i', yaxs='i', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
  #      xlab=xlabel, ylab=expression(-log[10](italic(p))), las=1, pch=20, ...)
  
  
  # The new way to initialize the plot.
  ## See http://stackoverflow.com/q/23922130/654296
  ## First, define your default arguments
  def_args <- list(xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                   xlim=c(xmin,xmax), ylim=c(ymin,ymax),
                   xlab=xlabel, ylab=ylabel)
  ## Next, get a list of ... arguments
  #dotargs <- as.list(match.call())[-1L]
  dotargs <- list(...)
  ## And call the plot function passing NA, your ... arguments, and the default
  ## arguments that were not defined in the ... arguments.
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  
  # If manually specifying chromosome labels, ensure a character vector and number of labels matches number chrs.
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs)==length(labs)) {
        labs <- chrlabs
      } else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    } else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  
  # Add an axis. 
  if (nchr==1) { #If single chromosome, ticks and labels automatic.
    axis(1, ...)
  } else { # if multiple chrs, use the ticks and labels you created above.
    axis(1, at=ticks, labels=labs, ...)
  }
  
  # Create a vector of alternating colors
  col=rep(col, max(d$CHR))
  
  # Add points to the plot
  if (nchr==1) {
    with(d, points(pos, logp, pch=20, col=rgb(matrix(col2rgb(col[1]),ncol=3), alpha=trans, maxColorValue=255), ...))
  } else {
    # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
    icol=1
    for (i in unique(d$index)) {
      with(d[d$index==unique(d$index)[i], ], points(pos, logp, pch=20, col=rgb(matrix(col2rgb(col[icol]),ncol=3), alpha=trans, maxColorValue=255), ...))
      icol=icol+1
    }
  }
  
  # Gather correct data for inferring quantile outliers
  if (is.null(outData)) {
    threshData=d$STAT
  } else {
    threshData=outData
  }

  # Add outlier threshold lines
  if (!is.null(outThresh1)) {
    if (outType1 == "two-tail") {
        low_thresh = quantile(threshData, outThresh1/2, na.rm=TRUE)
        high_thresh = quantile(threshData, 1-(outThresh1/2), na.rm=TRUE)
        abline(h=low_thresh, col="orange", lwd=2)
        abline(h=high_thresh, col="orange", lwd=2)
      } else if (outType1 == "upper") {
        thresh = quantile(threshData, 1-outThresh1, na.rm=TRUE)
        abline(h=thresh, col="orange", lwd=2)
      } else if (outType1 == "lower") {
        thresh = quantile(threshData, outThresh1, na.rm=TRUE)
        abline(h=thresh, col="orange", lwd=2)
      } else {
        stop("You have inputted an incorrect outlier threshold type.")
      }
  }
  
  if (!is.null(outThresh2)) {
      if (outType2 == "two-tail") {
          low_thresh = quantile(threshData, outThresh2/2, na.rm=TRUE)
          high_thresh = quantile(threshData, 1-(outThresh2/2), na.rm=TRUE)
          abline(h=low_thresh, col="blue", lwd=2)
          abline(h=high_thresh, col="blue", lwd=2)
      } else if (outType2 == "upper") {
          thresh = quantile(threshData, 1-outThresh2, na.rm=TRUE)
          abline(h=thresh, col="blue", lwd=2)
      } else if (outType2 == "lower") {
          thresh = quantile(threshData, outThresh2, na.rm=TRUE)
          abline(h=thresh, col="blue", lwd=2)
      } else {
          stop("You have inputted an incorrect outlier threshold type.")
      }
  }

  # Highlight regions from a dataframe of positions (3 columns)
  if (!is.null(high_reg)) {
    if (any(!(high_reg[,1] %in% d$CHR))) warning("You're trying to highlight regions on chromsomes that don't exist in your results.")
    for (i in 1:nrow(high_reg)) {
      row = high_reg[i,]
      d.region=d[which(d$CHR %in% row[,1]), ]
      d.region=subset(d.region, BP >= row[,2] & BP <= row[,3])
      with(d.region, points(pos, logp, col=high_reg_col, pch=17, cex=1.5, ...))
    }
  }
  
  # Highlight snps from a dataframe of positions (2 columns)
  if (!is.null(high_snp_pos)) {
    if (any(!(high_snp_pos[,1] %in% d$CHR)) & any(!(high_snp_pos[,2] %in% d$POS))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.high_snp_pos=d[which(d$CHR %in% high_snp_pos[,1]), ]
    d.high_snp_pos=subset(d.high_snp_pos, BP %in% high_snp_pos[,2])
    with(d.high_snp_pos, points(pos, logp, col=high_snp_pos_col, pch=17, cex=1.5, ...))
  }
  
  # Highlight snps from a character vector
  if (!is.null(high_snp_id)) {
    if (any(!(high_snp_id %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.high_snp_id=d[which(d$SNP %in% high_snp_id), ]
    with(d.high_snp_id, points(pos, logp, col=high_snp_id_col, pch=17, cex=1.5, ...)) 
  }
  
  # Highlight top SNPs
#   if (!is.null(annotateOutliers)) {
#     if (outlierType == "two-tail") {
#         low_thresh = quantile(d$STAT, annotateOutliers/2, na.rm=TRUE)
#         high_thresh = quantile(d$STAT, 1-(annotateOutliers/2), na.rm=TRUE)
#         topHits = subset(d, STAT <= low_thresh & STAT >= high_thresh)
#     } else if (outlierType == "upper") {
#         thresh = quantile(d$STAT, 1-annotateOutliers, na.rm=TRUE)
#         topHits = subset(d, STAT >= thresh)
#     } else if (outlierType == "lower") {
#         thresh = quantile(d$STAT, annotateOutliers, na.rm=TRUE)
#         topHits = subset(d, STAT <= thresh)
#     }
#     # extract top SNPs at given p-val
# #    topHits = subset(d, STAT >= annotatePval)
#     par(xpd = TRUE)
#     # annotate these SNPs
#     if (annotateTop == FALSE) {
#       with(subset(d, STAT >= annotateOutliers, 
#            textxy(pos, STAT, offset = 0.625, labs = topHits$SNP, cex = 0.45), ...))
#     }
#     else {
#       # could try alternative, annotate top SNP of each sig chr
#       topHits <- topHits[order(topHits$STAT),]
#       topSNPs <- NULL
#       
#       for (i in unique(topHits$CHR)) {
#         
#         chrSNPs <- topHits[topHits$CHR == i,]
#         topSNPs <- rbind(topSNPs, chrSNPs[1,])
#         
#       }
#       textxy(topSNPs$pos, -log10(topSNPs$STAT), offset = 0.625, labs = topSNPs$SNP, cex = 0.5, ...)
#     }
#   }  
#   par(xpd = FALSE)
}