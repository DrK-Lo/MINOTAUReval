setwd("/Users/katie/Desktop/MINOTAUReval")
source("evalsims/ComparePlot2.R")
source("evalsims/distanceFunctionsOther.R")
source("evalsims/getEmpiricalP2.R")
source("evalsims/Getdf.R")


if(!("devtools" %in% installed.packages())){install.packages("devtools", dependencies=TRUE)}
if(!("rrcovNA" %in% installed.packages())){install.packages("rrcovNA")}
if(!("qvalue" %in% installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite("qvalue")
  }

library(qvalue)
library(rrcovNA)
library(devtools)
install_github("NESCent/MINOTAUR")#, build_vignettes=TRUE)
library(MINOTAUR)

headdir <- "~/Google Drive/MultiOutlierVisualization/practiceData/KatieSims/"
resultsdir <- "~/Google Drive/MultiOutlierVisualization/practiceData/KatieSimsResults/"
filelist <- list.files(headdir)
filelist <- filelist[-55]
length(filelist)
if (!length(filelist)==72){print("Errr: missing simulation files")}


summarydf <- NULL
for (i in 1:length(filelist)){
  print(c(i, filelist[i]))
  infile <- paste(headdir, filelist[i], sep="")
  outfile <- paste(resultsdir, filelist[i], sep="")
  dat <- read.table(infile, header=TRUE)
  cols <- which(names(dat) %in% c("rho", "log.bf", "xtx", "TW.Zscore"))
  dat2 <- dat[dat$UseSNP,]
  
  ### Get stats based on covmat of all the data
  dat.out.defaultS <- Getdf(dat2, cols)
  colnames(dat.out.defaultS)[(ncol(dat.out.defaultS)-5):
                               ncol(dat.out.defaultS)] <-
    c("pcs.default","Hcd.default", "Md.default", 
      "Hd.default", "Kd.default", "Nd.default")
  #head(dat.out.defaultS)
  #dim(dat.out)
  
    # Covariance of neutral loci
    Sneut <- cov(dat2[dat2$s_high==0,cols])
    
    # Covariance of all loci
    Sall <- cov(dat2[,cols])
    
    # Covariance of MCD
    S_mcd <- CovNAMcd(dat2[,cols])
    
    # Compare cov of neutral data to cov_mcd of all data
    neut <- as.numeric(log(abs(Sneut)))
    all <- as.numeric(log(abs(Sall)))
    mcd <- as.numeric(log(abs(S_mcd@cov)))
    par(mar=c(4,1,1,1), bty="n")
     plot(neut, all, pch=19, col="blue",
          xlab="log(abs(covariance)):\nNeutral loci only", ylab="") 
      abline(0,1, lwd=2)
      abline(lm(all~neut), col="blue")
     points(neut, mcd, pch=17, col="orange")
     abline(lm(mcd~neut), col="orange")

    ### Get stats based on mcd covmat
    dat.out <- Getdf(dat.out.defaultS,cols,S=S_mcd@cov)
    head(dat.out)
    colnames(dat.out)[(ncol(dat.out)-3):
                               ncol(dat.out)] <-
    c( "Md.mcd", "Hd.mcd", "Kd.mcd", "Nd.mcd") 
     head(dat.out)


  ### Make plots ####
#    pdf(file = paste(outfile, ".pdf", sep=""),width = 4, height=6)
#      par(mfrow=c(2,1), mar=c(4,4,1,1), oma=c(0,0,2,1))
#      plot(dat.out[,cols[1]], dat.out[,cols[2]], col=factor(dat.out$s_high), xlab="Spearman's rho (GEA)", ylab="XTX (Fst analog)")
#      plot(dat.out[,cols[1]], dat.out[,cols[3]], col=factor(dat.out$s_high), xlab="Spearman's rho (GEA)", ylab="Z-score (LFMM, GEA)")
#      mtext(dat.out$demog[1] , side = 3, outer=TRUE)
#    dev.off()
#    pdf(file = paste(outfile, "mulitD.pdf", sep=""),width = 6, height=8)
#      ComparePlot(dat.out, colorVect=factor(dat$s_high), 9500:9996)
#    dev.off()

  ### Calculate Empirical Power ###
    out <- data.frame(infile = as.character(infile),
             demog = as.character(dat.out$demog[1]),
             xtx.ep = getEmpPower(dat.out$xtx, dat.out$s_high==0),
             bf.ep = getEmpPower(dat.out$log.bf, dat.out$s_high==0),
             rho.ep = getEmpPower(dat.out$rho, dat.out$s_high==0),
             lfmm.ep = getEmpPower(dat.out$TW.Zscore, dat.out$s_high==0),
             Hcd.ep = getEmpPower(dat.out$Hcd, dat.out$s_high==0),
             pcs.ep = getEmpPower(dat.out$pcs, dat.out$s_high==0),
             Md.default.ep = getEmpPower(dat.out$Md.default, dat.out$s_high==0),
             Hd.default.ep = getEmpPower(dat.out$Hd.default, dat.out$s_high==0),
             Kd.default.ep = getEmpPower(dat.out$Kd.default, dat.out$s_high==0),
             Nd.default.ep = getEmpPower(dat.out$Nd.default, dat.out$s_high==0),
             Md.mcd.ep = getEmpPower(dat.out$Md.mcd, dat.out$s_high==0),
             Hd.mcd.ep = getEmpPower(dat.out$Hd.mcd, dat.out$s_high==0),
             Kd.mcd.ep = getEmpPower(dat.out$Kd.mcd, dat.out$s_high==0),
             Nd.mcd.ep = getEmpPower(dat.out$Nd.mcd, dat.out$s_high==0),
    )
  summarydf <- rbind(summarydf,out)
}

write.table(summarydf,file = paste(resultsdir,"LandsharcSummary.txt", sep=""))

colnames(summarydf)<-sub(".ep", "",colnames(summarydf))

pdf(paste(resultsdir,"LandsharcSummary.pdf", sep=""), width = 6, height = 8)
    par(mfrow=c(4,1), mar=c(3,3,1,1),oma=c(1,3,1,0))
    x <- 0.5
    y <- 1.1
    stats <- c(3:12)
    colors <- c(rep("magenta",4), rep("blue",1), rep("lightblue",3), rep("darkblue",3))
    boxplot(summarydf[summarydf$demog=="IM",stats], col=colors,
            ylim=c(0,1.2))
    text(x,y, "IM", cex=2)
    boxplot(summarydf[summarydf$demog=="IBD",stats], col=colors,
            ylim=c(0,1.2))
    text(x,y, "IBD", cex=2)
    boxplot(summarydf[summarydf$demog=="1R",stats], col=colors,
            ylim=c(0,1.2))
    text(x,y, "1R", cex=2)
    boxplot(summarydf[summarydf$demog=="2R",stats], col=colors,
            ylim=c(0,1.2))
    text(x,y, "2R", cex=2)
    mtext("Empirical Power", side=2, outer=TRUE)
dev.off()
# colMeans(summarydf[summarydf$demog=="IM",3:10])
# colMeans(summarydf[summarydf$demog=="IBD",3:10])
# colMeans(summarydf[summarydf$demog=="1R",3:10])
# colMeans(summarydf[summarydf$demog=="2R",3:10])
