
#'  KE Lotterhos
#'  Feb 18, 2016
#'  Get a dataframe with all multivariate stats
#'  @param dfv2 is a dataframe containing observations in rows and statistics in columns
#'  @param colnums is a vector of column numbers on which to compute the multivar stat
#'  @param alpha is for the FastPCS function
#'  @author KE Lotterhos

# dfv <- read.table("data/OneRefSim.txt", header=TRUE)
# colnums <- c(10,11,12,15)
#source("misc/evalsims/distanceFunctionsOther.R")

Getdf <- function(dfv, colnums=1:ncol(dfv), S=NULL, M=NULL, subset=1:nrow(dfv), 
                  whichstats=c("md", "hd", "kd", "nd")){
 
  ### Check for duplicated rows and abort
  if (any(duplicated(dfv))) {
    writeLines("Error: Your data frame has duplicated rows")
    dfv[duplicated(dfv),]
    break()
  }
  
  dfv2 <- dfv[,colnums]
  rows.keep <- rep(TRUE, nrow(dfv2))
  ### Remove NAs
  if (any(is.na(dfv2))) {
    rows.keep <- !is.na(rowSums(dfv2))
    dfv2 <- dfv2[rows.keep,]
    writeLines(c("Rows with NAs were removed.  The data now has this many rows:", nrow(dfv2)))
  }

  if("md" %in% whichstats){
  writeLines("Calculating outlierliness based on Mahalanobis distance...")
    x<- system.time({
      tx <- try(Md <- Mahalanobis(dfv2, S=S, M=M))
      if("try-error" %in% class(tx)){
        Md <- NA}
      })
    print(x)
  }else{
    Md= "NotCalc"
  }

  if("hd" %in% whichstats){
  writeLines("Calculating outlierliness based on harmonic mean of euclidean distance correcting for covariances...")
    x<- system.time({
      tx <- try(Hd <- harmonicDist(dfv2,S=S, subset=subset))
      if("try-error" %in% class(tx)){
        Hd <- NA}
      })
    print(x)
  }else{
    Hd = "NotCalc"
  }
  
  if("kd" %in% whichstats){
  writeLines("Calculating outlierliness based on kernel density and given bandwith (assume covar)...")
    x<- system.time({
      
      #plot(bw, Kd.ML)
     tx <- try({
      Kd.ML <- kernelDeviance(dfv2, bandwidth = bw, S=S, subset=subset)
      #plot(bw, log(Kd.ML,10))
      #print(cbind(bw, log(Kd.ML,10)))
      bw.best <- bw[which(Kd.ML==min(Kd.ML))[1]]
      print(bw.best)
      Kd <- kernelDist(dfv2, bandwidth = bw.best, subset=subset)
      })
     if("try-error" %in% class(tx)){Kd <- NA}
    })
    print(x)
  }else{
    Kd = "NotCalc"
  }
  
  if("nd" %in% whichstats){
    writeLines("Calculating outlierliness based on euclidean distance to nearest neighbor (cov)...")
    x<- system.time({
      tx <- try(Nd <- neighborDist(dfv2, S=S, subset=subset))
      if("try-error" %in% class(tx)){
        Nd <- NA}
      })
    print(x)
  }else{
    Nd = "NotCalc"
  }

  dfv$pcs[rows.keep] <- pcs
  dfv$Hcd[rows.keep] <- Hcd
  dfv$Md[rows.keep] <- Md
  dfv$Hd[rows.keep] <- Hd
  dfv$Kd[rows.keep] <- Kd
  dfv$Nd[rows.keep] <- Nd
  return(dfv)
}
