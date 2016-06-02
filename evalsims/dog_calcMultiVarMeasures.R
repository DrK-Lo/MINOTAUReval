#!/usr/bin/env Rscript --vanilla --default-packages=base,utils,rrcovNA,MINOTAUR

#'  KE Lotterhos
#'  Feb 18, 2016
#'  Create a dataframe with all multivariate stats
#'  @param dfv2 A dataframe containing observations in rows and statistics in columns
#'  @param colnums A vector of column numbers on which to compute the multivar stat
#'  @param covar If TRUE, a covariance matrix will be calculated from the data using
#'  covNAmcd. If FALSE, the default covariance matrix will be calculated.
#'  @author KE Lotterhos

require(rrcovNA, quietly=TRUE, warn.conflicts=FALSE)
require(MINOTAUR, quietly=TRUE, warn.conflicts=FALSE)

args <- commandArgs(TRUE)

calcMultiVarMeasures <- function(dfv, colStart=1, colEnd=ncol(dfv), covar=FALSE){
    ### Check for duplicated rows and abort
    if (any(duplicated(dfv))) {
        writeLines("Error: Your data frame has duplicated rows")
        dfv[duplicated(dfv),]
        break()
    }
    
    colnums <- colStart:colEnd
    dfv2 <- dfv[,colnums]
    rows.keep <- rep(TRUE, nrow(dfv2))
    ### Remove NAs
    if (any(is.na(dfv2))) {
        rows.keep <- !is.na(rowSums(dfv2))
        dfv2 <- dfv2[rows.keep,]
        writeLines(c("Rows with NAs were removed.  The data now has this many rows:", nrow(dfv2)))
    }
    
    ### Sort out inputted covariance option
    if (covar==TRUE) {
        covCalc <- CovNAMcd(dfv2)
        covarMatrix <- covCalc@cov
        writeLines("Using a covariance matrix calculated using CovNAMcd.")
    } 
    else if (covar==FALSE) {
        covarMatrix <- NULL
        writeLines("Using the default covariance matrix calculation.")
    } 
    else {
        stop("Covariance matrix option must be TRUE or FALSE.")
    }
    
#     writeLines("Calculating outlierliness based on FastPCS...")
#     x<- system.time({
#         tx <- try(pcs<-FastPCS.out(dfv2))
#         if("try-error" %in% class(tx)){
#             pcs <- NA
#         }
#     })
#     print(x)
#     
#     writeLines("Calculating outlierliness based on clustering (DmWR)...")
#     x<- system.time({
#         tx <- try(Hcd <- hclust.ranking(dfv2))
#         if("try-error" %in% class(tx)){
#             Hcd <- NA
#         }
#     })
#     print(x)
    
    writeLines("Calculating outlierliness based on Mahalanobis distance...")
    x<- system.time({
        tx <- try(Md <- Mahalanobis(dfv2, S=covarMatrix))
        if("try-error" %in% class(tx)){
            Md <- NA}
    })
    print(x)
    
    writeLines("Calculating outlierliness based on harmonic mean of euclidean distance correcting for variances...")
    x<- system.time({
        tx <- try(Hd <- harmonicDist(dfv2, S=covarMatrix))
        if("try-error" %in% class(tx)){
            Hd <- NA}
    })
    print(x)
    
    writeLines("Calculating outlierliness based on kernel density and given bandwith...")
    x<- system.time({
        bw <- c(seq(0.01,0.1,by=0.01),seq(0.2,1,by=0.1))
        #plot(bw, Kd.ML)
        tx <- try({
            Kd.ML <- kernelDeviance(dfv2, bandwidth = bw, S=covarMatrix)
            bw.best <- bw[which(Kd.ML==min(Kd.ML))[1]]
            Kd <- kernelDist(dfv2, bandwidth = bw.best, S=covarMatrix)
        })
        if("try-error" %in% class(tx)){Kd <- NA}
    })
    print(x)
    
    writeLines("Calculating outlierliness based on euclidean distance to nearest neighbor...")
    x<- system.time({
        tx <- try(Nd <- neighborDist(dfv2, S=covarMatrix))
        if("try-error" %in% class(tx)){
            Nd <- NA}
    })
    print(x)
    
#     dfv$pcs[rows.keep] <- pcs
#     dfv$Hcd[rows.keep] <- Hcd
    dfv$Md[rows.keep] <- Md
    dfv$Hd[rows.keep] <- Hd
    dfv$Kd[rows.keep] <- Kd
    dfv$Nd[rows.keep] <- Nd
    return(dfv)
}

foo <- read.table(args[1], sep=",", header=TRUE)

outDF <- calcMultiVarMeasures(foo, colStart=args[2], colEnd=args[3], covar=args[4])

write.table(outDF, file=args[5], sep=",", col.names=TRUE, row.names=FALSE)