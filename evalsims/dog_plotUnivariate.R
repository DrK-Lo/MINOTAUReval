#!/usr/bin/env Rscript --vanilla --default-packages=base,utils

## plotUnivariate.R
##
## A script to plot univariate measures for each dog breed in the Schlamp (2016) dataset.
## All input files should be in CSV format with appropriate headers.
##
## Rscript plotUnivariate.R <input_df> <snp_pos> <region_pos> <plot_out>

source("~/Desktop/Schlamp_dog_scans/flex_manhattan.R")

plotUnivariate <- function(input_df, snp_pos, region_pos, plot_out) {
    
    input_data <- read.table(input_df, sep=",", header=TRUE)
    
    snp <- read.table(snp_pos, sep=",", header=TRUE)
    
    region <- read.table(region_pos, sep=",", header=TRUE)
    
    png(plot_out, width=8.5, height=11, units="in", res=300)
    
    par(mfrow=c(5,1), mar=c(2,4,1,2), oma=c(1,1,1,1))
    
    flex_manhattan(input_data, chr="chr", bp="bp", stat="hapFLK", snp="rs", 
                   outThresh1=0.05, outThresh2=0.01, high_reg=region, high_snp_pos=snp)
    
    flex_manhattan(input_data, chr="chr", bp="bp", stat="unstd_iHS", snp="rs", 
                   outThresh1=0.05, outThresh2=0.01, high_reg=region, high_snp_pos=snp)
    
    flex_manhattan(input_data, chr="chr", bp="bp", stat="win51_pi", snp="rs", 
                   outThresh1=0.05, outThresh2=0.01, high_reg=region, high_snp_pos=snp)
    
    flex_manhattan(input_data, chr="chr", bp="bp", stat="win51_TajD", snp="rs", 
                   outThresh1=0.05, outThresh2=0.01, high_reg=region, high_snp_pos=snp)
    
    flex_manhattan(input_data, chr="chr", bp="bp", stat="win51_H12", snp="rs", 
                   outThresh1=0.05, outThresh2=0.01, high_reg=region, high_snp_pos=snp)
    
    dev.off()
}

plotUnivariate(commandArgs(trailingOnly = TRUE)[1], commandArgs(trailingOnly = TRUE)[2], commandArgs(trailingOnly = TRUE)[3], commandArgs(trailingOnly = TRUE)[4])