


#source("/Users/katie/Desktop/CurrResearch/3-MINOTAUR/R/distanceFunctions.R")



#### 2 Refuge Simulation Example ######


  # install rrcovNA
  
  S_mcd <- CovNAMcd(dfv[,colnums])
  S_mcd
  Sneut_mcd <- CovNAMcd(dfv[dfv$s_high==0,colnums])
  Sneut_mcd 
  str(S_mcd@cov)

  plot(log(Sneut), log(S_mcd@cov)); abline(0,1)
  plot(log(Sneut), log(Sneut_mcd@cov)); abline(0,1)
  plot(log(Sneut_mcd@cov), log(S_mcd@cov)); abline(0,1)

 
   quartz()
 
 