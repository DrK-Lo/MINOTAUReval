
if(!("devtools" %in% installed.packages())){install.packages("devtools", dependencies=TRUE)}
if(!("rrcovNA" %in% installed.packages())){install.packages("rrcovNA")}
library(devtools)
install_github("NESCent/MINOTAUR", ref="develop")#, build_vignettes=TRUE)
library(MINOTAUR)
### Create variables
x <- rnorm(1000)
y <- x + rnorm(1000)
var(x)
var(y)
plot(x,y)
# Lets say in this case positive values of x and y are more significant
# regular Md won't capture this

# re-run with these lines of code if you want x and y to contribute equally
#x <- (x-mean(x))/sd(x) 
#y <- (y-mean(y))/sd(y)
#var(y)

#######################
### Choose one (right tailed or two tailed)
### Right-tailed x
  x.p <- (length(x)+1-rank(x))/(length(x)+1)
### Two-tailed X
  xright <- (length(x)+1-rank(x))/(length(x)+1)
  xleft <- rank(x)/(length(x)+1)
  x.p <- xleft
  x.p[x>mean(x)] <- xright[x>mean(x)]
  
### Right tailed Y
  y.p <- (length(y)+1-rank(y))/(length(y)+1)
### Two-tailed Y
  yright <- (length(y)+1-rank(y))/(length(y)+1)
  yleft <- rank(y)/(length(y)+1)
  y.p <- yleft
  y.p[y>mean(y)] <- yright[y>mean(y)]
#######################

plot(x, x.p)
plot(y, y.p)

log.xp <- -log(x.p, 10)
log.yp <- -log(y.p, 10)

#log.xp2 <- log((1-x.p)/x.p)
#hist(log.xp2)

plot(y, log.yp)
plot(x, log.xp)

hist(log.xp) #values of 0 are less significant

### Mahalanobis rank-P
length(x)
(S <- cov(data.frame(log.xp, log.yp))) #modified since I sent to Bob
M1.logp2<-Mahalanobis(data.frame(c(log.xp),c(log.yp)), S = S, M=c(0,0))

hist(M1.logp2)
col.vect <- rep("black", length(x))
pch.vect <- rep(19, length(x))
(M1.cutoff <- quantile(M1.logp2,probs = 0.95))
col.vect[M1.logp2>M1.cutoff] <- "blue"
pch.vect[M1.logp2>M1.cutoff] <- 4
plot(y, M1.logp2, col=col.vect, pch=pch.vect, cex=0.3)
plot(x, M1.logp2, col=col.vect, pch=pch.vect, cex=0.3)

plot(x,y, col=col.vect, pch=19, cex=0.3) # compare to DCMS

### Compare to DCMS 
DCMS <- DCMS(dfp=data.frame(x.p,y.p), dfv=data.frame(x,y))
hist(DCMS)
(D.cutoff <- quantile(DCMS,probs = 0.95))
pch.vect2 <- rep(19, length(x))
col.vect2 <- rep("black", length(x))
col.vect2[DCMS>D.cutoff] <- "blue"
pch.vect2[DCMS>D.cutoff] <- 4
plot(x, DCMS, col=col.vect2, pch=pch.vect, cex=0.3)
plot(y, DCMS, col=col.vect2, pch=pch.vect, cex=0.3)
plot(x, y, col=col.vect2, pch=pch.vect, cex=0.3) # compare to MD

### I'm interested in the difference between the Md and the DCMS plots.
plot(x,y, col=col.vect, pch=19, cex=0.3) # Md
plot(x, y, col=col.vect2, pch=19, cex=0.3) # DCMS

pdf("results/CompareMdRankPandDCMS.pdf", width=6, height=3)
par(mfrow=c(1,2), bty="l", mar=c(4,4,2,0.5))
  plot(log.xp, log.yp, col=col.vect, pch=pch.vect, cex=0.3, 
       xlab="-Log(p) x", ylab="-Log(p) y", main = "Md-rank-P") # Md
  plot(log.xp, log.yp, col=col.vect2, pch=pch.vect2, cex=0.3, 
       xlab="-Log(p) x", ylab="-Log(p) y", main = "DCMS") # DCMS
dev.off()
### By Md, there are more outliers along the outside of the distribution
### Also by Md, there are more outliers along the x variable than the y variable
### when the variance is higher in y than in x. It is important to standardize each variable
### by subtracting the mean and dividing by the variance before calculating Md
### (can rerun code above without comments to see this)
