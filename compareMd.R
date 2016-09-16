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

### Hack Mahalanobis for distance from 0
length(x)
(S <- cov(data.frame(log.xp, log.yp))) #modified since I sent to Bob
M1.logp2<-Mahalanobis(data.frame(c(log.xp,0,0),c(log.yp,0,0)), S = S, subset=1001:1002)
M1.logp2<-M1.logp2[-c(1001:1002)]

hist(M1.logp2)
col.vect <- rep("black", length(x))
(M1.cutoff <- quantile(M1.logp2,probs = 0.9))
col.vect[M1.logp2>M1.cutoff] <- "red"
plot(y, M1.logp2, col=col.vect, pch=19, cex=0.3)
plot(x, M1.logp2, col=col.vect, pch=19, cex=0.3)
plot(x,y, col=col.vect, pch=19, cex=0.3) # compare to DCMS

### Compare to DCMS 
DCMS <- DCMS(data.frame(x.p,y.p))
hist(DCMS)
(D.cutoff <- quantile(DCMS,probs = 0.95))
col.vect2 <- rep("black", length(x))
col.vect2[DCMS>D.cutoff] <- "red"
plot(x, DCMS, col=col.vect2, pch=19, cex=0.3)
plot(y, DCMS, col=col.vect2, pch=19, cex=0.3)
plot(x, y, col=col.vect2, pch=19, cex=0.3) # compare to MD

### I'm interested in the difference between the Md and the DCMS plots.
plot(x,y, col=col.vect, pch=19, cex=0.3) # Md
plot(x, y, col=col.vect2, pch=19, cex=0.3) # DCMS

plot(log.xp, log.yp, col=col.vect, pch=19, cex=0.3) # Md
plot(log.xp, log.yp, col=col.vect2, pch=19, cex=0.3) # DCMS

### By Md, there are more outliers along the outside of the distribution
### Also by Md, there are more outliers along the x variable than the y variable
### when the variance is higher in y than in x. It is important to standardize each variable
### by subtracting the mean and dividing by the variance before calculating Md
### (can rerun code above without comments to see this)
