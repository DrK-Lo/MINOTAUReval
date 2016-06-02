x <- c(rnorm(100), rnorm(10,mean=5))
> y <- c(rnorm(100), rnorm(10,mean=5))
library(rrcovNA)
z<-CovNAMcd(data.frame(x,y))
summary(z)
  # why can't I seem to get the robust distances?
z@best
  # are these the points chosen for the robust calculation?
summary(x)
summary(x[z@best])
summary(y)
summary(y[z@best])

plot(z)
plot(z, which="dd")
  # how were the horizontal and vertical lines calculated?

plot(z, which = "tolEllipsePlot", classic = TRUE)
