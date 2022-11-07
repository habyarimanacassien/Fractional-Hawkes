# Thinning simulation for (fractional) Hawkes processes
# This program estimates the distribution of number of events up to T

#################################################
### (1) Fractional Hawkes Process vs. Poisson Process

set.seed(1)

#install.packages("MittagLeffleR")
#install.packages("pracma")
#install.packages("arm")
library(MittagLeffleR)
library(pracma)
library(arm)

is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

T <- c(1,5,10)
Niter <- 10000 #Number of iterations
mu <- 1

x11(width=20.5,height=13.5) #size of PORTRAIT PAPER
par(mfrow=c(3,2))
for(s in 1:length(T))
{
m <- 0.5 # It can be a vector of values
for(i in 1:length(m))
{
  epsilon <- 1e-10
beta <- c(0.5,0.9)

for(j in 1:length(beta))
{
Niter <- Niter
NSimPoints2 <- c()
for (k in 1:Niter) {
  
  SimPoints <- c(0)
  t <- 0
  n <- 0
  
  while (t<T[s]){
    M <- mu + sum(m[i]*mlf(-(t+epsilon-SimPoints)^beta[j],beta[j],beta[j],1)*(t+epsilon-SimPoints)^(beta[j]-1))
    E <- rexp(1,M)
    t <- t+E
    U <- runif(1)
    if (U<(((mu+sum(m[i]*mlf(-(t-SimPoints)^beta[j],beta[j],beta[j],1)*(t-SimPoints)^(beta[j]-1))))/M)) {
      n <- n+1
      SimPoints <- c(SimPoints,t)		
    }
  }
  if (SimPoints[length(SimPoints)]>T[s]) {n <- n-1}
  NSimPoints2 <-  c(n,NSimPoints2)
  
}
x <- 0:max(NSimPoints2)
yp <- dpois(x,mu*T[s])
ymax <- max(yp)
mode2 <- length(which(Mode(NSimPoints2) == NSimPoints2))
y2 <- mode2/length(NSimPoints2)
ymax <- max(yp,y2)*1.05

discrete.histogram(NSimPoints2,ylim=c(0,ymax),xlab="")
points(x,dpois(x,mu*T[s]),pch=19,col="red",lwd = 3)

title(main="",xlab="Number of events")
# Add a legend
legend(x = "topright",
       inset = c(0, 0), 
       legend = c("Fractional Hawkes process", "Poisson process"),
       bty = "n",
       lty = c(1,0),
       pch=c(-1,1),
       col = c("blue", "red"),
       lwd = 6,
       cex = 0.8,
       xpd = TRUE)
}
}
}