# Thinning simulation for (fractional) Hawkes processes
# This program estimates the distribution of number of events up to T

### Fractional Hawkes Process vs. Poisson Process

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

#X11(width=7.5,height=2.7)  #size of PORTRAIT PAPER
x11(width=20.5,height=13.5) 
par(mfrow=c(3,2))
for(s in 1:length(T))
{
m <- 0.01 #c(0.01,0.05,0.10) # For stability
for(i in 1:length(m))
{
  epsilon <- 1e-10
beta <- c(0.5,0.9)# c(0.1,0.5,0.9)

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

#x11(width=20.5,height=11.5) 
discrete.histogram(NSimPoints2,ylim=c(0,ymax),xlab="")
points(x,dpois(x,mu*T[s]),pch=19,col="red",lwd = 3)

#Add title
#title(main=paste("T",paste("=",T[s]),"&","a",paste("=",m[i]),"&","b",paste("=",beta[j])),xlab="Number of events")#, ylab="Probability", line=1.7, cex.lab=1.0)
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
#################################################
### Fractional Hawkes Process vs. Hawkes Process with Exponential decay

set.seed(123)

T <- T
Niter <- Niter
mu <- mu

# X11(width=7.5,height=2.7)  #size of PORTRAIT PAPER
x11(width=20.5,height=12.5) 
par(mfrow=c(3,2))
for(s in 1:length(T))
{
  beta <- 0.99 #c(0.9,0.95,1.00) # For stability
  for(j in 1:length(beta))
  {
    epsilon <- 1e-10
    m <- c(0.1,0.5)# c(0.1,0.5,0.9)

    for(i in 1:length(m))
    {

NSimPoints1 <- c()
NSimPoints2 <- c()

for (k in 1:Niter) {
  
  SimPoints_1 <- c(0)
  t1 <- 0
  n1 <- 0
  
  while (t1<T[s]){
    M1 <- mu + sum(m[i]*exp(-1*(t1+epsilon-SimPoints_1)))
    E1 <- rexp(1,M1)
    t1 <- t1+E1
    U <- runif(1)
    if (U<((mu+sum(m[i]*exp(-1*(t1-SimPoints_1))))/M1)) {
      n1 <- n1+1
      SimPoints_1 <- c(SimPoints_1,t1)		
    }
  }
  
  if (SimPoints_1[length(SimPoints_1)]>T[s]) {n1 <- n1-1}
  
  NSimPoints1 <-  c(n1,NSimPoints1)
  
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
x1 <- 0:max(NSimPoints1)
x2 <- 0:max(NSimPoints2)
xmax <- max(x1,x2)+1
mode1 <- length(which(Mode(NSimPoints1) == NSimPoints1))
y1 <- mode1/length(NSimPoints1)
mode2 <- length(which(Mode(NSimPoints2) == NSimPoints2))
y2 <- mode2/length(NSimPoints2)
ymax <- max(y1,y2)*1.05
########
y <- sort(NSimPoints2)
v <- seq(min(y),max(y),by=1)
N2 <- c()
for(k in 1:length(v))
{
  N2[k] <- length(which(v[k] == NSimPoints2))/length(NSimPoints2)
}
N2
########
#x11(width=20.5,height=11.5) 
discrete.hist(NSimPoints1,freq=FALSE,pch=19,prob.col="blue", xlim=c(-0.5,xmax), ylim=c(0,ymax),xlab="")
points(v,N2,pch=19,col="red",lwd = 3)

#Add title
#title(main=paste("T",paste("=",T[s]),"&","b",paste("=",beta[j]),"&","a",paste("=",m[i])),xlab="Number of events")#, ylab="Probability", line=1.3, cex.lab=1.5)
title(main="",xlab="Number of events")
# Add a legend
legend(x = "topright",
       inset = c(0, 0), 
       legend = c("Fractional Hawkes process", "Exponential Hawkes process"),
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