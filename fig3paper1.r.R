
#' Monte Carlo Simulation & Integrate the Laplace Transform inversion
#' Plots the results of Monte Carlo Simulation & integral of the Laplace Transform inversion of HWP at multiple time values.

require(invLT)
require(NORMT3)

library(MittagLeffleR)
library(pracma)

set.seed(123)

lambda_0 <- 1 # Fix a value for the background intensity rate
beta <- 0.99  # Fix a value for $\beta$
deltax <- 0.01

alpha <- 0.1 # It can be a vector of values

W <- length(alpha)
for(w in 1:W){
  gamma <- c(0.1, 0.8) 
  G <- length(gamma)

  X11(width=7.5,height=4.5)  #size of PORTRAIT PAPER
  par(mfrow=c(1,2))
  for(j in 1:G){
    
    Time <- c(1, 30, 60, 90, 120)
ysim = NULL
###### Monte Carlo Simulation

for (i in 1:length(Time)){
  T <- Time[i]
  epsilon <- 1e-10
  NSimPoints2 <- c()
  N <- 1000
  for (k in 1:N) {
    SimPoints <- c(0)
    t <- 0
    n <- 0
    while (t<T){
      M <- lambda_0 + sum(alpha[w]*mlf(-gamma[j]*(t+epsilon-SimPoints)^beta,beta,beta,1)*gamma[j]*(t+epsilon-SimPoints)^(beta-1))
      E <- rexp(1,M)
      t <- t+E
      U <- runif(1)
      if (U<(((lambda_0+sum(alpha[w]*mlf(-gamma[j]*(t-SimPoints)^beta,beta,beta,1)*gamma[j]*(t-SimPoints)^(beta-1))))/M))
      {
        n <- n+1
        SimPoints <- c(SimPoints,t)
      }
    }
    if (SimPoints[length(SimPoints)]>T) {n <- n-1}
    NSimPoints2 <-  c(n,NSimPoints2)
  }
  ysim[i] <- mean(NSimPoints2)
}
ysim <- ysim
print(ysim)

###===== Numerical Integration

ivLT.value <- function(L.FUN, METHOD, tPts, nterms){
  FUN <- match.fun(L.FUN)
  METHOD <- match.fun(METHOD)
  iv.LT <- rep(0, length(tPts))
  for(tn in 1:length(tPts)){iv.LT[tn] <- METHOD(FUN, tPts[tn], nterms = nterms)}
  return(Re(iv.LT))
  }

ynum = NULL
L.Hawkes <- function(s){lambda_0*(gamma[j]+s^beta)/(s*((1-alpha[w])*gamma[j]+s^beta))}
for (z in 1:length(Time)){
  Z <- Time[z]
  yiv <- ivLT.value(L.Hawkes, METHOD = iv.opChalf, tPts = seq(0, Z, deltax), nterms = 100)
  time <- seq(0, Z, deltax)
  ynum[z] <- sum(yiv[2:length(time)])*deltax
}
ynum <- ynum
print(ynum)
ymax <- max(ysim,ynum)

# Plot
plot(Time,ysim, type="o", col = "red", lty=1, lwd=2, xaxt = "n", ylim=c(0,ymax),xlab="",ylab="")
axis(1, at = c(0, 40, 80, 120)) #X-axis
lines(Time,ynum, col = "blue", type="o", lty=1, lwd=2)

# Add a legend
title(xlab="t", ylab="E(N(t))", line=2.3, cex.lab=1.3)

legend(x = "bottomright",
       inset = c(0, 0),
       legend = c("Monte Carlo Simulation", "Numerical inversion"),
       bty = "n",
       lty = c(1, 1),
       pch=c(1,1),
       col = c("red", "blue"),
       lwd = 2,
       cex = 0.8,
       xpd = TRUE) 
  }
}
