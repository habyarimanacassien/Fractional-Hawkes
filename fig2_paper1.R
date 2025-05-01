
#' Monte Carlo Simulation & Integration of the Laplace Transform inversion
#' Plots the results of Monte Carlo Simulation & integral of the Laplace 
#' Transform inversion of FWP at multiple time values.
library(MittagLeffleR)
library(pracma)

set.seed(123)

b <- beta <- 0.5  # Fix a value for $\beta$
l0 <- lambda_0 <- 1  #Fix a value for the background intensity rate
a <- alpha <- 0.1 #It can be a vector of values
W <- length(a)
for(w in 1:W){
  g <- gamma <- c(0.1, 0.8, 1.7)
J <- length(g)

X11(width=7.5,height=3.2)  #size of PORTRAIT PAPER
par(mfrow=c(1,3))
for(j in 1:J){
###### Monte Carlo Simulation
  Time <- c(1, 30, 60, 90, 120)
  ysim = NULL
for (i in 1:length(Time)){
  T <- Time[i]
  e <- epsilon <- 1e-10
  NSimPoints2 <- c()
  N <- 1000
  #N <- 2
  for (k in 1:N) {
    SimPoints <- c(0)
    t <- 0
    n <- 0
    while (t<T){
      M <- l0 + sum(a[w]*mlf(-g[j]*(t+e-SimPoints)^b,b,b,1)*g[j]*(t+e-SimPoints)^(b-1))
      E <- rexp(1,M)
      t <- t+E
      U <- runif(1)
      if (U<(((l0+sum(a[w]*mlf(-g[j]*(t-SimPoints)^b,b,b,1)*g[j]*(t-SimPoints)^(b-1))))/M)) 
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

###Analytical Case
A <- l0/(1-a[w])
B <- a[w]*l0/(1-a[w])
C <- (1-a[w])*g[j]
P <- B/C^2
Q <- 2*B/C
yan <- Re(A*Time - P*exp(C^2*Time)*erfc(C*sqrt(Time)) - Q*sqrt(Time/pi) + P)
print(yan)
ymax <- max(ysim,yan)

# Plot
plot(Time,ysim, lty=1, type="p", col = "red", lwd=2, xaxt = "n", ylim=c(0,ymax),xlab="",ylab="")
axis(1, at = c(0, 40, 80, 120)) #X-axis
lines(Time,yan, col = "blue", type="o", lty=1, lwd=2)

# Add a legend
title(xlab="t", ylab="E(N(t))", line=2.3, cex.lab=1.3)
legend(x = "bottomright", inset = c(-0.0, 0), 
       legend = c("Monte Carlo Simulation", "Analytical result"),
       bty = "n", lty = c(0, 1), pch=c(1,-1), col = c("red", "blue"),
       lwd = 2, cex = 0.8, xpd = TRUE) 
}
}
