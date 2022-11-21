#Numerical inversions of the inverse Laplace Transform

require(invLT)
require(NORMT3)
library(MittagLeffleR)
library(pracma)
set.seed(123)

#Inversion formula
ivLT.plot <- function(L.FUN, METHOD = iv.opC, tPts = seq(0, 3600, 0.1), nterms = 100, pch = 19, ...)
  {
  FUN <- match.fun(L.FUN)
  METHOD <- match.fun(METHOD)
  iv.LT <- rep(0, length(tPts))
  for(tn in 1:length(tPts)){iv.LT[tn] <- METHOD(FUN, tPts[tn], nterms = nterms)}
  plot(tPts, Re(iv.LT), type = "p", xlab = "", ylab = "", xaxt = "n", ylim = c(lambda_0, ymax), ...) 
  axis(1, at = c(0, 900, 1800, 2700, 3600)) #X-axis
  title(xlab = "t", ylab=expression(lambda (t)), line=2.0, cex.lab=1.2)
  }

beta <- 0.9
T <- 3600  # Time (1 hour = 3600 seconds)
l0<- lambda_0 <- 1.0 # fix the background intensity rate
# Numerical inverse Laplace transform for fractional Hawkes process
a <- alpha <- 0.1 # It can be a vector of values
J <- length(a)
for(j in 1:J){
  X11(width=7.5,height=3.8)  #size of PORTRAIT PAPER

g <- gamma <- c(0.1, 0.8, 1.7)
n <- length(g) 
par(mfrow=c(1,3))
for(i in 1:n){
ymax <- lambda_0/(1-a[j])

# # Analytical function
# t <- seq(from=0, to=T, by=0.1)
# erfc_fxn <- erfc((1-a[j])*g[i]*sqrt(t))
# exp_fxn <- exp((1-a[j])^2*g[i]^2*t)
# z <- t[length(erfc_fxn[Re(erfc_fxn)>0])]
# lambda <- function(t){
#   f <- NULL  
#   f[t<=z] <- l0/(1-a[j])-a[j]*l0/(1-a[j])*exp((1-a[j])^2*g[i]^2*t[t<=z])*erfc((1-a[j])*g[i]*sqrt(t[t<=z]))
#   f[t>z] <- l0/(1-a[j])-a[j]*l0/(1-a[j])*( -1/((a[j]-1)*g[i]*sqrt(pi))*sqrt(1/t[t>z])+ 1/(2*(a[j]-1)^3*g[i]^3*sqrt(pi))*(1/t[t>z])^(3/2)- 3/(4*(a[j]-1)^5*g[i]^5*sqrt(pi))*(1/t[t>z])^(5/2)+
#                                              +15/(8*(a[j]-1)^7*g[i]^7*sqrt(pi))*(1/t[t>z])^(7/2)- 105/(16*(a[j]-1)^9*g[i]^9*sqrt(pi))*(1/t[t>z])^(9/2)+
#                                              +945/(32*(a[j]-1)^11*g[i]^11*sqrt(pi))*(1/t[t>z])^(11/2))
#   f
# } 


# General formula
lambda_gen <- (lambda_0/(1-a[j]))*(1-a[j]*mlf((a[j]-1)*g[i]*t^beta,beta,1,1))


# Intensity function Laplace
L.Hawkes <- function(s){lambda_0*(g[i]+s^beta)/(s*((1-a[j])*g[i]+s^beta))}

# Plot of Numerical inversion of Laplace transform
ivLT.plot(L.Hawkes, METHOD = iv.opC, tPts = seq(0, 3600, 0.1), nterms = 100,col = c("blue"))
axis(1, at = c(0, 900, 1800, 2700, 3600)) #X-axis

# # Plot of analytical result
# lines(t, lambda(t), col = "yellow", type="l", lty=1, lwd=2)

# Plot of general analytical result
lines(t, lambda_gen, col = "red", type="l", lty=1, lwd=2)

#Plot of asymptotic line
asymptot <- 0*t+ymax
lines(t, asymptot, col = "green", lty=2, lwd=1)

#title(main="Analytic inversion vs Numerical result")
#title(main = paste("Analytic inversion vs Numerical result for ",expression(beta=0.1)))
#Legend
legend(x = "bottomright",
       inset = c(0, 0, 0), 
       legend = c("Analytic inversion", "Numerical result", "Limiting line"),
       bty = "n",
       lty = c(1, -1, 2),
       pch=c(-1, 1, -1),
       col = c("red", "blue","green"),
       lwd = 1,
       cex = 0.9,
       xpd = TRUE)
}
}
#rm(list=ls())
