
#############################################################
###
### LUIS ALVARO CORREIA Student No. 1006508566
###
### STA2101 - Applied Statistics I - Prof. Jerry Brunner
###
### Assignment #1 - September 30th
###
#############################################################

### Ex. 1 ###

library (purrr)

SigLvl <- 0.95                      # Set significance level at 99%
zAlpha <- qnorm(1-(1-SigLvl)/2)     # Significance of normalized confidence interval 
N <- 10^3

n <- 10^3                           # Number of simulations (size of random sample)

set.seed(100650)                    # Set Seed to repsoduce results

Zn <- function(n, theta){
  x <- rbernoulli(n, theta)
  z <- 2*sqrt(n)*(asin(sqrt(mean(x))-asin(sqrt(theta))))
  return(z)
}

RealTheta <- 0.5

Z <- Zn(1:n,RealTheta)

lambda <- 3
X <- rpois(n, lambda)

Y <- (X - lambda)/sqrt(lambda)

var(X); var(Y)

hist(X); hist(Y)

plot(Z, xlab="Function Zn",ylab="",type="l", lwd=1)  # Plot the function of intetrest - illustrative

Mn <- cumsum(x)/(1:n)               # Calculate a series of estimation of Mn, which converges to the value of integral for large n (LLN)

Vn = cumsum((x-Mn)^2)/((1:n)^2)     # Calculate the variance of Mn estimations

I_hat <- mean(Mn)                   # Estimate of Integral Value

UppLim <- Mn+zAlpha*sqrt(Vn/n)      # Calculates the series of Upper Limites for Confidence Interval at given significance
LowLim <- Mn-zAlpha*sqrt(Vn/n)      # Calculates the series of Lower Limites for Confidence Interval at given significance

### Plot the Convergence of MC Simulation

plot(Mn, xlab="Convergence of Estimates",type="l",lwd=1,ylab="MC Estimate")
abline(h=I_hat, untf = FALSE, lty = 3)

Y1 <- min(Mn[(n/2):n])
Y2 <- max(Mn[(n/2):n])

plot(Mn[(n/2):n], ylim=c(min(I_hat,Y1),max(I_hat,Y2)), xlab="ZOOM - Estimates and Confidence Band",type="l",lwd=2,ylab="MC Estimate")
lines(UppLim[(n/2):n],col="red",lwd=1)
lines(LowLim[(n/2):n],col="blue",lwd=1)
abline(h=I_hat, untf = FALSE, lty = 3)

### OUTPUT of Results

sprintf("=========== RESULTS ===========")
sprintf("No. of MC Simulations = %5.0d", n)
sprintf("#1 - Integral Estimation = %2.5f", mean(Mn))
sprintf("#2 - Upper Limit for C.I. @%3.1f%% significance = %2.5f",SigLvl*100.0, mean(UppLim))
sprintf("#2 - Lower Limit for C.I. @%3.1f%% significance = %2.5f",SigLvl*100.0, mean(LowLim))

####################################


# Ex. 8
library(MASS)
# library(plot3D)
library(graphics)

mu <- cbind(c(7,0,6))
Sigma <- rbind( c(1,0,0), 
               c(1,2,0),
               c(0,-1,1) )
A  <-  rbind(c(1,1,0),
          c(0,1,1) ); A
muY <- A %*% mu # E(Y)
SigmaY <- A %*% Sigma %*% t(A) # cov(Y)

n = 10^2

Y <- mvrnorm(n, rep(0, 2), SigmaY)

# Generate sample from N(mu, Sigma)
head(Y)                                      
# Calculate kernel density estimate
Y.kde <- kde2d(Y[,1], Y[,2], 50)   # from MASS package

image(Y.kde)       # from base graphics package
contour(Y.kde, add = TRUE)  


#example(persp3D)  
#example(surf3D)  
#example(slice3D)

contour(Y.kde)
image(Y.kde)
persp(Y.kde, phi = 45, theta = 30)

# fancy contour with image
image(Y.kde); contour(Y.kde, add = T)

# fancy perspective
persp(Y.kde, phi = 45, theta = 30, shade = .1, border = NA)

##################
n <- 10^3
RealTheta <- 3

Theta <- RealTheta

RealMean <- 3.5
RealSD <- 1.5

X <- rexp(n, RealTheta)

# X <- rnorm(n, mean=RealMean, sd=RealSD)


hist(X,probability=TRUE)

A <- array(rep(-1/n,(n-1)*n),dim=c(n,n))           # Set A matrix
A[n,] <- rep(1/n,n);
for (i in 1:(n-1)) A[i,i] <- A[i,i]+1.0            # fill Diagonal

Y <- A %*% X
hist(Y)


f  <-  function(y) {
  theta1 <- mean(y)
  theta2 <- var(y)
  return(1/(sqrt(theta2)*sqrt(2*pi))*exp(-0.5*(y-theta1)^2/theta2))   # Univariate Normal
}

hist(Y,prob=T)
curve(f,add=T)



# EX <- RealTheta
# VarX <- RealTheta^2 


g <- function (x) { 1/ exp(x)}

Ybar <- mean(Y)

DPYbar <- as.double(sqrt(var(Y)))

Y1 <- sqrt(n)*(g(Y)-g(Ybar))/DPYbar

hist(Y1)
plot(Y1)

#################

library(mvtnorm)
x.points <- Y[,1]
y.points <- Y[,2]
z <- matrix(0,nrow=n,ncol=n)
#mu <- c(1,1)
#sigma <- matrix(c(2,1,1,1),nrow=2)
for (i in 1:n) {
  for (j in 1:n) {
    z[i,j] <- c(x.points[i],y.points[j])
  }
}
contour(x.points,y.points,z)

##################

var(mvrnorm(n = 1000, rep(0, 2), SigmaY))
var(mvrnorm(n = 1000, rep(0, 2), Sigma, empirical = TRUE))
