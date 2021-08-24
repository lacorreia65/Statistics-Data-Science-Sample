
#############################################################
###
### LUIS ALVARO CORREIA Student No. 1006508566
###
### STA2101 - Applied Statistics I - Prof. Jerry Brunner
###
### Assignment #2 - September 26th
###
#############################################################

### Ex. 17 ###

SigLvl <- 0.99                      # Set significance level at 99%
zAlpha <- qnorm(1-(1-SigLvl)/2)     # Significance of normalized confidence interval 
n <- 10^5                           # Number of simulations (size of random sample)

set.seed(100650)                    # Set Seed to repsoduce results

h <- function(x){exp(cos(1/x))}     # Define h(x) to be our function of interest and consider f(x) = 1

a <- 0.0                            # Lower limit of integration
b <- 0.5                            # Upper Limit of integration

I <- integrate(h,a,b)               # Integrate funcion *** FAIL ***

u <- runif(n,a,b)                   # Generate random sample of U, where U ~ Uniform[a,b]

x <- h(u)*(b-a)                     # calculates  h(U), on desired interval

curve(h,from=a, to=b, xlab="Function h(x)",ylab="",lwd=2)  # Plot the function of intetrest - illustrative

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
