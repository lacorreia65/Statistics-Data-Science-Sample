
#############################################################
###
### LUIS ALVARO CORREIA Student No. 1006508566
###
### STA2101 - Applied Statistics I - Prof. Jerry Brunner
###
### Assignment #1 - September 20th
###
#############################################################
library(purrr)

#No. of simulations
N   <- 100000

#Sample Size
n    <-  1000

#Real Value of theta
realtheta <- 0.50

#Matrix of simulations
Y <- array(rep(0,(N*n)),dim=c(N, n)) 
Z1 <- vector()

#Under H0: theta = theta0
theta0 <- 0.5

#Generating MC samples under H0
for (i in 1:N) {
  Y[i,] <- rbernoulli(n, p = realtheta)
  Z1[i] <- (sqrt(n)*(mean(Y[i,])-theta0))/sqrt(theta0*(1-theta0))
}

#Plotting the graphs for Z1 and p-value under H0
par(mfrow=c(2,1))
hist(Z1,prob=TRUE)
ft <- function(x) dnorm(x,mean = 0, sd = 1, log = FALSE)
curve(ft, add=TRUE)
p <- 1 - pnorm(Z1, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) 
hist(p, prob=TRUE, main="Histogram of p-value", xlab=expression(Z1(Y[1],...,Y[n])),xlim=c(0,1))
abline(h= 1)

### (b) NOW Changing the real value of theta to 0.55

#Real Value of theta
realtheta <- 0.55

#Under H0: theta = theta0
theta0 <- 0.50


#Matrix of simulations
Y1 <- array(rep(0,(N*n)),dim=c(N, n)) 

#Generating MC samples under H0
for (i in 1:N) {
  Y1[i,] <- rbernoulli(n, p = realtheta)
  Z1[i] <- (sqrt(n)*(mean(Y1[i,])-theta0))/sqrt(theta0*(1-theta0))
}

#Plotting the graphs for Z1 and p-value under H0
par(mfrow=c(2,1))
hist(Z1,prob=TRUE)
ft <- function(x) dnorm(x,mean = 0, sd = 1, log = FALSE)
curve(ft, add=TRUE)
p <- 1 - pnorm(Z1, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) 
hist(p, prob=TRUE, main="Histogram of p-value", xlab=expression(Z1(Y1[1],...,Y1[n])),xlim=c(0,1))
abline(h= 1)

mean(Z1)
