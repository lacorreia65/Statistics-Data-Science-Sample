#---
#title: "CHL5223 - Applied Bayesian Methods - Assignment 1"
#author: "Luis Correia - Student No. 1006508566"
#date: "January 18th 2020"
#---

library(ggplot2)
library(tidyverse)
#--------------------------------------------------------------------------------------
# CASE 1 - Non-Informative Prior
#--------------------------------------------------------------------------------------
alpha1 <- beta1 <- 1.0
set.seed(125)
hist(rbeta(10000,1, 1), # main=expression("Histogram of"~ list(theta[i]) ~ "|" ~ list(alpha[i],beta[i])),
     xlab = expression(list(theta[i]) ~ "|" ~ list(alpha[i],beta[i])))
#--------------------------------------------------------------------------------------
# CASE 2 - Informative Prior
#--------------------------------------------------------------------------------------
# --- Bisection Algorithm to obtain Alpha and Beta that safisfies the Informative Prior
#--------------------------------------------------------------------------------------

epsilon <- 10^(-6)  # Precision desired for the estimates of alpha and beta

beg <- 0.000001   # Initial Value for Alpha 
end <- 60.0       # Final value for Alpha
len <- 100        # Number of sections in the interval
Err <- FALSE

while(end - beg > epsilon) {
  aa<- seq(beg, end,length=len)
  bb<- aa
  tb <- cbind(aa, bb, cdfp95pct=pbeta(.60, aa,bb)-pbeta(.40,aa,bb))
  newbeg <- tail(which(tb[,3]<0.95),1)
  newend <- head(which(tb[,3]>0.95),1)
  if (newbeg == len) {
    # Bad choice of Beg and End - Stop condition
    cat("Input range to estimate Alpha is invalid.(Bisection Algorithm aborted!)")
    beg <- end
    Err <- TRUE
  }
  else {
    beg <- tb[newbeg,1]
    end <- tb[newend,1]
  }  
}

if (!Err) {
  alpha2 <- beta2 <- (beg+end)/2;
  names(alpha2) <- "Alpha"
  names(beta2) <- "Beta"
  cat("The prior distribution of Theta is Beta(", alpha2, ",", beta2,")\n")
  cat("This distribution guarantees that 95% of its mass is between 0.40 and 0.60 as we can see below:\n")
  cat("-> pbeta(0.4, 47.29982, 47.29982)=",pbeta(0.4, alpha2, beta2),"and ")
  cat("pbeta(0.6, 47.29982, 47.29982)=",pbeta(0.6, alpha2, beta2))
  cat("\n-> qbeta(c(0.025,0.975), 47.29982, 47.29982)=",qbeta(c(0.025,0.975), 47.29982, 47.29982))
  
}  
set.seed(1245)
hist(rbeta(10000,47.29982, 47.29982), # main=expression("Histogram of"~ list(theta[i]) ~ "|" ~ list(alpha[i],beta[i])),
     xlab = expression(list(theta[i]) ~ "|" ~ list(alpha[i],beta[i])))

# Setup Variables
N1 <- N2 <- N3 <- 5001    # No. of voters per district
ThrsVictoryD1 <- trunc(N1/2,0)+1  # Threshold for victory on district #1
ThrsVictoryD2 <- trunc(N2/2,0)+1  # Threshold for victory on district #2
ThrsVictoryD3 <- trunc(N3/2,0)+1  # Threshold for victory on district #3

# Sample for Vote Intention on each district
SD1 <- data.frame(Purple = 53, Brown = 45); SD1 <- cbind(District = 1, SD1, Total=sum(SD1))
SD2 <- data.frame(Purple = 72, Brown = 78); SD2 <- cbind(District = 2, SD2, Total=sum(SD2))
SD3 <- data.frame(Purple = 18, Brown = 22); SD3 <- cbind(District = 3, SD3, Total=sum(SD3))

kableExtra::kable(rbind(SD1, SD2, SD3),"latex", booktabs = T, caption = "Sampling of Vote Intention by District")

#--------------------------------------------------------------------------------------
# CASE 1 - Non-Informative Prior
#--------------------------------------------------------------------------------------

# District 1 - Calculate Posterior Distribution, Mean, Std Deviation and 95% C.I. 
posD1_Alpha_NI <- alpha1+SD1$Purple; #posD1_Alpha_NI
posD1_Beta_NI <- beta1+SD1$Total-SD1$Purple; #posD1_Beta_NI

cat("\n\n--- Summary for District 1 ---\n")

meanD1 <- posD1_Alpha_NI/(posD1_Alpha_NI+posD1_Beta_NI); 
cat("\nPosterior Probability for the percent o voters for Purple party in Dist.#1 is", format(meanD1, digits = 6, nsmall = 3))

cat("\nMean=", format(meanD1, digits = 6, nsmall = 3))

stdevD1 <- sqrt((posD1_Alpha_NI*posD1_Beta_NI)/((posD1_Alpha_NI+posD1_Beta_NI)^2*(posD1_Alpha_NI+posD1_Beta_NI+1))); 
cat("\nStd. Deviation=", format(stdevD1, digits = 6, nsmall = 3))

ci_D1_NI <- qbeta(c(.025, .975),posD1_Alpha_NI,posD1_Beta_NI)
cat("\n95% Credible Region (", format(ci_D1_NI[1], digits = 4, nsmall = 3), ",",format(ci_D1_NI[2], digits = 4, nsmall = 3),")")

ProbD1_NI <- 1-extraDistr::pbbinom(ThrsVictoryD1-1,N1,posD1_Alpha_NI,posD1_Beta_NI) # Beta-Binomial probability
cat("\nPosterior Probability Purple party win in Dist.#1 is", format(ProbD1_NI, digits = 6, nsmall = 3))

#--------------------------------------------------------------------------------------
# CASE 1 - Non-Informative Prior
#--------------------------------------------------------------------------------------

# District 2 - Calculate Posterior Distribution, Mean, Std Deviation and 95% C.I. 
posD2_Alpha_NI <- alpha1+SD2$Purple; #posD2_Alpha_NI
posD2_Beta_NI <- beta1+SD2$Total-SD2$Purple; #posD2_Beta_NI

cat("\n\n--- Summary for District 2 ---\n")

meanD2 <- posD2_Alpha_NI/(posD2_Alpha_NI+posD2_Beta_NI); 
cat("\nPosterior Probability for the percent o voters for Purple party in Dist.#2 is", format(meanD2, digits = 6, nsmall = 3))

cat("\nMean=", format(meanD2, digits = 6, nsmall = 3))

stdevD2 <- sqrt((posD2_Alpha_NI*posD2_Beta_NI)/((posD2_Alpha_NI+posD2_Beta_NI)^2*(posD2_Alpha_NI+posD2_Beta_NI+1))); 
cat("\nStd. Deviation=", format(stdevD2, digits = 6, nsmall = 3))

ci_D2_NI <- qbeta(c(.025, .975),posD2_Alpha_NI, posD2_Beta_NI)
cat("\n95% Credible Region (", format(ci_D2_NI[1], digits = 4, nsmall = 3), ",",format(ci_D2_NI[2], digits = 4, nsmall = 3),")")

ProbD2_NI <- 1-extraDistr::pbbinom(ThrsVictoryD2-1,N2,posD2_Alpha_NI,posD2_Beta_NI)  # Beta-Binomial probability
cat("\nPosterior Probability Purple party win in Dist.#2 is", format(ProbD2_NI, digits = 6, nsmall = 3))

#--------------------------------------------------------------------------------------
# CASE 1 - Non-Informative Prior
#--------------------------------------------------------------------------------------

# District 3 - Calculate Posterior Distribution, Mean, Std Deviation and 95% C.I. 
posD3_Alpha_NI <- alpha1+SD3$Purple; #posD3_Alpha_NI
posD3_Beta_NI <- beta1+SD3$Total-SD3$Purple; #posD3_Beta_NI

cat("\n\n--- Summary for District 3 ---\n")

meanD3 <- posD3_Alpha_NI/(posD3_Alpha_NI+posD3_Beta_NI); 
cat("\nPosterior Probability for the percent o voters for Purple party in Dist.#3 is", format(meanD3, digits = 6, nsmall = 3))

cat("\nMean=", format(meanD3, digits = 6, nsmall = 3))

stdevD3 <- sqrt((posD3_Alpha_NI*posD3_Beta_NI)/((posD3_Alpha_NI+posD3_Beta_NI)^2*(posD3_Alpha_NI+posD3_Beta_NI+1))); 
cat("\nStd. Deviation=", format(stdevD3, digits = 6, nsmall = 3))

ci_D3_NI <- qbeta(c(.025, .975),posD3_Alpha_NI, posD3_Beta_NI)
cat("\n95% Credible Region (", format(ci_D3_NI[1], digits = 4, nsmall = 3), ",",format(ci_D3_NI[2], digits = 4, nsmall = 3),")")

ProbD3_NI <- 1-extraDistr::pbbinom(ThrsVictoryD3-1,N3,posD3_Alpha_NI,posD3_Beta_NI)  # Beta-Binomial probability
cat("\nPosterior Probability Purple party win in Dist.#3 is", format(ProbD3_NI, digits = 6, nsmall = 3))

#--------------------------------------------------------------------------------------
# CASE 2 - Informative Prior
#--------------------------------------------------------------------------------------

# District 1 - Calculate Posterior Distribution, Mean, Std Deviation and 95% C.I. 
posD1_Alpha <- alpha2+SD1$Purple; #posD1_Alpha
posD1_Beta <- beta2+SD1$Total-SD1$Purple; #posD1_Beta

cat("\n\n--- Summary for District 1 ---\n")

meanD1 <- posD1_Alpha/(posD1_Alpha+posD1_Beta); 
cat("\nPosterior Probability for the percent o voters for Purple party in Dist.#1 is", format(meanD1, digits = 6, nsmall = 3))

cat("\nMean=", format(meanD1, digits = 6, nsmall = 3))

stdevD1 <- sqrt((posD1_Alpha*posD1_Beta)/((posD1_Alpha+posD1_Beta)^2*(posD1_Alpha+posD1_Beta+1))); 
cat("\nStd. Deviation=", format(stdevD1, digits = 6, nsmall = 3))

ci_D1 <- qbeta(c(.025, .975),posD1_Alpha,posD1_Beta)
cat("\n95% Credible Region (", format(ci_D1[1], digits = 4, nsmall = 3), ",",format(ci_D1[2], digits = 4, nsmall = 3),")")

ProbD1 <- 1-extraDistr::pbbinom(ThrsVictoryD1-1,N1,posD1_Alpha,posD1_Beta)  # Beta-Binomial probability
cat("\nPosterior Probability Purple party win in Dist.#1 is", format(ProbD1, digits = 6, nsmall = 3))

#--------------------------------------------------------------------------------------
# CASE 2 - Informative Prior
#--------------------------------------------------------------------------------------

# District 2 - Calculate Posterior Distribution, Mean, Std Deviation and 95% C.I. 
posD2_Alpha <- alpha2+SD2$Purple; #posD2_Alpha
posD2_Beta <- beta2+SD2$Total-SD2$Purple; #posD2_Beta

cat("\n\n--- Summary for District 2 ---\n")

meanD2 <- posD2_Alpha/(posD2_Alpha+posD2_Beta); 
cat("\nPosterior Probability for the percent o voters for Purple party in Dist.#2 is", format(meanD2, digits = 6, nsmall = 3))

cat("\nMean=", format(meanD2, digits = 6, nsmall = 3))

stdevD2 <- sqrt((posD2_Alpha*posD2_Beta)/((posD2_Alpha+posD2_Beta)^2*(posD2_Alpha+posD2_Beta+1))); 
cat("\nStd. Deviation=", format(stdevD2, digits = 6, nsmall = 3))

ci_D2 <- qbeta(c(.025, .975),posD2_Alpha, posD2_Beta)
cat("\n95% Credible Region (", format(ci_D2[1], digits = 4, nsmall = 3), ",",format(ci_D2[2], digits = 4, nsmall = 3),")")

ProbD2 <- 1-extraDistr::pbbinom(ThrsVictoryD2-1,N2,posD2_Alpha,posD2_Beta)  # Beta-Binomial probability
cat("\nPosterior Probability Purple party win in Dist.#2 is", format(ProbD2, digits = 6, nsmall = 3))

#--------------------------------------------------------------------------------------
# CASE 2 - Informative Prior
#--------------------------------------------------------------------------------------

# District 3 - Calculate Posterior Distribution, Mean, Std Deviation and 95% C.I. 
posD3_Alpha <- alpha2+SD3$Purple; #posD3_Alpha
posD3_Beta <- beta2+SD3$Total-SD3$Purple; #posD3_Beta

cat("\n\n--- Summary for District 3 ---\n")

meanD3 <- posD3_Alpha/(posD3_Alpha+posD3_Beta); 
cat("\nPosterior Probability for the percent o voters for Purple party in Dist.#3 is", format(meanD3, digits = 6, nsmall = 3))

cat("\nMean=", format(meanD3, digits = 6, nsmall = 3))

stdevD3 <- sqrt((posD3_Alpha*posD3_Beta)/((posD3_Alpha+posD3_Beta)^2*(posD3_Alpha+posD3_Beta+1))); 
cat("\nStd. Deviation=", format(stdevD3, digits = 6, nsmall = 3))

ci_D3 <- qbeta(c(.025, .975),posD3_Alpha, posD3_Beta)
cat("\n95% Credible Region (", format(ci_D3[1], digits = 4, nsmall = 3), ",",format(ci_D3[2], digits = 4, nsmall = 3),")")

ProbD3 <- 1-extraDistr::pbbinom(ThrsVictoryD3-1,N3,posD3_Alpha,posD3_Beta)  # Beta-Binomial probability
cat("\nPosterior Probability Purple party win in Dist.#3 is", format(ProbD3, digits = 6, nsmall = 3))

# Set seed for DEBUG
set.seed(1965)

N=10000                   # No. of Elections to Simulate

# Variable to store X_i (No. of votes for Purple Party on each district)
Sim_X1_NI <- rep(NA,N)
Sim_X2_NI <- rep(NA,N)
Sim_X3_NI <- rep(NA,N)

# Case 1 - Non-Informative Prior
for (n in 1:N) {
  # Generate a Theta_i for each district
  thetaD1 <- rbeta(1,posD1_Alpha_NI,posD1_Beta_NI);
  thetaD2 <- rbeta(1,posD2_Alpha_NI,posD2_Beta_NI) 
  thetaD3 <- rbeta(1,posD3_Alpha_NI,posD3_Beta_NI)
  
  #Calculate The total No. of Votes for Purple Party on each district
  Sim_X1_NI[n] <- rbinom(1,N1,thetaD1)
  Sim_X2_NI[n] <- rbinom(1,N2,thetaD2)
  Sim_X3_NI[n] <- rbinom(1,N3,thetaD3)
}

# Plot the probabilities for each District - Case 1 - Non-Informative Prior

ProbPurpWinD1_NI <- length(which(Sim_X1_NI>(ThrsVictoryD1)))/N; 
ProbPurpWinD2_NI <- length(which(Sim_X2_NI>(ThrsVictoryD2)))/N; 
ProbPurpWinD3_NI <- length(which(Sim_X3_NI>(ThrsVictoryD3)))/N; 

# Calculating probabilities of win 2 districts

ProbPurpWinD1D2_NI <- length(which(Sim_X1_NI>(ThrsVictoryD1) & Sim_X2_NI>(ThrsVictoryD2)))/N; 
ProbPurpWinD1D3_NI <- length(which(Sim_X1_NI>(ThrsVictoryD1) & Sim_X3_NI>(ThrsVictoryD3)))/N; 
ProbPurpWinD2D3_NI <- length(which(Sim_X2_NI>(ThrsVictoryD2) & Sim_X3_NI>(ThrsVictoryD3)))/N; 

# Calculating probabilities of win 3 districts
ProbPurpWinD1D2D3_NI <- length(which(Sim_X1_NI>(ThrsVictoryD1) & Sim_X2_NI>(ThrsVictoryD2) & Sim_X3_NI>(ThrsVictoryD3)))/N;  

SD11 <- cbind(District = 1, Freq=length(which(Sim_X1_NI>(ThrsVictoryD1))), Prop=ProbPurpWinD1_NI)
SD12 <- cbind(District = 2, Freq=length(which(Sim_X2_NI>(ThrsVictoryD2))), Prop=ProbPurpWinD2_NI)
SD13 <- cbind(District = 3, Freq=length(which(Sim_X3_NI>(ThrsVictoryD3))), Prop=ProbPurpWinD3_NI)

kableExtra::kable(rbind(SD11, SD12, SD13),"latex", booktabs = T, caption = "Purple Party wins in 10,000 simulations")

SD21 <- cbind(District = "1 & 2", Freq=length(which(Sim_X1_NI>(ThrsVictoryD1) & Sim_X2_NI>(ThrsVictoryD2))), Prop=ProbPurpWinD1D2_NI)
SD22 <- cbind(District = "1 & 3", Freq=length(which(Sim_X1_NI>(ThrsVictoryD1) & Sim_X3_NI>(ThrsVictoryD3))), Prop=ProbPurpWinD1D3_NI)
SD23 <- cbind(District = "2 & 3", Freq=length(which(Sim_X2_NI>(ThrsVictoryD2) & Sim_X3_NI>(ThrsVictoryD3))), Prop=ProbPurpWinD2D3_NI)

kableExtra::kable(rbind(SD21, SD22, SD23),"latex", booktabs = T, 
                  caption = "Purple Party wins at least 02 districts in 10,000 simulations")

#cat("\nProbability of Purple Win District #1 and #2 =",ProbPurpWinD1D2_NI)
#cat("\nProbability of Purple Win District #1 and #3 =",ProbPurpWinD1D3_NI)
#cat("\nProbability of Purple Win District #2 and #3 =",ProbPurpWinD2D3_NI)

F1 <- as.integer(SD21[1,2])-length(which(Sim_X1_NI>(ThrsVictoryD1) & Sim_X2_NI>(ThrsVictoryD2) & Sim_X3_NI>(ThrsVictoryD3)))
F2 <- as.integer(SD22[1,2])-length(which(Sim_X1_NI>(ThrsVictoryD1) & Sim_X2_NI>(ThrsVictoryD2) & Sim_X3_NI>(ThrsVictoryD3)))
F3 <- as.integer(SD23[1,2])-length(which(Sim_X1_NI>(ThrsVictoryD1) & Sim_X2_NI>(ThrsVictoryD2) & Sim_X3_NI>(ThrsVictoryD3)))
F4 <- length(which(Sim_X1_NI>(ThrsVictoryD1) & Sim_X2_NI>(ThrsVictoryD2) & Sim_X3_NI>(ThrsVictoryD3)))

SD31 <- cbind(District = "Only 1 & 2", Freq=F1, Prop=F1/N)
SD32 <- cbind(District = "Only 1 & 3", Freq=F2, Prop=F2/N)
SD33 <- cbind(District = "Only 2 & 3", Freq=F3, Prop=F3/N)
SD34 <- cbind(District = "All", Freq=F4, Prop=F4/N)

kableExtra::kable(rbind(SD31, SD32, SD33, SD34),"latex", 
                  booktabs = T, caption = "Purple Party wins exactly 02 districts / All districts in 10,000 simulations")

# Plot the probabilities for each District - Case 1 - Non-Informative Prior

cat("\nProbability of Purple Win District #1=",ProbPurpWinD1_NI)
cat("\nProbability of Purple Win District #2=",ProbPurpWinD2_NI)
cat("\nProbability of Purple Win District #3=",ProbPurpWinD3_NI)

# Determining intersections with Full Victory

cat("\n\nProbability that the Purple Party have majority in Town Council=",sum(F1+F2+F3+F4)/N)

# Variable to store X_i (No. of votes for Purple Party on each district)
Sim_X1 <- rep(NA,N)
Sim_X2 <- rep(NA,N)
Sim_X3 <- rep(NA,N)

# Case 1 - Non-Informative Prior
for (n in 1:N) {
  # Generate a Theta_i for each district
  thetaD1 <- rbeta(1,posD1_Alpha,posD1_Beta) 
  thetaD2 <- rbeta(1,posD2_Alpha,posD2_Beta) 
  thetaD3 <- rbeta(1,posD3_Alpha,posD3_Beta)
  
  #Calculate The total No. of Votes for Purple Party on each district
  Sim_X1[n] <- rbinom(1,N1,thetaD1)
  Sim_X2[n] <- rbinom(1,N2,thetaD2)
  Sim_X3[n] <- rbinom(1,N3,thetaD3)
}

# Plot the probabilities for each District - Case 1 - Non-Informative Prior

ProbPurpWinD1 <- length(which(Sim_X1>(ThrsVictoryD1)))/N; 
ProbPurpWinD2 <- length(which(Sim_X2>(ThrsVictoryD2)))/N; 
ProbPurpWinD3 <- length(which(Sim_X3>(ThrsVictoryD3)))/N; 

# Calculating probabilities of win 2 districts

ProbPurpWinD1D2 <- length(which(Sim_X1>(ThrsVictoryD1) & Sim_X2>(ThrsVictoryD2)))/N; 
ProbPurpWinD1D3 <- length(which(Sim_X1>(ThrsVictoryD1) & Sim_X3>(ThrsVictoryD3)))/N; 
ProbPurpWinD2D3 <- length(which(Sim_X2>(ThrsVictoryD2) & Sim_X3>(ThrsVictoryD3)))/N; 

# Calculating probabilities of win 3 districts
ProbPurpWinD1D2D3 <- length(which(Sim_X1>(ThrsVictoryD1) & Sim_X2>(ThrsVictoryD2) & Sim_X3>(ThrsVictoryD3)))/N;  

SD11 <- cbind(District = 1, Freq=length(which(Sim_X1>(ThrsVictoryD1))), Prop=ProbPurpWinD1)
SD12 <- cbind(District = 2, Freq=length(which(Sim_X2>(ThrsVictoryD2))), Prop=ProbPurpWinD2)
SD13 <- cbind(District = 3, Freq=length(which(Sim_X3>(ThrsVictoryD3))), Prop=ProbPurpWinD3)

kableExtra::kable(rbind(SD11, SD12, SD13),"latex", booktabs = T, caption = "Purple Party wins in 10,000 simulations")

SD21 <- cbind(District = "1 & 2", Freq=length(which(Sim_X1>(ThrsVictoryD1) & Sim_X2>(ThrsVictoryD2))), Prop=ProbPurpWinD1D2)
SD22 <- cbind(District = "1 & 3", Freq=length(which(Sim_X1>(ThrsVictoryD1) & Sim_X3>(ThrsVictoryD3))), Prop=ProbPurpWinD1D3)
SD23 <- cbind(District = "2 & 3", Freq=length(which(Sim_X2>(ThrsVictoryD2) & Sim_X3>(ThrsVictoryD3))), Prop=ProbPurpWinD2D3)

kableExtra::kable(rbind(SD21, SD22, SD23),"latex", booktabs = T, 
                  caption = "Purple Party wins at least 02 districts in 10,000 simulations")
                  
F1 <- as.integer(SD21[1,2])-length(which(Sim_X1>(ThrsVictoryD1) & Sim_X2>(ThrsVictoryD2) & Sim_X3>(ThrsVictoryD3)))
F2 <- as.integer(SD22[1,2])-length(which(Sim_X1>(ThrsVictoryD1) & Sim_X2>(ThrsVictoryD2) & Sim_X3>(ThrsVictoryD3)))
F3 <- as.integer(SD23[1,2])-length(which(Sim_X1>(ThrsVictoryD1) & Sim_X2>(ThrsVictoryD2) & Sim_X3>(ThrsVictoryD3)))
F4 <- length(which(Sim_X1>(ThrsVictoryD1) & Sim_X2>(ThrsVictoryD2) & Sim_X3>(ThrsVictoryD3)))

SD31 <- cbind(District = "Only 1 & 2", Freq=F1, Prop=F1/N)
SD32 <- cbind(District = "Only 1 & 3", Freq=F2, Prop=F2/N)
SD33 <- cbind(District = "Only 2 & 3", Freq=F3, Prop=F3/N)
SD34 <- cbind(District = "All", Freq=F4, Prop=F4/N)

kableExtra::kable(rbind(SD31, SD32, SD33, SD34),"latex", booktabs = T, 
                  caption = "Purple Party wins exactly 02 districts / All districts in 10,000 simulations")

# Plot the probabilities for each District - Case 1 - Non-Informative Prior

cat("\nProbability of Purple Win District #1=",ProbPurpWinD1)
cat("\nProbability of Purple Win District #2=",ProbPurpWinD2)
cat("\nProbability of Purple Win District #3=",ProbPurpWinD3)

# Determining intersections with Full Victory

cat("\n\nProbability that the Purple Party have majority in Town Council=",sum(F1+F2+F3+F4)/N)

kableExtra::kable(rbind(cbind(District="1",Simul=ProbPurpWinD1_NI, BetaBin=format(ProbD1_NI, digits = 4, nsmall = 4)), 
                        cbind(District="2",Simul=ProbPurpWinD2_NI, BetaBin=format(ProbD2_NI, digits = 4, nsmall = 4)), 
                        cbind(District="3",Simul=ProbPurpWinD3_NI, BetaBin=format(ProbD3_NI, digits = 4, nsmall = 4))),
                        "latex", booktabs = T, 
                        caption = "Simulated vs. and Beta-Binomial Probability (Case 1: Non-Informative Prior)")

kableExtra::kable(rbind(cbind(District="1",Simul=ProbPurpWinD1, BetaBin=format(ProbD1, digits = 4, nsmall = 4)), 
                        cbind(District="2",Simul=ProbPurpWinD2, BetaBin=format(ProbD2, digits = 4, nsmall = 4)), 
                        cbind(District="3",Simul=ProbPurpWinD3, BetaBin=format(ProbD3, digits = 4, nsmall = 4))),
                        "latex", booktabs = T, 
                        caption = "Simulated vs. and Beta-Binomial Probability (Case 2: Informative Prior)")

# Plot Sample densities to Compare Non-Informative vs. Informative Prior
D_Work <- rbind(data.frame(Distr1=Sim_X1_NI,Distr2=Sim_X2_NI,Distr3=Sim_X3_NI,PType="Non-Informative"),
                data.frame(Distr1=Sim_X1,Distr2=Sim_X2,Distr3=Sim_X3,PType="Informative"))
#colnames(D_Work) <- c("Distr1","Distr2","Distr3","PType" )
D_Work %>% 
  pivot_longer(cols=Distr1:Distr3, names_to = "District", values_to = "Votes" ) %>% 
  ggplot(mapping = aes(x = Votes, group = District))+
  geom_density(aes(colour=District), size=1.2)+  #linetype=District, 
  labs(x = "Votes", y = "Density") +
  facet_grid(~PType)+
  theme_bw()

TbCI <- data.frame(
  District=c(1, 2, 3),
  NonInf=c(format(ci_D1_NI[2]-ci_D1_NI[1], digits = 4, nsmall = 4),
           format(ci_D2_NI[2]-ci_D2_NI[1], digits = 4, nsmall = 4),
           format(ci_D3_NI[2]-ci_D3_NI[1], digits = 4, nsmall = 4)),
  Inform=c(c(format(ci_D1[2]-ci_D1[1], digits = 4, nsmall = 4),
             format(ci_D2[2]-ci_D2[1], digits = 4, nsmall = 4),
             format(ci_D3[2]-ci_D3[1], digits = 4, nsmall = 4))))

kableExtra::kable(TbCI, "latex", booktabs = T, 
                  caption = "Comparison of Amplitude of 95\\% Credible Interval")

# The Data
#------------------------------------------
D1 <- c(53,49,63,72,55,65)
D2 <- c(28,27,36,42,25,35)

# Prior 1 - Calculate the Posterior
#------------------------------------------
# Belief 0 - Tau is known
tau <- 1/36

# Belief 1 - Average height is 66 inches
mu0 <- 66

# Belief 3 - Mu0 is between 63 and 69 with 95% probability
tau0 <- 4/9

# Calculate the Posterior for each data-set
#------------------------------------------
# Data-set 1 - 
cat("-----------------------------------------")
cat("\nSummary Statistics for Data-Set1- (",D1,")\n")
mu_l1 <- (tau0*mu0+length(D1)*tau*mean(D1))/(tau0+length(D1)*tau); 
cat("\nPosterior Mean is ",format(mu_l1, digits = 6, nsmall = 3))

tau_l1 <- (tau0+length(D1)*tau);
cat("\nPosterior Std. Deviation is ",format(sqrt(1/tau_l1), digits = 6, nsmall = 3))

ci_DS1 <- qnorm(c(.025, .975),mean=mu_l1, sd=sqrt(1/tau_l1))
cat("\n95% Credible Region for mean (", format(ci_DS1[1], digits = 4, nsmall = 3), ",",format(ci_DS1[2], digits = 4, nsmall = 3),")")

#xP1D1 <- rnorm(N,mean=mu_l1,sd=sqrt(1/tau_l1))
#cat("\n\nBy simulation (N=",N,") >>> Mean=",mean(xP1D1)," | Std.Dev.=",sd(xP1D1))

# Data-set 2 - 
cat("\n\n-----------------------------------------")
cat("\nSummary Statistics for Data-Set2- (",D2,")\n")
mu_l2 <- (tau0*mu0+length(D2)*tau*mean(D2))/(tau0+length(D2)*tau); 
cat("\nPosterior Mean is ",format(mu_l2, digits = 6, nsmall = 3))

tau_l2 <- (tau0+length(D2)*tau);
cat("\nPosterior Std. Deviation is ",format(sqrt(1/tau_l2), digits = 6, nsmall = 3))

ci_DS2 <- qnorm(c(.025, .975),mean=mu_l2, sd=sqrt(1/tau_l2))
cat("\n95% Credible Region for mean (", format(ci_DS2[1], digits = 4, nsmall = 3), ",",
    format(ci_DS2[2], digits = 4, nsmall = 3),")")

#xP1D2 <- rnorm(N,mean=mu_l2,sd=sqrt(1/tau_l2))
#cat("\n\nBy simulation (N=",N,") >>> Mean=",mean(xP1D2)," | Std.Dev.=",sd(xP1D2))

# Prior 2 - Calculate the Posterior
#------------------------------------------
# Belief 1 - Average height is 66 inches
mu0 <- 66

# Belief 2 - Tau is Gamma(alpha,beta)
thetaP2 <- 4
alphaP2 <- 1
betaP2 <- 36

# Calculate the Posterior for each data-set
#------------------------------------------
# Data-set 1 - 
cat("-----------------------------------------")
cat("\nSummary Statistics for Data-Set1- (",D1,")\n")
alpha_l1P2 <- alphaP2+length(D1)/2
beta_l1P2 <- betaP2+0.5*sum((D1-mean(D1))^2)+(thetaP2*length(D1)*(mean(D1)-mu0)^2)/(2*(thetaP2+length(D1)))

## Analytically calculating the statistics

mu_l1P2 <- (thetaP2*mu0+length(D1)*mean(D1))/(thetaP2+length(D1))
cat("\nPosterior Mean is ",format(mu_l1P2, digits = 6, nsmall = 3))

tau_l1P2A <- 1/((beta_l1P2/((thetaP2+length(D1))*alpha_l1P2))*((2*alpha_l1P2)/(2*alpha_l1P2-2)))
cat("\n(Analytically)Posterior Std. Deviation is ",format(sqrt(1/tau_l1P2A), digits = 6, nsmall = 3)) # Analytically

TauG <- rgamma(10*N,alpha_l1P2, beta_l1P2)
m1 <- rnorm(10*N, mean=mu_l1P2, sd=1/sqrt(TauG*(thetaP2+length(D1))))
tau_l1P2S <- 1/var(m1)
cat("\n(Simulated)Posterior Std. Deviation is ",format(sqrt(1/tau_l1P2S), digits = 6, nsmall = 3))

tau_l1P2 <- tau_l1P2S  # Using the Simulated output to be compatible w/ C.I.

ci_DS1P2 <- qnorm(c(.025, .975),mean=mu_l1P2, sd=sqrt(1/tau_l1P2))
cat("\n95% Credible Region for mean (", format(ci_DS1P2[1], digits = 4, nsmall = 3), 
    ",",format(ci_DS1P2[2], digits = 4, nsmall = 3),")")

ci_DS1P2SD <- quantile(1/sqrt(TauG*(thetaP2+length(D1))), c(.025, .975))
cat("\n95% Credible Region for Std.Deviation (", format(ci_DS1P2SD[1], digits = 4, nsmall = 3), 
    ",",format(ci_DS1P2SD[2], digits = 4, nsmall = 3),")")

#------------------------------------------
# Data-set 2 - 
cat("\n\n-----------------------------------------")
cat("\nSummary Statistics for Data-Set2- (",D2,")\n")
alpha_l2P2 <- alphaP2+length(D2)/2
beta_l2P2 <- betaP2+0.5*sum((D2-mean(D2))^2)+(thetaP2*length(D2)*(mean(D2)-mu0)^2)/(2*(thetaP2+length(D2)))

## Analytically calculating the statistics

mu_l2P2 <- (thetaP2*mu0+length(D2)*mean(D2))/(thetaP2+length(D2))
cat("\nPosterior Mean is ",format(mu_l2P2, digits = 6, nsmall = 3))

#--- NEW
tau_l2P2A <- 1/((beta_l2P2/((thetaP2+length(D2))*alpha_l2P2))*((2*alpha_l2P2)/(2*alpha_l2P2-2)))
cat("\n(Analytically) Posterior Std. Deviation is ",format(sqrt(1/tau_l2P2A), digits = 6, nsmall = 3)) # Analytically

TauG <- rgamma(10*N,alpha_l2P2, beta_l2P2)
m2 <- rnorm(10*N, mean=mu_l2P2, sd=1/sqrt(TauG*(thetaP2+length(D2))))
tau_l2P2S <- 1/var(m2)
cat("\n(Simulated) Posterior Std. Deviation is ",format(sqrt(1/tau_l2P2S), digits = 6, nsmall = 3))

tau_l2P2 <- tau_l2P2S  # Using the Simulated output to be compatible w/ C.I.
#----

ci_DS2P2 <- qnorm(c(.025, .975),mean=mu_l2P2, sd=sqrt(1/tau_l2P2))
cat("\n95% Credible Region for mean (", format(ci_DS2P2[1], digits = 4, nsmall = 3), 
    ",",format(ci_DS2P2[2], digits = 4, nsmall = 3),")")

#--- NEW
ci_DS2P2SD <- quantile(1/sqrt(TauG*(thetaP2+length(D2))), c(.025, .975))
cat("\n95% Credible Region for Std.Deviation (", format(ci_DS2P2SD[1], digits = 4, nsmall = 3), 
    ",",format(ci_DS2P2SD[2], digits = 4, nsmall = 3),")")
#---

# Prior 3 - Calculate the Posterior
#------------------------------------------
# Belief 1 - Average height is 66 inches
mu0 <- 66

# Belief 2 - Tau is Gamma(alpha,beta)
thetaP3 <- 0.1
alphaP3 <- 0.001
betaP3 <- 0.001

# Calculate the Posterior for each data-set
#------------------------------------------
# Data-set 1 - 
cat("-----------------------------------------")
cat("\nSummary Statistics for Data-Set1- (",D1,")\n")
alpha_l1P3 <- alphaP3+length(D1)/2
beta_l1P3 <- betaP3+0.5*sum((D1-mean(D1))^2)+(thetaP3*length(D1)*(mean(D1)-mu0)^2)/(2*(thetaP3+length(D1)))

## Analytically calculating the statistics

mu_l1P3 <- (thetaP3*mu0+length(D1)*mean(D1))/(thetaP3+length(D1))
cat("\nPosterior Mean is ",format(mu_l1P3, digits = 6, nsmall = 3))

tau_l1P3A <- 1/((beta_l1P3/((thetaP3+length(D1))*alpha_l1P3))*((2*alpha_l1P3)/(2*alpha_l1P3-2))) # v18
# tau_l1P3A <- 1/(alpha_l1P3/beta_l1P3^2)  # v19
cat("\n(Analytically) Posterior Std. Deviation is ",format(sqrt(1/tau_l1P3A), digits = 6, nsmall = 3)) # Analytically

TauG <- rgamma(10*N,alpha_l1P3, beta_l1P3); # mean(1/sqrt(TauG)); quantile(1/sqrt(TauG),c(.025, .5, .975))
m1 <- rnorm(10*N, mean=mu_l1P3, sd=1/sqrt(TauG*(thetaP3+length(D1))))
tau_l1P3S <- 1/var(m1)
cat("\n(Simulated) Posterior Std. Deviation is ",format(sqrt(1/tau_l1P3S), digits = 6, nsmall = 3)) # v18
# cat("\n(Simulated) Posterior Std. Deviation is ",format(sqrt(1/tau_l1P3S), digits = 6, nsmall = 3)) # Reviewed v19

tau_l1P3 <- tau_l1P3S  # Using the Simulated output to be compatible w/ C.I.

ci_DS1P3 <- qnorm(c(.025, .975),mean=mu_l1P3, sd=sqrt(1/tau_l1P3))
cat("\n95% Credible Region for mean (", format(ci_DS1P3[1], digits = 4, nsmall = 3), 
    ",",format(ci_DS1P3[2], digits = 4, nsmall = 3),")")

ci_DS1P3SD <- quantile(1/sqrt(TauG*(thetaP3+length(D1))), c(.025, .975))
cat("\n95% Credible Region for Std.Deviation (", format(ci_DS1P3SD[1], digits = 4, nsmall = 3), 
    ",",format(ci_DS1P3SD[2], digits = 4, nsmall = 3),")")

#------------------------------------------
# Data-set 2 - 
cat("\n\n-----------------------------------------")
cat("\nSummary Statistics for Data-Set2- (",D2,")\n")
alpha_l2P3 <- alphaP3+length(D2)/2
beta_l2P3 <- betaP3+0.5*sum((D2-mean(D2))^2)+(thetaP3*length(D2)*(mean(D2)-mu0)^2)/(2*(thetaP3+length(D2)))

## Analytically calculating the statistics

mu_l2P3 <- (thetaP3*mu0+length(D2)*mean(D2))/(thetaP3+length(D2))
cat("\nPosterior Mean is ",format(mu_l2P3, digits = 6, nsmall = 3))

#--- NEW
tau_l2P3A <- 1/((beta_l2P3/((thetaP3+length(D2))*alpha_l2P3))*((2*alpha_l2P3)/(2*alpha_l2P3-2)))
cat("\n(Analytically)Posterior Std. Deviation is ",format(sqrt(1/tau_l2P3A), digits = 6, nsmall = 3)) # Analytically

TauG <- rgamma(10*N,alpha_l2P3, beta_l2P3)
m2 <- rnorm(10*N, mean=mu_l2P3, sd=1/sqrt(TauG*(thetaP3+length(D2))))
tau_l2P3S <- 1/var(m2)
cat("\n(Simulated)Posterior Std. Deviation is ",format(sqrt(1/tau_l2P3S), digits = 6, nsmall = 3))

tau_l2P3 <- tau_l2P3S  # Using the Simulated output to be compatible w/ C.I.
#----

ci_DS2P3 <- qnorm(c(.025, .975),mean=mu_l2P3, sd=sqrt(1/tau_l2P3))
cat("\n95% Credible Region for mean (", format(ci_DS2P3[1], digits = 4, nsmall = 3), 
    ",",format(ci_DS2P3[2], digits = 4, nsmall = 3),")")

#--- NEW
ci_DS2P3SD <- quantile(1/sqrt(TauG*(thetaP3+length(D2))), c(.025, .975))
cat("\n95% Credible Region for Std.Deviation (", format(ci_DS2P3SD[1], digits = 4, nsmall = 3), 
    ",",format(ci_DS2P3SD[2], digits = 4, nsmall = 3),")")
#---

# Plot Sample densities to Compare Non-Informative vs. Informative Prior

xP1D1 <- rnorm(N,mean=mu_l1,sd=sqrt(1/tau_l1))
xP1D2 <- rnorm(N,mean=mu_l2,sd=sqrt(1/tau_l2))
xP2D1 <- rnorm(N,mean=mu_l1P2,sd=sqrt(1/((thetaP2+length(D1))*rgamma(N,alpha_l1P2,beta_l1P2))))
xP2D2 <- rnorm(N,mean=mu_l2P2,sd=sqrt(1/((thetaP2+length(D2))*rgamma(N,alpha_l2P2,beta_l2P2))))
xP3D1 <- rnorm(N,mean=mu_l1P3,sd=sqrt(1/((thetaP3+length(D1))*rgamma(N,alpha_l1P3,beta_l1P3))))
xP3D2 <- rnorm(N,mean=mu_l2P3,sd=sqrt(1/((thetaP3+length(D2))*rgamma(N,alpha_l2P3,beta_l2P3))))

D_Work <- as.tibble(rbind(data.frame(Height=xP1D1,Prior="1",Dataset="D1"),
                          data.frame(Height=xP1D2,Prior="1",Dataset="D2"),
                          data.frame(Height=xP2D1,Prior="2",Dataset="D1"),
                          data.frame(Height=xP2D2,Prior="2",Dataset="D2"),
                          data.frame(Height=xP3D1,Prior="3",Dataset="D1"),
                          data.frame(Height=xP3D2,Prior="3",Dataset="D2")))
D_Work %>%
  ggplot(aes(x = Height))+
  geom_density(aes(x=Height,group=Prior, colour=Prior), size=1.2)+ 
  labs(x = "Height", y = "Density") +
  scale_linetype_manual(name = "Dataset", values = c("solid", "dashed"))+
  facet_grid(~Dataset)+
  theme_bw()

# We will use the MArginal distribution of mu (t-Student as described) 

cat("-----------------------------------------")
cat("\nSummary Statistics of Predictive Distribution\n")

cat("\n#>>> Analytically ---\n")
cat("\nPosterior Mean is ",format(mu_l1P3, digits = 6, nsmall = 3))  # Revise V19

VarNewX <- (beta_l1P3/((thetaP3+length(D1))*alpha_l1P3))*((2*alpha_l1P3)/(2*alpha_l1P3-2))  # From Supplemental paper

cat("\nPosterior Std. Deviation is ",format(sqrt(VarNewX), digits = 6, nsmall = 3))

ci_NewX <- qt(c(.025, .975), 2*alpha_l1P3)*sqrt(beta_l1P3/((thetaP3+length(D1))*alpha_l1P3))+mu_l1P3  # Analytcially (v19)
cat("\n95% Credible Region for New Observation (", format(ci_NewX[1], digits = 4, nsmall = 3), 
    ",",format(ci_NewX[2], digits = 4, nsmall = 3),")")

# Simulating from t Distribution
cat("\n\n#>>> Simulated ---\n")
NewX <- rt(N,2*alpha_l1P3)*sqrt(beta_l1P3/((thetaP3+length(D1))*alpha_l1P3))+mu_l1P3  # Revised v19

cat("\nPosterior Mean is ",format(mean(NewX), digits = 6, nsmall = 3))  # V18
cat("\nPosterior Std. Deviation is ",format(sd(NewX), digits = 6, nsmall = 3))
# cat("\n\n",quantile(NewX, c(.025, .975)),"\n\n")
ci_NewX <- quantile(NewX, c(.025, .975))
cat("\n95% Credible Region for New Observation (", format(ci_NewX[1], digits = 4, nsmall = 3), 
    ",",format(ci_NewX[2], digits = 4, nsmall = 3),")")

ggplot(data=data.frame(Height=NewX),aes(x = Height))+
  geom_density(aes(x=Height), size=1.2, color="deepskyblue")+ 
  labs(x = "Height", y = "Density") +
  theme_bw()

# Build Dataframes to display C.I.s

TbCID1 <- data.frame(Prior = c("Prior #1","Prior #2","Prior #3"),
                     Min = c(ci_DS1[1],ci_DS1P2[1],ci_DS1P3[1]),
                     Max= c(ci_DS1[2],ci_DS1P2[2],ci_DS1P3[2]),
                     MM = c(mu_l1,mu_l1P2,mu_l1P3),
                     DS = c("D1","D1","D1"))

TbCID2 <- data.frame(Prior = c("Prior #1","Prior #2","Prior #3"),
                     Min = c(ci_DS2[1],ci_DS2P2[1],ci_DS2P3[1]),
                     Max= c(ci_DS2[2],ci_DS2P2[2],ci_DS2P3[2]),
                     MM = c(mu_l2,mu_l2P2,mu_l2P3),
                     DS = c("D2","D2","D2"))

TbCID1SD <- data.frame(Prior = c("Prior #1","Prior #2","Prior #3"),
                       Min = c(1/sqrt(tau_l1), ci_DS1P2SD[1],ci_DS1P3SD[1]),
                       Max= c(1/sqrt(tau_l1), ci_DS1P2SD[2],ci_DS1P3SD[2]),
                       MM = c(1/sqrt(tau_l1), 1/sqrt(tau_l1P2),1/sqrt(tau_l1P3)),
                       DS = c("D1", "D1","D1"))

TbCID2SD <- data.frame(Prior = c("Prior #1","Prior #2","Prior #3"),
                       Min = c(1/sqrt(tau_l2), ci_DS2P2SD[1],ci_DS2P3SD[1]),
                       Max= c(1/sqrt(tau_l2), ci_DS2P2SD[2],ci_DS2P3SD[2]),
                       MM = c(1/sqrt(tau_l2), 1/sqrt(tau_l2P2),1/sqrt(tau_l2P3)),
                       DS = c("D2", "D2","D2"))
# Plot C.I.s
ggplot(data = rbind(TbCID1,TbCID2)) +
  geom_point(aes(x=DS, y=MM, color = "Mean")) +
  geom_errorbar(aes(x=DS, ymin=Min, ymax=Max,
                    color = "C.I."), width = 0.5) +
  xlab("Dataset") +
  ylab("Mean") +
  scale_color_manual(name = "",
                     values = c("Mean" = "Blue",
                                "C.I." = "Black"))+
  facet_grid(~Prior)+
  theme_bw()

# Plot C.I.s
ggplot(data = rbind(TbCID1SD,TbCID2SD)) +
  geom_point(aes(x=DS, y=MM, color = "Std.Dev.")) +
  geom_errorbar(aes(x=DS, ymin=Min, ymax=Max,
                    color = "C.I."), width = 0.5) +
  xlab("Dataset") +
  ylab("Std. Deviation") +
  scale_color_manual(name = "",
                     values = c("Std.Dev." = "Blue",
                                "C.I." = "Black"))+
  facet_grid(~Prior)+
  theme_bw()


# Calculus of Theta0p and Tau0p

# Sampled Data
Ybar1 <- -1.82
S1 <- .21
n1 <- 209

theta0p <- (Ybar1*S1^(-2))/(S1^(-2))
tau0p <- (S1^(-2))

cat("\nParameters of Posterior are: Theta0p=",format(theta0p, digits = 6, nsmall = 3), 
    " and Tau0p=",format(tau0p, digits = 6, nsmall = 3)) 

# Calculus of Theta1p and Tau1p

# New Sampled Data
Ybar2 <- -1.02
S2 <- .28
n2 <- 79

theta1p <- (Ybar1*S1^(-2)+Ybar2*S2^(-2))/(S1^(-2)+S2^(-2))
tau1p <- (S1^(-2)+S2^(-2))

cat("\nParameters of Posterior are: Theta1p=",format(theta1p, digits = 6, nsmall = 3), 
    " and Tau1p=",format(tau1p, digits = 6, nsmall = 3)) 

# Calculus of Theta0p and Tau0p

# Sampled Data in reversed order
theta0p <- (Ybar2*S2^(-2))/(S2^(-2))
tau0p <- (S2^(-2))

cat("\nParameters of Posterior are: Theta0p=",format(theta0p, digits = 6, nsmall = 3), 
    " and Tau0p=",format(tau0p, digits = 6, nsmall = 3)) 

# Calculus of Theta1p and Tau1p
# New Sampled Data from p.o.v. of my colleague

theta1p <- (Ybar2*S2^(-2)+Ybar1*S1^(-2))/(S2^(-2)+S1^(-2))
tau1p <- (S2^(-2)+S1^(-2))

cat("\nParameters of Posterior are: Theta1p=",format(theta1p, digits = 6, nsmall = 3), 
    " and Tau1p=",format(tau1p, digits = 6, nsmall = 3)) 

# General Formula to calculate posterior belief of Theta_K

# Setup variables
DataHbA1c <- data.frame (
  Ybar = c(Ybar1, Ybar2, -1.9, -2.0, -1.21),
  S = c(S1, S2, 0.945, 0.285, 0.545)
)
K = nrow(DataHbA1c)

# Calculate Overal Theta_Kp

theta_p <- rep(0,K)
tau_p <- rep(0,K)

pCrit <- qnorm(.975,0,1)

# Calculating the Posterior parameters for each cumulative entry
theta_p <- cumsum((DataHbA1c$Ybar[1:K]*DataHbA1c$S[1:K]^(-2)))/cumsum(DataHbA1c$S[1:K]^(-2))
tau_p <- cumsum(DataHbA1c$S[1:K]^(-2))

# Calculating the 95% C.I. for each belief of theta
ci_p = data.frame (
  Dt = c("Y1","Y1toY2","Y1toY3","Y1toY4","Y1toY5"),
  Min = theta_p-sqrt(1/tau_p)*pCrit, 
  MM = theta_p,
  Max = theta_p+sqrt(1/tau_p)*pCrit
)

cat("\nParameters of Posterior are: Theta5p=",format(theta_p[K], digits = 6, nsmall = 3), 
    " and Tau5p=",format(tau_p[K], digits = 6, nsmall = 3)) 

# Plot C.I.s
ggplot(data = ci_p) +
  geom_point(aes(x=Dt, y=MM, color = "Theta")) +
  geom_errorbar(aes(x=Dt, ymin=Min, ymax=Max,
                    color = "C.I."), width = 0.5) +
  xlab("Data Incorporated to Model") +
  ylab(expression("theta")) +
  scale_color_manual(name = "",
                     values = c("Theta" = "Blue",
                                "C.I." = "Black"))+
  theme_bw()

# Simulating The change in the density of posterior belief as we add data
Y1 <- rnorm(N,mean=theta_p[1],sd=sqrt(1/tau_p[1]))
Y1toY2 <- rnorm(N,mean=theta_p[2],sd=sqrt(1/tau_p[2]))
Y1toY3 <- rnorm(N,mean=theta_p[3],sd=sqrt(1/tau_p[3]))
Y1toY4 <- rnorm(N,mean=theta_p[4],sd=sqrt(1/tau_p[4]))
Y1toY5 <- rnorm(N,mean=theta_p[5],sd=sqrt(1/tau_p[5]))

# Plot Sample densities to Compare Non-Informative vs. Informative Prior
D_WorkHbA1c <- as.tibble(rbind(data.frame(Theta=Y1, Dt="Y1"),
                          data.frame(Theta=Y1toY2,Dt="Y1toY2"),
                          data.frame(Theta=Y1toY3,Dt="Y1toY3"),
                          data.frame(Theta=Y1toY4,Dt="Y1toY4"),
                          data.frame(Theta=Y1toY5,Dt="Y1toY5")))
# The palette with black:
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp3 <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", 
          "#C3D7A4", "#52854C", "#4E84C4", "#293352")
D_WorkHbA1c %>%
  ggplot(aes(x = Theta))+
  geom_density(aes(x=Theta,group=Dt, colour=Dt), size=1.2)+ 
  labs(x = "Theta", y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()



