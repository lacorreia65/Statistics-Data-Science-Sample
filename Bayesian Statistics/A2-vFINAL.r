knitr::opts_chunk$set(echo = FALSE, message = FALSE, fig.width=5, fig.height=4.0)
library(tidyverse)
library(R2OpenBUGS)
library(kableExtra)
library(ggplot2)
library(forecast)

# The palette with black - Used in Graphs with :
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp3 <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", 
          "#C3D7A4", "#52854C", "#4E84C4", "#293352")

# Read database
bb_data <- read.csv("Data/brnbdy.csv")

# Introducing Scale Power transformation - Weisberg(2013) pp.189-190
scalePower <- function (Z,Lambda=0, Mod = FALSE) {
  # Introduced by Box, Cox (1964) for strictly positive Z
  if (Mod)
    gm = exp(sum(log(Z)/length(Z))) 
  else
    gm = 1.0
  if (Lambda == 0)
    return(gm*log(Z))
  else
    return(gm^(1-Lambda)*(Z^Lambda-1)/Lambda)
}

# Setup Transformed Variables
Y  <-  scalePower(bb_data$brain)
X  <-  scalePower(bb_data$body)


# Prepare dataframe and plot brain vs. body weight
df <- data.frame (lBrW = Y, lBdW = X)

# Plot Graph for 
df %>% 
  ggplot(aes(x=X, y=Y))+
  geom_point(size=1.8)+
  labs(x = "log(Body Weight)", y = "log(Brain Weight)") +
  ggtitle(NULL)+
  theme_bw()

# Remove Temporary Variable 
remove(df)

# Setup Initial Variables
Mu_alpha <- Mu_beta <- 0
Tau_alpha <- Tau_beta <- 0.0001
a_Tau<- b_Tau <- 0.0001  # v.16

# Setup Iteration Variables
N <- 20000

alpha <- beta <- tau <- rep(N,0)

# Setup \theta^{(0)}
alpha[1] <- beta[1] <- 0
tau[1] <- 1

# Calculates the value of parameters, given everything else
sample_Alpha <- function (X, Y, beta, tau, mu_alpha, tau_alpha) {
  tau_alphaStar <- length(Y)*tau+tau_alpha
  mu_alphaStar <- (tau*sum(Y-beta*X)+tau_alpha*mu_alpha)/tau_alphaStar
  alpha <- rnorm(1,mean=mu_alphaStar, sd=sqrt(1/tau_alphaStar))
  return(alpha)
}

sample_Beta <- function (X, Y, alpha, tau, mu_beta, tau_beta) {
  tau_betaStar<- tau*sum(X^2)+tau_beta
  mu_betaStar <- (tau*sum(X*(Y-alpha))+tau_beta*mu_beta)/tau_betaStar
  beta <- rnorm(1,mean=mu_betaStar, sd=sqrt(1/tau_betaStar))
  return(beta)
}

sample_Tau <- function (X, Y, alpha, beta, a_r, b_r) {
  a_rStar <- a_r+length(Y)/2
  b_rStar <- b_r+0.5*sum((Y-(alpha+beta*X))^2)
  tau <- rgamma(1,shape=a_rStar, rate=b_rStar)
  return(tau)
}

Gibbs <- function(X, Y, alpha, beta, tau, N) {
  for (i in 2:N) {
    alpha[i] <- sample_Alpha(X, Y, beta[i-1], tau[i-1], Mu_alpha, Tau_alpha)
    beta[i] <- sample_Beta(X, Y, alpha[i], tau[i-1], Mu_beta, Tau_beta)
    tau[i] <- sample_Tau(X, Y, alpha[i], beta[i], a_Tau, b_Tau)   # v.17
  }
  return(list(Alpha = alpha, Beta = beta, Tau = tau))
}

# Set Seed for reproducibility
set.seed(737)

#Calculates 3 chains
ThetaC1 <- Gibbs(X, Y, alpha, beta, tau, N)
ThetaC2 <- Gibbs(X, Y, alpha, beta, tau, N)
ThetaC3 <- Gibbs(X, Y, alpha, beta, tau, N)


# Plot ACF * ALPHA, BETA and TAU *

# Determine Burnin based on the convergence of ACF plot

# Build working dataframe
D_WorkAcfABT <- data.frame(ValCh1=ThetaC1$Alpha,
                           ValCh2=ThetaC1$Beta,
                           ValCh3=ThetaC1$Tau)


P1 <- ggAcf(D_WorkAcfABT$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(alpha)) +
  ggtitle(NULL)+
  theme_bw()
P2 <- ggAcf(D_WorkAcfABT$ValCh2, lag.max = 100)+
  labs(x = "Lag", y = expression(beta)) +
  ggtitle(NULL)+
  theme_bw()
P3 <- ggAcf(D_WorkAcfABT$ValCh3, lag.max = 100)+
  labs(x = "Lag", y = expression(tau)) +
  ggtitle(NULL)+
  theme_bw()

ggpubr::ggarrange(P1, P2, P3, ncol = 2, nrow = 2)

# Remove Working Variable to free memory
remove(D_WorkAcfABT)

# Set Burn-in parameter
Burnin <- 100

# Prepare Data-Frame for next steps
D1_Work <- rbind(data.frame(alpha=ThetaC1$Alpha[Burnin:N],
                            beta=ThetaC1$Beta[Burnin:N],
                            tau=ThetaC1$Tau[Burnin:N],Chain=factor(1)),
                 data.frame(alpha=ThetaC2$Alpha[Burnin:N],
                            beta=ThetaC2$Beta[Burnin:N],
                            tau=ThetaC2$Tau[Burnin:N],Chain=factor(2)),
                 data.frame(alpha=ThetaC3$Alpha[Burnin:N],
                            beta=ThetaC3$Beta[Burnin:N],
                            tau=ThetaC3$Tau[Burnin:N],Chain=factor(3)))


# Plot TracePlots * ALPHA, BETA and TAU *

set.seed(340)

Sz <- 5000

L <- N-Burnin+1

# Sampling 5000 points to generate "thinner" traceplots
S <- if (Sz <= L) sort(sample(Burnin:N, Sz, replace = FALSE)) else Burnin:N

# Build working dataframe
D_WorkTpAlpha <- data.frame(ValCh1=ThetaC1$Alpha[S],
                            ValCh2=ThetaC2$Alpha[S],
                            ValCh3=ThetaC3$Alpha[S])
P1 <- D_WorkTpAlpha %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(alpha)) +
  theme_bw()

# Build working dataframe
D_WorkTpBeta <- data.frame(ValCh1=ThetaC1$Beta[S],
                           ValCh2=ThetaC2$Beta[S],
                           ValCh3=ThetaC3$Beta[S])
P2 <- D_WorkTpBeta %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta)) +
  theme_bw()

# Build working dataframe
D_WorkTpTau <- data.frame(ValCh1=ThetaC1$Tau[S],
                          ValCh2=ThetaC2$Tau[S],
                          ValCh3=ThetaC3$Tau[S])
P3 <- D_WorkTpTau %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(tau)) +
  theme_bw()

# Build working dataframe
D_WorkTpAlphaBeta <- data.frame(ValAlpha=c(ThetaC1$Alpha[S],
                                           ThetaC2$Alpha[S],
                                           ThetaC3$Alpha[S]),
                                ValBeta=c(ThetaC1$Beta[S],
                                          ThetaC2$Beta[S],
                                          ThetaC3$Beta[S]))

P4 <- D_WorkTpAlphaBeta %>%
  ggplot(aes(x=ValAlpha, y=ValBeta))+
  geom_point(size=0.6)+ 
  labs(x = expression(alpha), y = expression(beta)) +
  theme_bw()

ggpubr::ggarrange(P1+theme(axis.title.x=element_blank(),
                           legend.position="none"), 
                  P2+theme(axis.title.x=element_blank(),
                           legend.position="none"), 
                  P3+theme(axis.title.x=element_blank(),
                           legend.position="none"), 
                  P4+theme(legend.position="none"), ncol = 2, nrow = 2)

# Remove Working Variable to free memory
remove(D_WorkTpAlpha, D_WorkTpBeta, D_WorkTpTau, D_WorkTpAlphaBeta)

# Plot the Graphs of Alpha
P1 <- D1_Work %>%
  ggplot(mapping = aes(x = alpha, group = Chain))+
  geom_density(aes(colour=Chain), size=0.5)+
  labs(x = expression(alpha), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

# Plot the Graphs of Beta
P2 <- D1_Work %>%
  ggplot(mapping = aes(x = beta, group = Chain))+
  geom_density(aes(colour=Chain), size=0.5)+
  labs(x = expression(beta), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

# Plot the Graphs of Tau
P3 <-D1_Work %>%
  ggplot(mapping = aes(x = tau, group = Chain))+
  geom_density(aes(colour=Chain), size=0.5)+
  labs(x = expression(tau), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

ggpubr::ggarrange(P1, P2, P3,ncol = 1, nrow = 3)


# Summary Statistics - alpha
prt_SummaryStat <- function (sample_, parameter) {
  m_ <- mean(sample_)
  sd_ <- sd(sample_)
  cat("\n---- Summary Statistics for ",parameter," -----")
  cat("\nMean",parameter,": ",format(m_, digits = 2, nsmall = 3)) 
  cat("\nStd.Dev.",parameter,": ",format(sd_, digits = 2, nsmall = 4))
  CI <- qnorm(p = c(0.025,0.975), mean = m_, sd = sd_)
  # cat("\nC.I. for",parameter,": [",CI[1],",",CI[2],"]\n")
  cat("\n95% Credible Region for",parameter,"(", format(CI[1], digits = 2, nsmall = 3), 
      ",",format(CI[2], digits = 4, nsmall = 3),")\n")
}
prt_SummaryStat(D1_Work %>% pull(alpha), "alpha")
prt_SummaryStat(D1_Work %>% pull(beta), "beta")
prt_SummaryStat(D1_Work %>% pull(tau), "tau")


# Using the Alphas and Betas found on item (a) We will generate a sample of exp(Y_55)
D2_Work <- D1_Work %>% 
  mutate(Y_55 = alpha+beta*log(55))

D2_Work %>%
  ggplot(mapping = aes(x = exp(Y_55), group = Chain))+
  geom_density(aes(colour=Chain), size=0.8)+
  labs(x = expression("Estimated Brain Weight | Body Weight=55Kg"), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

# Summary Statistics for the Distribution of an animal with average bodyweight of 55
m_a <- mean(exp(D2_Work$Y_55))
sd_a <- sd(exp(D2_Work$Y_55))
cat("\n---- Summary Statistics for Brain Weight | Body Weight = 55kg -----")
cat("\nMean : ",format(m_a, digits = 2, nsmall = 3)) 
cat("\nStd.Dev. : ",format(sd_a, digits = 2, nsmall = 4))
CI_a <- qnorm(p = c(0.025,0.975), mean = m_a, sd = sd_a)
cat("\n95% Credible Region (", format(CI_a[1], digits = 2, nsmall = 3), 
    ",",format(CI_a[2], digits = 2, nsmall = 3),")\n")
cat("\n>>> From this model, our belief of Brain Weight for an animal with \n")
cat("    average Body Weight of 55Kg is",format(m_a, digits = 2, nsmall = 2),
    "grams, with a 95% credible\n")
cat("    region of (",format(CI_a[1], digits = 2, nsmall = 2),",",
    format(CI_a[2], digits = 2, nsmall = 2),") grams.\n")

# Setup Data-set - Read Data
shist_data <- read.table(file="data/SmokeAgeDeath.csv",header=TRUE,sep=",")
colnames(shist_data) <- c("Smoking", "Age", "Deaths", "PersonYears")

# Setup Variables
NSmok=nrow(shist_data)  # Number of Items in data-set
NCSmok=4   # Number of classes of Smoking
NCAge=5    # Number of Classes of Age

# Setup OpenBugs running parameters
NSim <- 30000   # No. of simulations for productio
NChain <- 5     # No. of chains for production
NThin <- 5      # n.thin parameter for production
Burnin <- 10000 # Burn-In parameter for production
Sz <- 5000      # Size of samples for trace/acf plots


# Printing the Data-Frame
shist_data %>%
  kbl(booktabs = TRUE,
      caption = "Deaths by Lung Cancer, per Age, Smoking Level and Person-Years") %>% 
  kable_styling(latex_options = "striped")


cat("
model{
for(i in 1:NSmok)
{
#Age   Smoking   Deaths  PersonYears

Deaths[i]~dpois(lambda[i])
log(lambda[i]) <-  log(PersonYears[i]) + beta0 + beta.Age[Age[i]]+ beta.Smoking[Smoking[i]] + bias[i] 
bias[i] ~dnorm(0,tau)
bias.adj[i] <- bias[i] - mean(bias[]) 
}

# Loop for Age
for(ia in 1:NCAge){         
beta.Age[ia]~dnorm(0,tau.Age)
beta.Age.adj[ia] <- beta.Age[ia] - mean(beta.Age[]) 
}

# Loop for Smoking Level
for(is in 1:NCSmok){         
beta.Smoking[is]~dnorm(0,tau.Smoking)
beta.Smoking.adj[is] <- beta.Smoking[is] - mean(beta.Smoking[])
}

# Initial Priors - Verify
beta0 ~ dnorm(0, 0.008)   # Overall rate of disease
beta0.adj <- beta0 + mean(bias[]) + mean(beta.Age[])+ mean(beta.Smoking[])

# Initial Priors - Verify
std ~ dunif(0, 4)
tau <- 1/pow(std,2)

# Initial Priors - Verify
std.Age ~dunif(0, 4)
tau.Age <- 1/pow(std.Age,2)
std.Smoking ~ dunif(0,4) 
tau.Smoking <- 1/pow(std.Smoking,2)

}", file="LungCancerCentered.txt")

# Setup Parameters
params=c("beta.Age.adj",  "std.Age", "beta.Smoking.adj", "std.Smoking", 
         "beta0.adj", "std", "tau", "tau.Age", "tau.Smoking")

# Setup Initial Values
init.fun=function(){list(
  beta.Age=rnorm(5), std.Age=runif(1,0,4),
  beta.Smoking=rnorm(4), std.Smoking=runif(1,0,4), 
  std=runif(1,0,4), beta0=rnorm(1),
  bias=rnorm(20,0,.1))}


# Run Open Bugs
set.seed(2602)
model.data <- list("Age", "Smoking","Deaths","PersonYears","NSmok", 
                   "NCSmok", "NCAge")
attach(shist_data)
LungCanc_v0=bugs(model.data, init.fun, params, model.file="LungCancerCentered.txt",
                 working.directory = ".", n.chains=NChain, n.iter=NSim, 
                 n.burnin=Burnin, n.thin=NThin)
detach(shist_data)


# Get Simulation from OpenBugs
SArray= LungCanc_v0$sims.array   # Data Arrays
vname=attr(SArray,"dimnames")[3][[1]]  # Variable Names


# Print Summary statistics of parameters oif Interest
RNames <- rownames(LungCanc_v0[["summary"]][,c(1:3,5,7)]) # List of parameters
df_Prt <- data.frame(Parameter = RNames)
df_Prt <- cbind(df_Prt, as_tibble(LungCanc_v0[["summary"]][,c(1:3,5,7)]))
df_Prt %>%
  kbl(booktabs = TRUE, 
      caption = "OpenBugs Summary - Deaths from Lung Cancer") %>% 
  kable_styling(latex_options = "striped")

# Plot Sample densities to Compare Distributions over *BETA-0*
D_WorkBeta0 <- as_tibble(data.frame(LvlBeta=as.vector(SArray[,,"beta0.adj"])))

D_WorkBeta0 %>%
  ggplot(aes(x = LvlStd))+
  geom_density(aes(x=LvlBeta), size=1.2)+ 
  labs(x = expression(list(beta[0])), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

# Remove Working Variable to free memory
remove(D_WorkBeta0)


# Plot Sample densities to Compare Distributions over *AGE*
D_WorkAge <- as_tibble(rbind(data.frame(LvlAge=as.vector(SArray[,,paste0("beta.Age.adj[", 1, "]")]),Age="0-45"),
                             data.frame(LvlAge=as.vector(SArray[,,paste0("beta.Age.adj[", 2, "]")]),Age="45-54"),
                             data.frame(LvlAge=as.vector(SArray[,,paste0("beta.Age.adj[", 3, "]")]),Age="55-64"),
                             data.frame(LvlAge=as.vector(SArray[,,paste0("beta.Age.adj[", 4, "]")]),Age="65-74"),
                             data.frame(LvlAge=as.vector(SArray[,,paste0("beta.Age.adj[", 5, "]")]),Age="74+")))


P0 <- D_WorkAge %>%
  ggplot(aes(x = LvlAge))+
  geom_density(aes(x=LvlAge,group=Age, colour=Age), size=1.2)+ 
  labs(x = expression(list(beta[Age_][i])), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

P0 + geom_segment(aes(x = -2, y = 2, xend = 1.5, yend = 4), size = 0.3,
                arrow = arrow(length = unit(0.5, "cm")))

# Remove Working Variable to free memory
remove(D_WorkAge)

# Plot Sample densities to Compare Distributions over *SMOKE LEVEL*
D_WorkSMoke <- as_tibble(rbind(data.frame(LvlSmoke=as.vector(SArray[,,paste0("beta.Smoking.adj[", 1, "]")]),Level="Lv1 - Never Smoked"),
                               data.frame(LvlSmoke=as.vector(SArray[,,paste0("beta.Smoking.adj[", 2, "]")]),Level="Lv2 - Past Smoker"),
                               data.frame(LvlSmoke=as.vector(SArray[,,paste0("beta.Smoking.adj[", 3, "]")]),Level="Lv3 - 1-20 Cig/Day"),
                               data.frame(LvlSmoke=as.vector(SArray[,,paste0("beta.Smoking.adj[", 4, "]")]),Level="Lv4 - +20 Cig/Day")))


D_WorkSMoke %>%
  ggplot(aes(x = LvlSmoke))+
  geom_density(aes(x=LvlSmoke,group=Level, colour=Level), size=1.2)+ 
  labs(x = expression(list(beta[Smoke_][i])), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

# Remove Working Variable to free memory
remove(D_WorkSMoke)

# Plot Sample densities to Compare Distributions over *STD DEVIATION* and *PRECISION*
D_WorkStd <- as_tibble(rbind(data.frame(LvlStd=as.vector(SArray[,,"std.Age"]),Level="Age"),
                             data.frame(LvlStd=as.vector(SArray[,,"std.Smoking"]),Level="Smoking")))

D_WorkPrc <- as_tibble(rbind(data.frame(LvlStd=as.vector(SArray[,,"tau.Age"]),Level="Age"),
                             data.frame(LvlStd=as.vector(SArray[,,"tau.Smoking"]),Level="Smoking")))


P1 <- D_WorkStd %>%
  ggplot(aes(x = LvlStd))+
  geom_density(aes(x=LvlStd,group=Level, colour=Level), size=1.2)+ 
  labs(x = expression(list(sigma[Age]~","~sigma[Smoke])), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

P2 <- D_WorkPrc %>%
  ggplot(aes(x = LvlStd))+
  geom_density(aes(x=LvlStd,group=Level, colour=Level), size=1.2)+ 
  labs(x = expression(list(tau[Age]~","~tau[Smoke])), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

ggpubr::ggarrange(P1, P2, ncol = 1, nrow = 2)

# Remove Working Variable to free memory
remove(D_WorkStd, D_WorkPrc)


# Plot TracePlots * BETA.AGE *

# Sampling 5000 points to generate "thinner" traceplots
set.seed(312)

L <- NSim-Burnin

S <- if (Sz <= L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe

P <- list()  

for (i in 1:5) {
  # Build working dataframe
  D_WorkTpAge <- data.frame(ValCh1=SArray[S,1,paste0("beta.Age.adj[", i, "]")],
                            ValCh2=SArray[S,2,paste0("beta.Age.adj[", i, "]")],
                            ValCh3=SArray[S,3,paste0("beta.Age.adj[", i, "]")],
                            ValCh4=SArray[S,4,paste0("beta.Age.adj[", i, "]")],
                            ValCh5=SArray[S,5,paste0("beta.Age.adj[", i, "]")])
  
  p <- D_WorkTpAge %>%
    ggplot(aes(seq(from=1,to=Sz)))+
    geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
    geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
    geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
    geom_line(aes(y=ValCh4, colour=4), size=0.8)+ 
    geom_line(aes(y=ValCh5, colour=5), size=0.8)+ 
    labs(y = eval(bquote(expression(beta[Age_][.(i)])))) +
    theme_bw()
  
  P <- c(P, list(p))
}

ggpubr::ggarrange(P[[1]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"), 
                  P[[2]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[3]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[4]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[5]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"), ncol = 3, nrow = 2)

# Remove Working Variable to free memory
remove(p, P, D_WorkTpAge)

# Plot TracePlots * BETA.SMOKING *

# Sampling 5000 points to generate "thinner" traceplots
set.seed(351)

L <- NSim-Burnin

S <- if (Sz <= L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe

P <- list()  

for (i in 1:4) {
  # Build working dataframe
  D_WorkTpSmok <- data.frame(ValCh1=SArray[S,1,paste0("beta.Smoking.adj[", i, "]")],
                             ValCh2=SArray[S,2,paste0("beta.Smoking.adj[", i, "]")],
                             ValCh3=SArray[S,3,paste0("beta.Smoking.adj[", i, "]")],
                             ValCh4=SArray[S,4,paste0("beta.Smoking.adj[", i, "]")],
                             ValCh5=SArray[S,5,paste0("beta.Smoking.adj[", i, "]")])
  
  p <- D_WorkTpSmok %>%
    ggplot(aes(seq(from=1,to=Sz)))+
    geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
    geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
    geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
    geom_line(aes(y=ValCh4, colour=4), size=0.8)+ 
    geom_line(aes(y=ValCh5, colour=5), size=0.8)+ 
    labs(y = eval(bquote(expression(beta[Smoke_][.(i)])))) +
    theme_bw()
  P <- c(P, list(p))
}

ggpubr::ggarrange(P[[1]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"), 
                  P[[2]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[3]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[4]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"), ncol = 2, nrow = 2)

# Remove Working Variable to free memory
remove(p, P, D_WorkTpSmok)

# Plot TracePlots * BETA.0 *

# Sampling 5000 points to generate "thinner" traceplots
set.seed(761)

L <- NSim-Burnin

S <- if (Sz <= L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe
D_WorkTpBeta0 <- data.frame(ValCh1=SArray[S,1,paste0("beta0.adj")],
                            ValCh2=SArray[S,2,paste0("beta0.adj")],
                            ValCh3=SArray[S,3,paste0("beta0.adj")],
                            ValCh4=SArray[S,4,paste0("beta0.adj")],
                            ValCh5=SArray[S,5,paste0("beta0.adj")])

P0 <- D_WorkTpBeta0 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  geom_line(aes(y=ValCh4, colour=4), size=0.8)+ 
  geom_line(aes(y=ValCh5, colour=5), size=0.8)+ 
  labs(x = "Iteration", y = expression(beta[0])) +
  theme_bw()

P0 + theme(axis.title.x=element_blank(),
           legend.position="none")

# Remove Working Variable to free memory
remove(D_WorkTpBeta0)

# Plot TracePlots * STD.AGE and STD.SMOKING *

# Sampling 5000 points to generate "thinner" traceplots
set.seed(390)

L <- NSim-Burnin

S <- if (Sz <= L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe
D_WorkTpStdAge <- data.frame(ValCh1=SArray[S,1,paste0("std.Age")],
                             ValCh2=SArray[S,2,paste0("std.Age")],
                             ValCh3=SArray[S,3,paste0("std.Age")],
                             ValCh4=SArray[S,4,paste0("std.Age")],
                             ValCh5=SArray[S,5,paste0("std.Age")])

P1 <- D_WorkTpStdAge %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  geom_line(aes(y=ValCh4, colour=4), size=0.8)+ 
  geom_line(aes(y=ValCh5, colour=5), size=0.8)+ 
  labs(x = "Iteration", y = expression(sigma[Age])) +
  theme_bw()

# Build working dataframe
D_WorkTpStdSmok <- data.frame(ValCh1=SArray[S,1,paste0("std.Smoking")],
                              ValCh2=SArray[S,2,paste0("std.Smoking")],
                              ValCh3=SArray[S,3,paste0("std.Smoking")],
                              ValCh4=SArray[S,4,paste0("std.Smoking")],
                              ValCh5=SArray[S,5,paste0("std.Smoking")])

P2 <- D_WorkTpStdSmok %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  geom_line(aes(y=ValCh4, colour=4), size=0.8)+ 
  geom_line(aes(y=ValCh5, colour=5), size=0.8)+ 
  labs(x = "Iteration", y = expression(sigma[Smoke])) +
  theme_bw()

ggpubr::ggarrange(P1+theme(axis.title.x=element_blank(),
                           legend.position="none"), 
                  P2+theme(axis.title.x=element_blank(),
                           legend.position="none"), ncol = 1, nrow = 2)

# Remove Working Variable to free memory
remove(D_WorkTpStdAge, D_WorkTpStdSmok)

# Plot TracePlots * BETA.AGE *

P <- list()  

for (i in 1:5) {
  # Build working dataframe
  D_WorkAcfAge <- data.frame(ValCh1=SArray[,1,paste0("beta.Age.adj[", i, "]")])
  
  p <- ggAcf(D_WorkAcfAge$ValCh1, lag.max = 100)+
    labs(x = "Lag", y = eval(bquote(expression(beta[Age_][.(i)])))) +
    ggtitle(NULL)+
    theme_bw()
  P <- c(P, list(p))
}

ggpubr::ggarrange(P[[1]]+theme(axis.title.x=element_blank(),
                               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
                               legend.position="none"), 
                  P[[2]]+theme(axis.title.x=element_blank(),
                               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
                               legend.position="none"),
                  P[[3]]+theme(axis.title.x=element_blank(),
                               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
                               legend.position="none"),
                  P[[4]]+theme(axis.title.x=element_blank(),
                               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
                               legend.position="none"),
                  P[[5]]+theme(axis.title.x=element_blank(),
                               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
                               legend.position="none"), ncol = 3, nrow = 2)

# Remove Working Variable to free memory
remove(p, P, D_WorkAcfAge)

# Plot ACF * BETA.SMOKING *

# Group 1-4 - Elaborate a List of Graphs

P <- list()  

for (i in 1:4) {
  # Build working dataframe
  D_WorkAcfSmok <- data.frame(ValCh1=SArray[,1,paste0("beta.Smoking.adj[", i, "]")])
  
  p <- ggAcf(D_WorkAcfSmok$ValCh1, lag.max = 100)+
    labs(x = "Lag", y = eval(bquote(expression(beta[Smoke_][.(i)])))) +
    ggtitle(NULL)+
    theme_bw()
  P <- c(P, list(p))
}

ggpubr::ggarrange(P[[1]]+theme(axis.title.x=element_blank(),
                               legend.position="none"), 
                  P[[2]]+theme(axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[3]]+theme(axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[4]]+theme(axis.title.x=element_blank(),
                               legend.position="none"), ncol = 2, nrow = 2)

# Remove Working Variable to free memory
remove(p, P, D_WorkAcfSmok)


# Plot TracePlots * BETA.0 *

# Build working dataframe
D_WorkAcfBeta0 <- data.frame(ValCh1=SArray[,1,paste0("beta0.adj")])

P1 <- ggAcf(D_WorkAcfBeta0$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(beta[0])) +
  ggtitle(NULL)+
  theme_bw()

P1+theme(axis.title.x=element_blank(),
         legend.position="none")

# Remove Working Variable to free memory
remove(D_WorkAcfBeta0)

# Comparison of Smoking Categories
SMokCat <- rbind(LungCanc_v0[["summary"]]["beta.Smoking.adj[1]",c(1:3,5,7)],
                 LungCanc_v0[["summary"]]["beta.Smoking.adj[4]",c(1:3,5,7)])
RNames <- c("beta.Smoking.adj[1]", "beta.Smoking.adj[4]")
df_Prt <- data.frame(Parameter = RNames)
df_Prt <- cbind(df_Prt, as_tibble(SMokCat))
df_Prt %>%
  kbl(booktabs = TRUE, 
      caption = "Smoking Categories - Never Smoked vs. >20 cigarettes per day") %>% 
  kable_styling(latex_options = "striped")

```

```{r}
# Calculate Probability Increase between categories
ProbInc <- exp(SArray[,,"beta.Smoking.adj[4]"]-SArray[,,"beta.Smoking.adj[1]"])
m_Smoke <- mean(ProbInc)
sd_Smoke <- sd(ProbInc)
cat("\n---- Summary Statistics for Increase in Probability -----")
cat("\nMean : ",format(m_Smoke, digits = 2, nsmall = 3)) 
cat("\nStd.Dev. : ",format(sd_Smoke, digits = 2, nsmall = 4))
CI_Smoking <- qnorm(c(0.025,0.975), mean=m_Smoke, sd=sd_Smoke)
cat("\n95% Credible Region (", format(CI_Smoking[1], digits = 2, nsmall = 3), 
    ",",format(CI_Smoking[2], digits = 2, nsmall = 3),")\n")
cat("\n>>> From this model, our belief of increase probability of death by lung cancer\n")
cat("    of someone who smokes >20 cigarettes per day vs. who never smoked is",format(m_Smoke, digits = 2, nsmall = 2))
cat("\n    with corresponding 95% credible region of (",    format(CI_Smoking[1], digits = 2, nsmall = 2),",",
    format(CI_Smoking[2], digits = 2, nsmall = 2),")\n")

# Setup - Read Data Diabetics
diab_data <- read.table(file="data/DiabetesDrugEffect.csv",header=TRUE,sep=",")
NDiab <- nrow(diab_data)

# Setup OpenBugs running parameters

NSim <- 30000   # No. of simulations for productio
NChain <- 5     # No. of chains for production
NThin <- 5      # n.thin parameter for production
Burnin <- 10000 # Burn-In parameter for production
Sz <- 5000      # Size of samples for trace/acf plots


### **************** M O D E L  -  1 ****************** ###

### ---- SETUP MODELS ---- ###

# MODEL 1 - Splitted
cat("
model{
for(i in 1:NDiab)
{
#StudyID  StudyN  diff  Sediff

diff[i]~dnorm(delta[i], tau[i]) #v.19
delta[i]~dnorm(theta, tau0)
tau[i] <- 1/pow(Sediff[i],2)

}

# Initial Priors - Verify
theta ~dnorm(-1.5, 30)  # v.19
sigma0 ~dgamma(2, 1)    # v.19
tau0 <- 1/pow(sigma0,2)

}", file="DiabetDrug_M1.txt")

### ---- PARAM SETTINGS ---- ###

# Setup Parameters
paramsM1=c("theta", "tau0", "sigma0", "delta")

### ---- INITIALIZATION PROCEDURES ---- ###

# Setup Initial Values for Model 1
init.funM1=function(){list(
  theta = rnorm(1, mean=-1.5, sd=sqrt(.3)),
  sigma0=rgamma(1, 2, 1))}  # v.19

### ---- PROCESS MODELS ON OPENBUGS ---- ###

# Run Open Bugs - Model 1
set.seed(2602)
model1.data <- list("StudyID", "StudyN", "diff", "Sediff","NDiab")
attach(diab_data)
DiabetDrug_M1=bugs(model1.data, init.funM1, paramsM1, model.file="DiabetDrug_M1.txt",
                   working.directory = ".", n.chains=NChain, n.iter=NSim, n.burnin=Burnin, n.thin=NThin)
detach(diab_data)

### ---- ORGANIZE OUTPUTS ---- ###

# Get Simulation from OpenBugs - Model 1
DArrayM1= DiabetDrug_M1$sims.array   # Data Arrays
vname=attr(DArrayM1,"dimnames")[3][[1]]  # Variable Names


### **************** M O D E L  -  2 ****************** ###

### ---- SETUP MODELS ---- ###

cat("
model{
for(i in 1:NDiab)
{
#StudyID  StudyN  diff  Sediff

diff[i]~dnorm(theta, tau[i])
tau[i] <- 1/(pow(sigma0,2)+pow(Sediff[i],2))
}

# Initial Priors - Verify
theta ~dnorm(-1.5, 30)  # v.19
sigma0 ~dgamma(2, 1)    # v.19

}", file="DiabetDrug_M2.txt")

### ---- PARAM SETTINGS ---- ###

# Setup Parameters
paramsM2=c("theta", "sigma0")


### ---- INITIALIZATION PROCEDURES ---- ###

# Setup Initial Values for Model 2
init.funM2=function(){list(
  theta = rnorm(1, mean=-1.5, sd=sqrt(.3)),
  sigma0=rgamma(1, 2, 1))}  # v.19

### ---- PROCESS MODELS ON OPENBUGS ---- ###

# Run Open Bugs - Model 2
set.seed(8421)
model2.data <- list("StudyID", "StudyN", "diff", "Sediff","NDiab")
attach(diab_data)
DiabetDrug_M2=bugs(model2.data, init.funM2, paramsM2, model.file="DiabetDrug_M2.txt",
                   working.directory = ".", n.chains=NChain, n.iter=NSim, n.burnin=Burnin, n.thin=NThin)
detach(diab_data)

### ---- ORGANIZE OUTPUTS ---- ###

# Get Simulation from OpenBugs - Model 2
DArrayM2= DiabetDrug_M2$sims.array   # Data Arrays
vname=attr(DArrayM2,"dimnames")[3][[1]]  # Variable Names


# Print Summary statistics of parameters oif Interest
RNames <- rownames(DiabetDrug_M1[["summary"]][,c(1:3,5,7)]) # List of parameters
df_Prt <- data.frame(Parameter = RNames)
df_Prt <- cbind(df_Prt, as_tibble(DiabetDrug_M1[["summary"]][,c(1:3,5,7)]))
df_Prt %>%
  kbl(booktabs = TRUE, 
      caption = "Model 1 - Difference between Drug A and B") %>% 
  kable_styling(latex_options = "striped")


# Print Summary statistics of parameters oif Interest
RNames <- rownames(DiabetDrug_M2[["summary"]][,c(1:3,5,7)]) # List of parameters
df_Prt <- data.frame(Parameter = RNames)
df_Prt <- cbind(df_Prt, as_tibble(DiabetDrug_M2[["summary"]][,c(1:3,5,7)]))
df_Prt %>%
  kbl(booktabs = TRUE, 
      caption = "Model 2 - Difference between Drug A and B") %>% 
  kable_styling(latex_options = "striped")

# Plot TracePlots and ACF of *THETA *

# Sampling 5000 points to generate "thinner" traceplots
set.seed(3961)

L <- NSim - Burnin

S <- if (Sz <= L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe
D_WorkTpTheta <- data.frame(ValCh1=DArrayM1[S,1,paste0("theta")],
                            ValCh2=DArrayM1[S,2,paste0("theta")],
                            ValCh3=DArrayM1[S,3,paste0("theta")],
                            ValCh4=DArrayM1[S,4,paste0("theta")],
                            ValCh5=DArrayM1[S,5,paste0("theta")])

P1 <- D_WorkTpTheta %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  geom_line(aes(y=ValCh4, colour=4), size=0.8)+ 
  geom_line(aes(y=ValCh5, colour=5), size=0.8)+ 
  labs(y = expression(theta)) +
  theme_bw()

P2 <- ggAcf(DArrayM1[,1,paste0("theta")], lag.max = 100)+
  labs(x = "Lag", y = expression(theta)) +
  ggtitle(NULL)+
  theme_bw()

ggpubr::ggarrange(P1+theme(axis.title.x=element_blank(),
                           legend.position="none"), P2, ncol = 1, nrow = 2)

# Remove Working Variable to free memory
remove(D_WorkTpTheta)

# Plot TracePlots * SIGMA0*

# Sampling 5000 points to generate "thinner" traceplots
set.seed(2391)

L <- NSim - Burnin

S <- if (Sz <= L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe
D_WorkTpSigma0 <- data.frame(ValCh1=DArrayM1[S,1,paste0("sigma0")],
                             ValCh2=DArrayM1[S,2,paste0("sigma0")],
                             ValCh3=DArrayM1[S,3,paste0("sigma0")],
                             ValCh4=DArrayM1[S,4,paste0("sigma0")],
                             ValCh5=DArrayM1[S,5,paste0("sigma0")])

P1 <- D_WorkTpSigma0 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  geom_line(aes(y=ValCh4, colour=4), size=0.8)+ 
  geom_line(aes(y=ValCh5, colour=5), size=0.8)+ 
  labs(y = expression(sigma[0])) +
  theme_bw()

P2 <- ggAcf(DArrayM1[,1,paste0("sigma0")], lag.max = 100)+
  labs(x = "Lag", y = expression(sigma[0])) +
  ggtitle(NULL)+
  theme_bw()

ggpubr::ggarrange(P1+theme(axis.title.x=element_blank(),
                           legend.position="none"), P2, ncol = 1, nrow = 2)

# Remove Working Variable to free memory
remove(D_WorkTpSigma0)

# Plot TracePlots for Deltas

# Sampling 5000 points to generate "thinner" traceplots
set.seed(746)

L <- NSim - Burnin

S <- if (Sz <= L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Group 1-4 - Elaborate a List of Graphs

P <- list()  

for (i in 1:4) {
  # Build working dataframe
  D_WorkTpDelta <- data.frame(ValCh1=DArrayM1[S,1,paste0("delta[", i, "]")],
                              ValCh2=DArrayM1[S,2,paste0("delta[", i, "]")],
                              ValCh3=DArrayM1[S,3,paste0("delta[", i, "]")],
                              ValCh4=DArrayM1[S,4,paste0("delta[", i, "]")],
                              ValCh5=DArrayM1[S,5,paste0("delta[", i, "]")])
  p <- D_WorkTpDelta %>%
    ggplot(aes(seq(from=1,to=Sz)))+
    geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
    geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
    geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
    geom_line(aes(y=ValCh4, colour=4), size=0.8)+ 
    geom_line(aes(y=ValCh5, colour=5), size=0.8)+ 
    labs(y = eval(bquote(expression(delta[.(i)])))) +
    theme_bw()
  P <- c(P, list(p))
}

ggpubr::ggarrange(P[[1]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"), 
                  P[[2]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[3]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[4]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"), ncol = 2, nrow = 2)

# Remove Working Variable to free memory
remove(p, P, D_WorkTpDelta)

# Plot TracePlots for Deltas

# Sampling 5000 points to generate "thinner" traceplots
set.seed(456)

L <- NSim - Burnin

S <- if (Sz <= L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Group 5-8 - Elaborate a List of Graphs

P <- list()  

for (i in 5:8) {
  # Build working dataframe
  D_WorkTpDelta <- data.frame(ValCh1=DArrayM1[S,1,paste0("delta[", i, "]")],
                              ValCh2=DArrayM1[S,2,paste0("delta[", i, "]")],
                              ValCh3=DArrayM1[S,3,paste0("delta[", i, "]")],
                              ValCh4=DArrayM1[S,4,paste0("delta[", i, "]")],
                              ValCh5=DArrayM1[S,5,paste0("delta[", i, "]")])
  p <- D_WorkTpDelta %>%
    ggplot(aes(seq(from=1,to=Sz)))+
    geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
    geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
    geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
    geom_line(aes(y=ValCh4, colour=4), size=0.8)+ 
    geom_line(aes(y=ValCh5, colour=5), size=0.8)+ 
    labs(y = eval(bquote(expression(delta[.(i)])))) +
    theme_bw()
  P <- c(P, list(p))
}

ggpubr::ggarrange(P[[1]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"), 
                  P[[2]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[3]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[4]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"), ncol = 2, nrow = 2)

# Remove Working Variable to free memory
remove(p, P, D_WorkTpDelta)

# Plot TracePlots for Deltas

# Sampling 5000 points to generate "thinner" traceplots
set.seed(198)

L <- NSim - Burnin

S <- if (Sz <= L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Group 9-12 - Elaborate a List of Graphs

P <- list()  

for (i in 9:12) {
  # Build working dataframe
  D_WorkTpDelta <- data.frame(ValCh1=DArrayM1[S,1,paste0("delta[", i, "]")],
                              ValCh2=DArrayM1[S,2,paste0("delta[", i, "]")],
                              ValCh3=DArrayM1[S,3,paste0("delta[", i, "]")],
                              ValCh4=DArrayM1[S,4,paste0("delta[", i, "]")],
                              ValCh5=DArrayM1[S,5,paste0("delta[", i, "]")])
  p <- D_WorkTpDelta %>%
    ggplot(aes(seq(from=1,to=Sz)))+
    geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
    geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
    geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
    geom_line(aes(y=ValCh4, colour=4), size=0.8)+ 
    geom_line(aes(y=ValCh5, colour=5), size=0.8)+ 
    labs(y = eval(bquote(expression(delta[.(i)])))) +
    theme_bw()
  P <- c(P, list(p))
}

ggpubr::ggarrange(P[[1]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"), 
                  P[[2]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[3]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[4]]+theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank(),
                               legend.position="none"), ncol = 2, nrow = 2)

# Remove Working Variable to free memory
remove(p, P, D_WorkTpDelta)

# Plot ACF plots for Deltas

# Sampling 5000 points to generate "thinner" traceplots

# Group 1-4 - Elaborate a List of Graphs

P <- list()  

for (i in 1:4) {
  # Build working dataframe
  D_WorkTpDelta <- data.frame(ValCh1=DArrayM1[,1,paste0("delta[", i, "]")])
  
  p <- ggAcf(D_WorkTpDelta$ValCh1, lag.max = 100)+
    labs(x = "Lag", y = eval(bquote(expression(delta[.(i)])))) +
    ggtitle(NULL)+
    theme_bw()
  P <- c(P, list(p))
}

ggpubr::ggarrange(P[[1]]+theme(axis.title.x=element_blank(),
                               legend.position="none"), 
                  P[[2]]+theme(axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[3]]+theme(axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[4]]+theme(axis.title.x=element_blank(),
                               legend.position="none"), ncol = 2, nrow = 2)

# Remove Working Variable to free memory
remove(p, P, D_WorkTpDelta)

# Plot ACF plots for Deltas

# Sampling 5000 points to generate "thinner" traceplots

# Group 5-8 - Elaborate a List of Graphs

P <- list()  

for (i in 5:8) {
  # Build working dataframe
  D_WorkTpDelta <- data.frame(ValCh1=DArrayM1[,1,paste0("delta[", i, "]")])
  
  p <- ggAcf(D_WorkTpDelta$ValCh1, lag.max = 100)+
    labs(x = "Lag", y = eval(bquote(expression(delta[.(i)])))) +
    ggtitle(NULL)+
    theme_bw()
  P <- c(P, list(p))
}

ggpubr::ggarrange(P[[1]]+theme(axis.title.x=element_blank(),
                               legend.position="none"), 
                  P[[2]]+theme(axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[3]]+theme(axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[4]]+theme(axis.title.x=element_blank(),
                               legend.position="none"), ncol = 2, nrow = 2)

# Remove Working Variable to free memory
remove(p, P, D_WorkTpDelta)

# Plot ACF plots for Deltas

# Group 9-12 - Elaborate a List of Graphs

P <- list()  

for (i in 9:12) {
  # Build working dataframe
  D_WorkTpDelta <- data.frame(ValCh1=DArrayM1[,1,paste0("delta[", i, "]")])
  
  p <- ggAcf(D_WorkTpDelta$ValCh1, lag.max = 100)+
    labs(x = "Lag", y = eval(bquote(expression(delta[.(i)])))) +
    ggtitle(NULL)+
    theme_bw()
  P <- c(P, list(p))
}

ggpubr::ggarrange(P[[1]]+theme(axis.title.x=element_blank(),
                               legend.position="none"), 
                  P[[2]]+theme(axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[3]]+theme(axis.title.x=element_blank(),
                               legend.position="none"),
                  P[[4]]+theme(axis.title.x=element_blank(),
                               legend.position="none"), ncol = 2, nrow = 2)
# Remove Working Variable to free memory
remove(p, P, D_WorkTpDelta)

# Plot TracePlots and ACF of *THETA *

# Sampling 5000 points to generate "thinner" traceplots
set.seed(3871)

L <- NSim - Burnin

S <- if (Sz <= L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe
D_WorkTpTheta <- data.frame(ValCh1=DArrayM2[S,1,paste0("theta")],
                            ValCh2=DArrayM2[S,2,paste0("theta")],
                            ValCh3=DArrayM2[S,3,paste0("theta")],
                            ValCh4=DArrayM2[S,4,paste0("theta")],
                            ValCh5=DArrayM2[S,5,paste0("theta")])

P1 <- D_WorkTpTheta %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  geom_line(aes(y=ValCh4, colour=4), size=0.8)+ 
  geom_line(aes(y=ValCh5, colour=5), size=0.8)+ 
  labs(y = expression(theta)) +
  theme_bw()

P2 <- ggAcf(DArrayM2[,1,paste0("theta")], lag.max = 100)+
  labs(x = "Lag", y = expression(theta)) +
  ggtitle(NULL)+
  theme_bw()

ggpubr::ggarrange(P1+theme(axis.title.x=element_blank(),
                           legend.position="none"), P2, ncol = 1, nrow = 2)

# Remove Working Variable to free memory
remove(D_WorkTpTheta)

# Plot TracePlots * SIGMA0*

# Sampling 5000 points to generate "thinner" traceplots
set.seed(8463)

L <- NSim - Burnin

S <- if (Sz <= L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe
D_WorkTpSigma0 <- data.frame(ValCh1=DArrayM2[S,1,paste0("sigma0")],
                             ValCh2=DArrayM2[S,2,paste0("sigma0")],
                             ValCh3=DArrayM2[S,3,paste0("sigma0")],
                             ValCh4=DArrayM2[S,4,paste0("sigma0")],
                             ValCh5=DArrayM2[S,5,paste0("sigma0")])

P1 <- D_WorkTpSigma0 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  geom_line(aes(y=ValCh4, colour=4), size=0.8)+ 
  geom_line(aes(y=ValCh5, colour=5), size=0.8)+ 
  labs(y = expression(sigma[0])) +
  theme_bw()

P2 <- ggAcf(DArrayM2[,1,paste0("sigma0")], lag.max = 100)+
  labs(x = "Lag", y = expression(sigma[0])) +
  ggtitle(NULL)+
  theme_bw()

ggpubr::ggarrange(P1+theme(axis.title.x=element_blank(),
                           legend.position="none"), P2, ncol = 1, nrow = 2)

# Remove Working Variable to free memory
remove(D_WorkTpSigma0)

# Plot Sample densities to Compare Distributions of *THETA* and *SIGMA0*

set.seed(9023)

L <- NSim-Burnin

S <- if (Sz <= L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

D_WorkThetas <- as_tibble(rbind(data.frame(LvlStd=as.vector(DArrayM1[S,,paste0("theta")]), Model = "M1"),
                                data.frame(LvlStd=as.vector(DArrayM2[S,,paste0("theta")]), Model = "M2")))
D_WorkSigmas <- as_tibble(rbind(data.frame(LvlStd=as.vector(DArrayM1[S,,paste0("sigma0")]), Model = "M1"),
                                data.frame(LvlStd=as.vector(DArrayM2[S,,paste0("sigma0")]), Model = "M2")))


P1 <- D_WorkThetas %>%
  ggplot(aes(x = LvlStd))+
  geom_density(aes(x=LvlStd, group=Model, colour=Model), size=0.8)+ 
  labs(x = expression(list(theta[Model_][1]~"vs."~theta[Model_][2])), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

P2 <- D_WorkSigmas %>%
  ggplot(aes(x = LvlStd))+
  geom_density(aes(x=LvlStd, group=Model, colour=Model), size=0.8)+ 
  labs(x = expression(list((sigma[0])[Model_][1]~"vs."~(sigma[0])[Model_][2])), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

ggpubr::ggarrange(P1, P2, ncol = 1, nrow = 2)

# Remove Working Variable to free memory
remove(D_WorkThetas, D_WorkSigmas)


# Prepare Data - Get delta[1-12]
D3_Work <- as.data.frame(DiabetDrug_M1[["summary"]][4:15,c(1:3,5,7)])

# Generates Studies labels
StudyLst <- rapply(list(diab_data$StudyID), sprintf, fmt = "%02d", how = "replace")

# Collects C.I's of observed values from studies
CI_Study <- data.frame(ID = StudyLst[[1]],   
                       lower = diab_data$diff-diab_data$Sediff, 
                       estim = diab_data$diff, 
                       upper = diab_data$diff+diab_data$Sediff, 
                       Source = "Study")

# Collects C.R's from Summary report 
CI_Estim <- data.frame(ID = StudyLst[[1]], 
                       lower = D3_Work$`2.5%`, 
                       estim = D3_Work$mean, 
                       upper = D3_Work$`97.5%`, 
                       Source = "Model-1")

# Organize data to generate graphs & analysis
D4_Work <- rbind(CI_Study, CI_Estim)

#  (i) Identifying the direction of change of estimates
# (ii) Identifying differences o length of CR and CI for delta and "diff"
P1 <- D4_Work %>% 
  ggplot() +
  geom_point(aes(x=ID, y=estim, color = Source, shape = Source, size = 0.5)) +
  geom_errorbar(aes(x=ID, ymin=lower, ymax=upper,
                    color = Source), width = 1) +
  xlab("Study No.") +
  ylab(expression("C.R.+"~delta[i]~" and 'diff'")) +
  scale_color_manual(values = cbp2)+
  facet_wrap(~Source)+
  theme_bw()

P1 + theme(legend.position="none",
           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))

# (i) Identifying the direction of posterior mean of delta and "diff"
P2 <- D4_Work %>% 
  ggplot() +
  geom_point(aes(x=ID, y=estim, colour = Source, shape = Source), size=5) +
  xlab(expression("Study No.")) +
  ylab(expression(delta[i]~" and 'diff'")) +
  scale_color_manual(values = cbp2)+
  theme_bw()

P2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))

remove(D3_Work, D4_Work, StudyLst)

