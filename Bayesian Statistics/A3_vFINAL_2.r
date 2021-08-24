library(tidyverse)
library(R2OpenBUGS)
library(kableExtra)
library(ggplot2)
library(forecast)
library(coda)
# library(boa)

# The palette with black - Used in Graphs with :
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp3 <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", 
          "#C3D7A4", "#52854C", "#4E84C4", "#293352")

#
#data from Healy, page 90.
#  (MJR Healy, 1988, Glim: An Introduction, Clarendon Press: Oxford.)
#  Looking to see if Smoking is a risk factor for hypertension, controlling for obesity, snoring, and gender
#  Note 1: there was no males or females who were smokers and obese and who did not snore (so 1 1 0 had no exposures)
#  Note 2: here we are simply looking at the effect of smoking given the other factors.  We are ignoring the possibility that 
#  obesity might be related to smoking or that snoring might be strongly effected by smoking and obesity.  
#  In modern epi, these factors might be consider to be in the <<causal path>> and perhaps you might not control for them in this way.
cat(
  "smoke  obese  snore male hypoten n
0 0 0 1 5 60
0 0 0 0 10 149
1 0 0 1 2 17
1 0 0 0 6 16
0 1 0 1 1 12
0 1 0 0 2 9
0 0 1 1 36 187
0 0 1 0 28 138
1 0 1 1 13 85
1 0 1 0 4 39
0 1 1 1 15 51
0 1 1 0 11 28
1 1 1 1 8 23
1 1 1 0 4 12
", file= "Data/SmokeHyperData.txt")

# Setup Data-set - Read Data
SmokeHyper=read.table("Data/SmokeHyperData.txt",header=TRUE,sep = "")

# Setup Variables
N=nrow(SmokeHyper)  # Number of Items in data-set

# Setup OpenBugs running parameters
NSim <- 30000    # No. of simulations for productio
NChain <- 3      # No. of chains for production
NThin <- 50      # n.thin parameter for production
Burnin <- 10000  # Burn-In parameter for production
Sz <- 5000       # Size of samples for trace/acf plots


# Printing the Data-Frame
SmokeHyper %>%
  kbl(booktabs = TRUE, digits = 4,
      caption = "Data - Relationship between Hypertension and Smoking") %>% 
  kable_styling(latex_options = "striped")


cat("
model{
  for(i in 1:N){
   hypoten[i] ~ dbin(mu[i], n[i])
   logit(mu[i]) <- b0 + b.smok*smoke[i]+ b.ob*obese[i]+ b.sn*snore[i] + 
     b.male*male[i] + b.smsn*smoke[i]*snore[i] + b[i]
    b[i] ~dnorm(0, tau.b)
   }
  b.smok ~ dnorm(0, .04) # so, sd =5.  exp(5) ~ 148 which is huge
  b.ob ~ dnorm(0, .04) 
  b.sn ~ dnorm(0, .04) 
  b.male ~ dnorm(0, .04) 
  b0 ~ dnorm(0, .04) 
  b.smsn ~dnorm(0, .04)
  sd.b ~ dunif(0, 5)
  tau.b <- 1/pow(sd.b,2)
  }
  ", file="SmokeHyperMod3.txt")

# Setup Parameters
paramsM0=c("b0", "b.smok", "b.ob", "b.sn", "b.male", "b.smsn" , "sd.b")

bugM0.dat=list("hypoten", "n", "smoke", "obese", "snore", "male", "N")  # what variable you need in the model


# Setup Initial Values
initM0.fun=function(){ list(  b=runif(N,-.8,-.2), 
                              b0=runif(1,-.8,-.2),
                              b.smok=runif(1,-.8,-.2),b.ob=runif(1,-.8,-.2), b.sn=runif(1,-.8,-.2),
                              b.male=runif(1,-.8,-.2), b.smsn=runif(1, -8,-.2), sd.b=runif(1,.2,.8)	
) }

#### Could change the code below...
# Run Open Bugs - NO BURNIN / NO THINN
set.seed(2602)
attach(SmokeHyper)
SmokeHypeBaseM0=bugs(bugM0.dat, initM0.fun, paramsM0, model.file="SmokeHyperMod3.txt",
                     n.chains=NChain, n.iter=NSim, n.burnin=0, n.thin=1, debug=FALSE)
detach(SmokeHyper)


# Get Simulation from OpenBugs
SArrayM0= SmokeHypeBaseM0$sims.array   # Data Arrays
vname=attr(SArrayM0,"dimnames")[3][[1]]  # Variable Names


# Print Summary statistics of parameters oif Interest
RNames <- rownames(SmokeHypeBaseM0[["summary"]][,c(1:3,5,7)]) # List of parameters
df_Prt <- data.frame(Parameter = RNames)
df_Prt <- cbind(df_Prt, as_tibble(SmokeHypeBaseM0[["summary"]][,c(1:3,5,7)]))
df_Prt %>%
  kbl(booktabs = TRUE, digits = 4, 
      caption = "OpenBugs Summary - Hypertension and Smoke Study (Run-0)") %>% 
  kable_styling(latex_options = "striped")

# Plot TracePlots * SIGMA, BETA0, BETA.SMOKE and BETA.OBESITY *

# Sampling 5000 points to generate "thinner" traceplots
set.seed(312)

L <- NSim-Burnin

OldSz <- Sz  
Sz <- L       # Override Original value

S <- if (Sz < L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe
D_WorkSigma <- data.frame(ValCh1=SArrayM0[S,1,paste0("sd.b")],
                          ValCh2=SArrayM0[S,2,paste0("sd.b")],
                          ValCh3=SArrayM0[S,3,paste0("sd.b")])

D_WorkBeta0 <- data.frame(ValCh1=SArrayM0[S,1,paste0("b0")],
                          ValCh2=SArrayM0[S,2,paste0("b0")],
                          ValCh3=SArrayM0[S,3,paste0("b0")])

D_WorkBetaSmok <- data.frame(ValCh1=SArrayM0[S,1,paste0("b.smok")],
                             ValCh2=SArrayM0[S,2,paste0("b.smok")],
                             ValCh3=SArrayM0[S,3,paste0("b.smok")])

D_WorkBetaObes <- data.frame(ValCh1=SArrayM0[S,1,paste0("b.ob")],
                             ValCh2=SArrayM0[S,2,paste0("b.ob")],
                             ValCh3=SArrayM0[S,3,paste0("b.ob")])

# Plot Traceplots for selected data-frames
P1 <- D_WorkSigma %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(sigma)) +
  theme_bw()
P2 <- D_WorkBeta0 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta[0])) +
  theme_bw()
P3 <- D_WorkBetaSmok %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta[Smoke])) +
  theme_bw()
P4 <- D_WorkBetaObes %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta[Obesity])) +
  theme_bw()

# Plot all Graphs in a same frame
ggpubr::ggarrange(P1+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"), 
                  P2+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"),
                  P3+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"),
                  P4+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"), ncol = 2, nrow = 2)

# Remove Working Variable to free memory
Sz <- OldSz # Restore original value
remove(OldSz, D_WorkSigma, D_WorkBeta0, D_WorkBetaSmok, D_WorkBetaObes )

# Plot ACF Plots * SIGMA, BETA0, BETA.SMOKE and BETA.OBESITY *

# Build working dataframe
D_WorkSigma <- data.frame(ValCh1=SArrayM0[,1,paste0("sd.b")])
D_WorkBeta0 <- data.frame(ValCh1=SArrayM0[,1,paste0("b0")])
D_WorkBetaSmok <- data.frame(ValCh1=SArrayM0[,1,paste0("b.smok")])
D_WorkBetaObes <- data.frame(ValCh1=SArrayM0[,1,paste0("b.ob")])

# Plot ACF Plots for selected data-frames
P1 <- ggAcf(D_WorkSigma$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(sigma)) +
  ggtitle(NULL)+
  theme_bw()
P2 <- ggAcf(D_WorkBeta0$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(beta[0])) +
  ggtitle(NULL)+
  theme_bw()
P3 <- ggAcf(D_WorkBetaSmok$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(beta[Smoke])) +
  ggtitle(NULL)+
  theme_bw()
P4 <- ggAcf(D_WorkBetaObes$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(beta[Obesity])) +
  ggtitle(NULL)+
  theme_bw()

# Plot all Graphs in a same frame
ggpubr::ggarrange(P1+theme(axis.title.x=element_blank(),
                           legend.position="none"), 
                  P2+theme(axis.title.x=element_blank(),
                           legend.position="none"),
                  P3+theme(axis.title.x=element_blank(),
                           legend.position="none"),
                  P4+theme(axis.title.x=element_blank(),
                           legend.position="none"), ncol = 2, nrow = 2)
# Remove Working Variable to free memory
remove(D_WorkSigma, D_WorkBeta0, D_WorkBetaSmok, D_WorkBetaObes )

# Setup Parameters for Run-01

paramsM1 <- paramsM0

bugM1.dat <- bugM0.dat

initM1.fun <- initM0.fun

# Run Open Bugs - BURNIN=10,000 / THINN=50
set.seed(2212)
attach(SmokeHyper)
SmokeHypeBaseM1=bugs(bugM1.dat, initM1.fun, paramsM1, model.file="SmokeHyperMod3.txt",
                     n.chains=NChain, n.iter=NSim, n.burnin=Burnin, n.thin=NThin , debug=FALSE)
detach(SmokeHyper)


# Get Simulation from OpenBugs
SArrayM1= SmokeHypeBaseM1$sims.array   # Data Arrays
vname=attr(SArrayM1,"dimnames")[3][[1]]  # Variable Names


# Print Summary statistics of parameters oif Interest
RNames <- rownames(SmokeHypeBaseM1[["summary"]][,c(1:3,5,7)]) # List of parameters
df_Prt <- data.frame(Parameter = RNames)
df_Prt <- cbind(df_Prt, as_tibble(SmokeHypeBaseM1[["summary"]][,c(1:3,5,7)]))
df_Prt %>%
  kbl(booktabs = TRUE, digits = 4, 
      caption = "OpenBugs Summary - Hypertension and Smoke Study (Run-1)") %>% 
  kable_styling(latex_options = "striped")

# Plot TracePlots * SIGMA, BETA0, BETA.SMOKE and BETA.OBESITY *

# Sampling 5000 points to generate "thinner" traceplots
set.seed(312)

L <- NSim-Burnin

OldSz <- Sz  
Sz <- L       # Override Original value

S <- if (Sz < L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe
D_WorkSigma <- data.frame(ValCh1=SArrayM1[S,1,paste0("sd.b")],
                          ValCh2=SArrayM1[S,2,paste0("sd.b")],
                          ValCh3=SArrayM1[S,3,paste0("sd.b")])

D_WorkBeta0 <- data.frame(ValCh1=SArrayM1[S,1,paste0("b0")],
                          ValCh2=SArrayM1[S,2,paste0("b0")],
                          ValCh3=SArrayM1[S,3,paste0("b0")])

D_WorkBetaSmok <- data.frame(ValCh1=SArrayM1[S,1,paste0("b.smok")],
                             ValCh2=SArrayM1[S,2,paste0("b.smok")],
                             ValCh3=SArrayM1[S,3,paste0("b.smok")])

D_WorkBetaObes <- data.frame(ValCh1=SArrayM1[S,1,paste0("b.ob")],
                             ValCh2=SArrayM1[S,2,paste0("b.ob")],
                             ValCh3=SArrayM1[S,3,paste0("b.ob")])

# Plot Traceplots for selected data-frames
P1 <- D_WorkSigma %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(sigma)) +
  theme_bw()
P2 <- D_WorkBeta0 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta[0])) +
  theme_bw()
P3 <- D_WorkBetaSmok %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta[Smoke])) +
  theme_bw()
P4 <- D_WorkBetaObes %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta[Obesity])) +
  theme_bw()

# Plot all Graphs in a same frame
ggpubr::ggarrange(P1+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"), 
                  P2+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"),
                  P3+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"),
                  P4+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"), ncol = 2, nrow = 2)

# Remove Working Variable to free memory
Sz <- OldSz # Restore original value
remove(OldSz, D_WorkSigma, D_WorkBeta0, D_WorkBetaSmok, D_WorkBetaObes )

# Plot ACF Plots * SIGMA, BETA0, BETA.SMOKE and BETA.OBESITY *

# Build working dataframe
D_WorkSigma <- data.frame(ValCh1=SArrayM1[,1,paste0("sd.b")])
D_WorkBeta0 <- data.frame(ValCh1=SArrayM1[,1,paste0("b0")])
D_WorkBetaSmok <- data.frame(ValCh1=SArrayM1[,1,paste0("b.smok")])
D_WorkBetaObes <- data.frame(ValCh1=SArrayM1[,1,paste0("b.ob")])

# Plot ACF Plots for selected data-frames
P1 <- ggAcf(D_WorkSigma$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(sigma)) +
  ggtitle(NULL)+
  theme_bw()
P2 <- ggAcf(D_WorkBeta0$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(beta[0])) +
  ggtitle(NULL)+
  theme_bw()
P3 <- ggAcf(D_WorkBetaSmok$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(beta[Smoke])) +
  ggtitle(NULL)+
  theme_bw()
P4 <- ggAcf(D_WorkBetaObes$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(beta[Obesity])) +
  ggtitle(NULL)+
  theme_bw()

# Plot all Graphs in a same frame
ggpubr::ggarrange(P1+theme(axis.title.x=element_blank(),
                           legend.position="none"), 
                  P2+theme(axis.title.x=element_blank(),
                           legend.position="none"),
                  P3+theme(axis.title.x=element_blank(),
                           legend.position="none"),
                  P4+theme(axis.title.x=element_blank(),
                           legend.position="none"), ncol = 2, nrow = 2)
# Remove Working Variable to free memory
remove(D_WorkSigma, D_WorkBeta0, D_WorkBetaSmok, D_WorkBetaObes )

# Prepare data to plot densities to Compare Distributions of 
# * SIGMA, BETA0, BETA.SMOKE and BETA.OBESITY *

set.seed(9023)

L <- NSim-Burnin

S <- if (Sz < L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe
D_WorkSigma <- as_tibble(rbind(data.frame(LvlStd=as.vector(SArrayM1[S,1,paste0("sd.b")]), Chain = "C1"),
                               data.frame(LvlStd=as.vector(SArrayM1[S,2,paste0("sd.b")]), Chain = "C2"),
                               data.frame(LvlStd=as.vector(SArrayM1[S,3,paste0("sd.b")]), Chain = "C3")))
D_WorkBeta0 <- as_tibble(rbind(data.frame(LvlStd=as.vector(SArrayM1[S,1,paste0("b0")]), Chain = "C1"),
                               data.frame(LvlStd=as.vector(SArrayM1[S,2,paste0("b0")]), Chain = "C2"),
                               data.frame(LvlStd=as.vector(SArrayM1[S,3,paste0("b0")]), Chain = "C3")))
D_WorkBetaSmok <- as_tibble(rbind(data.frame(LvlStd=as.vector(SArrayM1[S,1,paste0("b.smok")]), Chain = "C1"),
                                  data.frame(LvlStd=as.vector(SArrayM1[S,2,paste0("b.smok")]), Chain = "C2"),
                                  data.frame(LvlStd=as.vector(SArrayM1[S,3,paste0("b.smok")]), Chain = "C3")))
D_WorkBetaObes <- as_tibble(rbind(data.frame(LvlStd=as.vector(SArrayM1[S,1,paste0("b.ob")]), Chain = "C1"),
                                  data.frame(LvlStd=as.vector(SArrayM1[S,2,paste0("b.ob")]), Chain = "C2"),
                                  data.frame(LvlStd=as.vector(SArrayM1[S,3,paste0("b.ob")]), Chain = "C3")))

# Calculates Batch Means and respective SE - Author: Michael Escobar (thank you!)
CalcBatchMeans <- function (x, Batn = 50) {
  BigN=length(x)
  BatInc=ceiling((1:BigN)/(BigN/Batn) )
  BM=tapply(x,BatInc,mean)
  return(list(MCE=(sd(BM)/sqrt(length(BM))), BM=BM))
}

# Calculates the Standard Error via Auto-Correlation Function - Author: Michael Escobar (thank you!)
CalcAcSe=function(x,lag.max=50){
  autoc=(acf(x,lag.max=lag.max,plot=FALSE))$acf
  sd(x)/sqrt(length(x))*sqrt(-1+2*sum(autoc))
}

CalcChainStats <- function (parm, chain, TxtParm) {
  L <- CalcBatchMeans(SArrayM1[,chain,paste0(parm)])
  S <- CalcAcSe(SArrayM1[,chain,paste0(parm)],120)
  cat("\n---- Summary Statistics for",TxtParm,"(Chain ",chain,") -----")
  cat("\nB-Mean : ",format(mean(L$BM), digits = 2, nsmall = 4)) 
  cat("\nB-S.E.: ",format(L$MCE, digits = 2, nsmall = 6))
  CI_up <- mean(L$BM)+2*L$MCE; CI_lo <- mean(L$BM)-2*L$MCE; 
  cat("\n95% Credible Region (", 
      format(CI_lo, digits = 2, nsmall = 4), ",",
      format(CI_up, digits = 2, nsmall = 4),")\n")
  cat("\n(ACF) S.E.: ",format(S, digits = 2, nsmall = 6))
  CI_up <- mean(L$BM)+2*S; CI_lo <- mean(L$BM)-2*S; 
  cat("\n(ACF) 95% Credible Region (", 
      format(CI_lo, digits = 2, nsmall = 4), ",",
      format(CI_up, digits = 2, nsmall = 4),")\n")
  cat("---------------------------------------------------\n")
}

# Calculates Posterior Means using Batch-Means and ACF - SIGMA
CalcChainStats("sd.b",1,"SIGMA")
CalcChainStats("sd.b",2,"SIGMA")
CalcChainStats("sd.b",3,"SIGMA")

D_WorkSigma %>%
  ggplot(aes(x = LvlStd))+
  geom_density(aes(x=LvlStd, group=Chain, colour=Chain), size=0.8)+ 
  labs(x = expression(sigma), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

# Calculates Posterior Means using Batch-Means and ACF - BETA0
CalcChainStats("b0",1,"BETA0")
CalcChainStats("b0",2,"BETA0")
CalcChainStats("b0",3,"BETA0")


D_WorkBeta0 %>%
  ggplot(aes(x = LvlStd))+
  geom_density(aes(x=LvlStd, group=Chain, colour=Chain), size=0.8)+ 
  labs(x = expression(beta[0]), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

# Calculates Posterior Means using Batch-Means and ACF - BETA-SMOKE
CalcChainStats("b.smok",1,"BETA-SMOKE")
CalcChainStats("b.smok",2,"BETA-SMOKE")
CalcChainStats("b.smok",3,"BETA-SMOKE")



D_WorkBetaSmok %>%
  ggplot(aes(x = LvlStd))+
  geom_density(aes(x=LvlStd, group=Chain, colour=Chain), size=0.8)+ 
  labs(x = expression(beta[Smoke]), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

# Calculates Posterior Means using Batch-Means and ACF - BETA-OBESITY
CalcChainStats("b.ob",1,"BETA-OBESITY")
CalcChainStats("b.ob",2,"BETA-OBESITY")
CalcChainStats("b.ob",3,"BETA-OBESITY")


D_WorkBetaObes %>%
  ggplot(aes(x = LvlStd))+
  geom_density(aes(x=LvlStd, group=Chain, colour=Chain), size=0.8)+ 
  labs(x = expression(beta[Obesity]), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

# Remove Working Variables and workdataframes used in this question to free memory
remove(D_WorkSigma, D_WorkBeta0, D_WorkBetaSmok, D_WorkBetaObes )

# Prepare Structures for CODA
Chn.sigma <- list( C1=mcmc(SArrayM1[,1,paste0("sd.b")]),
                   C2=mcmc(SArrayM1[,2,paste0("sd.b")]),
                   C3=mcmc(SArrayM1[,3,paste0("sd.b")]))
Chn.b0 <- list( C1=mcmc(SArrayM1[,1,paste0("b0")]),
                C2=mcmc(SArrayM1[,2,paste0("b0")]),
                C3=mcmc(SArrayM1[,3,paste0("b0")]))
Chn.bsmok <- list( C1=mcmc(SArrayM1[,1,paste0("b.smok")]),
                   C2=mcmc(SArrayM1[,2,paste0("b.smok")]),
                   C3=mcmc(SArrayM1[,3,paste0("b.smok")]))
Chn.bob <- list( C1=mcmc(SArrayM1[,1,paste0("b.ob")]),
                 C2=mcmc(SArrayM1[,2,paste0("b.ob")]),
                 C3=mcmc(SArrayM1[,3,paste0("b.ob")]))

# Calculate Summary MCMC of estimates
CalcDiags <- function (MC, TxtParm) {
  MMC <- cbind(MC$C1,MC$C2,MC$C3)
  colnames(MMC) <- c("Chain1","Chain2","Chain3")
  cat("\n--------- CODA Summary Statistics for",TxtParm,"---------\n")
  print(summary(mcmc(MMC)))
  cat(">>> Effective Size:\n")
  print(effectiveSize(MMC))
  cat("\n>>> GEWEKE Diagnostics:")
  print(geweke.diag(MMC))
  geweke.plot(mcmc(MMC))
  cat(">>> GELMAN Diagnostics:\n")
  print(gelman.diag(MC))
  gelman.plot(MC)
}


CalcDiags(Chn.sigma, "SIGMA")

CalcDiags(Chn.b0, "BETA-0")

CalcDiags(Chn.bsmok, "BETA-SMOKE")

CalcDiags(Chn.bob, "BETA-OBESITY")

# Setup Data-set - Read Data
HarvestingData <- data.frame( x=c(16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46),
                              y=c(2508,2518,3304,3423,3057,3190,3500,3883,3823,3646,3708,
                                  3333,3517,3241,3103,2776))

# Setup Variables
N=nrow(HarvestingData)  # Number of Items in data-set

# Setup OpenBugs running parameters
NSim <- 30000    # No. of simulations for productio
NChain <- 3      # No. of chains for production
NThin <- 10      # n.thin parameter for production
Burnin <- 10000  # Burn-In parameter for production
Sz <- 5000       # Size of samples for trace/acf plots


# Printing the Data-Frame
HarvestingData %>%
  kbl(booktabs = TRUE, digits = 1,
      caption = "(x):Harvesting Date (in days) vs. (y):Yield (in Kg/Ha)") %>% 
  kable_styling(latex_options = "striped")


HarvestingData %>%
  ggplot(aes(x=x, y=y))+
  geom_point(size=2.5)+ 
  labs(x = "X (# of days after flowering)", y = "Y (Yield in Kg/ha)") +
  theme_bw()

# MODEL Q2_1
cat("### Model 1 -
model{
  for (i in 1:N) {
    y[i]~dnorm(mu[i],tau)
    mu[i]<- b[1] + b[2]*(x[i]-31)
 
    # Get the residuals for the observed value
    res[i] <- (y[i]-mu[i])           # Estimate the residual for this model  - Item (1)
    stdres[i] <- res[i]*sqrt(tau)    # Calculate the standard residuals - Item (2)
    
    dev1.obs[i] <- pow(res[i],2)
    dev2.obs[i] <- pow(stdres[i],2)
    
    # Get a replicated Sample - sample of the predictive distribution
    y.rep[i]~dnorm(mu[i],tau)
    p.smaller[i] <- step(y[i]-y.rep[i])  # Check to see the probability of getting extreme value - Item (3)
    
    # Residual and Moments replicated - this gives the predicted distribution for these values
    res.rep[i] <- y.rep[i]-mu[i]
    stdres.rep[i] <- res.rep[i]*sqrt(tau)
    
    dev1.rep[i] <- pow(res.rep[i],2)
    dev2.rep[i] <- pow(stdres.rep[i],2)
    
    # Likelihood for each observed and replicated data
    # note: Need to know the density function of the probability model
    
    loglike[i] <- (0.5)*log(tau/6.283185) + (-0.5)*tau*pow((y[i]-mu[i]),2)
    loglike.rep[i]<-  (0.5)*log(tau/6.283185) + (-0.5)*tau*pow((y.rep[i]-mu[i]),2)
    
    p.inv[i]<- 1/exp(loglike[i])      #  This is to find the predictive ordinate of the observations
  }
  
  # Prior definitions
  b[1]~dnorm(0,.000001)
  b[2]~dnorm(0,.000001)
  tau~dgamma(.0001,.0001)

  # Summing Diagnostics Values
  chidev1.obs <- sum(dev1.obs[])
  chidev2.obs <- sum(dev2.obs[])

  chidev1.rep <- sum(dev1.rep[] )
  chidev2.rep <- sum(dev2.rep[] )

  chidev1.pval<-step(chidev1.obs-chidev1.rep)
  chidev2.pval<-step(chidev2.obs-chidev2.rep)

  #   Deviance statistic
  dev<-   -2*sum(loglike[])
  dev.rep <-  -2*sum(loglike.rep[])
  dev.pval<-step(dev-dev.rep)

}", file="HarvestModQ21.txt")


# Setup Parameters
bugMQ21.dat<-list("x", "y", "N")

initMQ21.fun<-function(){ list(b=runif(2,-.8,-.2), tau=runif(1,.2,.8), y.rep=rnorm(N))} 

paramsMQ21<-c("b", "tau",          # the rest are for the model checking
              "mu", "res", "stdres", "res.rep", "stdres.rep", "p.smaller",
              "p.inv", "chidev1.pval", "chidev2.pval", "chidev1.obs", "chidev2.obs",
              "chidev1.rep", "chidev2.rep", "dev", "dev.rep", "dev.pval")


#### Could change the code below...
# Run Open Bugs - 
set.seed(2157)
attach(HarvestingData)
HarvestingMQ21=bugs(bugMQ21.dat, initMQ21.fun, paramsMQ21, model.file="HarvestModQ21.txt",
                    n.chains=NChain, n.iter=NSim, n.burnin=Burnin, n.thin=NThin)
detach(HarvestingData)


# Get Simulation from OpenBugs
SArrayMQ21= HarvestingMQ21$sims.array   # Data Arrays
vname=attr(SArrayMQ21,"dimnames")[3][[1]]  # Variable Names
MQ21.coda <- as.mcmc.list(HarvestingMQ21)


# Print Summary statistics of parameters oif Interest
RNames <- rownames(HarvestingMQ21[["summary"]][,c(1:3,5,7)]) # List of parameters
df_Prt <- data.frame(Parameter = RNames)
df_Prt <- cbind(df_Prt, as_tibble(HarvestingMQ21[["summary"]][,c(1:3,5,7)]))
df_Prt %>%
  filter(substr(Parameter,1,2)=="b[" |
           substr(Parameter,1,3)=="tau") %>% 
  kbl(booktabs = TRUE, digits = 4,
      caption = "OpenBugs Summary - Harvesting Date vs. Yield (Model-1)") %>% 
  kable_styling(latex_options = "striped")

# Plot TracePlots * b[1], b[2] and tau *

# Sampling 5000 points to generate "thinner" traceplots
set.seed(312)

L <- NSim-Burnin

#OldSz <- Sz  
#Sz <- L       # Override Original value

S <- if (Sz < L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe
D_Workb1<- data.frame(ValCh1=SArrayMQ21[S,1,paste0("b[1]")],
                      ValCh2=SArrayMQ21[S,2,paste0("b[1]")],
                      ValCh3=SArrayMQ21[S,3,paste0("b[1]")])

D_Workb2 <- data.frame(ValCh1=SArrayMQ21[S,1,paste0("b[2]")],
                       ValCh2=SArrayMQ21[S,2,paste0("b[2]")],
                       ValCh3=SArrayMQ21[S,3,paste0("b[2]")])

D_Worktau <- data.frame(ValCh1=SArrayMQ21[S,1,paste0("tau")],
                        ValCh2=SArrayMQ21[S,2,paste0("tau")],
                        ValCh3=SArrayMQ21[S,3,paste0("tau")])

# Plot Traceplots for selected data-frames
P1 <- D_Workb1 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta[1])) +
  theme_bw()
P2 <- D_Workb2 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta[2])) +
  theme_bw()
P3 <- D_Worktau %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(tau)) +
  theme_bw()

# Plot all Graphs in a same frame
ggpubr::ggarrange(P1+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"), 
                  P2+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"),
                  P3+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"), ncol = 2, nrow = 2)

# Remove Working Variable to free memory
# Sz <- OldSz # Restore original value
remove(D_Workb1, D_Workb2, D_Worktau )

# Plot ACF Plots * b[1], b[2] and tau *

# Build working dataframe
D_Workb1 <- data.frame(ValCh1=SArrayMQ21[,1,paste0("b[1]")])
D_Workb2 <- data.frame(ValCh1=SArrayMQ21[,1,paste0("b[2]")])
D_Worktau <- data.frame(ValCh1=SArrayMQ21[,1,paste0("tau")])

# Plot ACF Plots for selected data-frames
P1 <- ggAcf(D_Workb1$ValCh1, lag.max = 50)+
  labs(x = "Lag", y = expression(beta[1])) +
  ggtitle(NULL)+
  theme_bw()
P2 <- ggAcf(D_Workb2$ValCh1, lag.max = 50)+
  labs(x = "Lag", y = expression(beta[2])) +
  ggtitle(NULL)+
  theme_bw()
P3 <- ggAcf(D_Worktau$ValCh1, lag.max = 50)+
  labs(x = "Lag", y = expression(tau)) +
  ggtitle(NULL)+
  theme_bw()

# Plot all Graphs in a same frame
ggpubr::ggarrange(P1+theme(axis.title.x=element_blank(),
                           legend.position="none"), 
                  P2+theme(axis.title.x=element_blank(),
                           legend.position="none"),
                  P3+theme(axis.title.x=element_blank(),
                           legend.position="none"), ncol = 2, nrow = 2)
# Remove Working Variable to free memory
remove(D_Workb1, D_Workb2, D_Worktau )

# MODEL Q2_2
cat("### Model 2 -
model{
  for (i in 1:N) {
    y[i]~dnorm(mu[i],tau)
    mu[i]<- b[1] + b[2]*(x[i]-31)+ b[3]*pow((x[i]-31),2)
 
    # Get the residuals for the observed value
    res[i] <- (y[i]-mu[i])           # Estimate the residual for this model  - Item (1)
    stdres[i] <- res[i]*sqrt(tau)    # Calculate the standard residuals - Item (2)
    
    dev1.obs[i] <- pow(res[i],2)
    dev2.obs[i] <- pow(stdres[i],2)
    
    # Get a replicated Sample - sample of the predictive distribution
    y.rep[i]~dnorm(mu[i],tau)
    p.smaller[i] <- step(y[i]-y.rep[i])  # Check to see the probability of getting extreme value - Item (3)
    
    # Residual and Moments replicated - this gives the predicted distribution for these values
    res.rep[i] <- y.rep[i]-mu[i]
    stdres.rep[i] <- res.rep[i]*sqrt(tau)
    
    dev1.rep[i] <- pow(res.rep[i],2)
    dev2.rep[i] <- pow(stdres.rep[i],2)
    
    # Likelihood for each observed and replicated data
    # note: Need to know the density function of the probability model
    
    loglike[i] <- (0.5)*log(tau/6.283185) + (-0.5)*tau*pow((y[i]-mu[i]),2)
    loglike.rep[i]<-  (0.5)*log(tau/6.283185) + (-0.5)*tau*pow((y.rep[i]-mu[i]),2)
    
    p.inv[i]<- 1/exp(loglike[i])      #  This is to find the predictive ordinate of the observations
  }
  
  # Prior definitions
  b[1]~dnorm(0,.000001)
  b[2]~dnorm(0,.000001)
  b[3]~dnorm(0,.01)
  tau~dgamma(.0001,.0001)

  # Summing Diagnostics Values
  chidev1.obs <- sum(dev1.obs[])
  chidev2.obs <- sum(dev2.obs[])

  chidev1.rep <- sum(dev1.rep[] )
  chidev2.rep <- sum(dev2.rep[] )

  chidev1.pval<-step(chidev1.obs-chidev1.rep)
  chidev2.pval<-step(chidev2.obs-chidev2.rep)

  #   Deviance statistic
  dev<-   -2*sum(loglike[])
  dev.rep <-  -2*sum(loglike.rep[])
  dev.pval<-step(dev-dev.rep)

}", file="HarvestModQ22.txt")

# Setup Parameters
paramsMQ22 <- c("b", "tau",          # the rest are for the model checking
                "mu", "res", "stdres", "res.rep", "stdres.rep", "p.smaller",
                "p.inv", "chidev1.pval", "chidev2.pval", "chidev1.obs", "chidev2.obs",
                "chidev1.rep", "chidev2.rep", "dev", "dev.rep", "dev.pval")

bugMQ22.dat=list("x", "y", "N")  # what variable you need in the model

# Setup Initial Values
initMQ22.fun=function(){ 
  list(b=runif(3,-.8,-.2), tau=runif(1,.2,.8), y.rep=rnorm(N)) 
}

#### Run Openbugs Simulation
set.seed(3508)
attach(HarvestingData)
HarvestingMQ22=bugs(bugMQ22.dat, initMQ22.fun, paramsMQ22, model.file="HarvestModQ22.txt",
                    n.chains=NChain, n.iter=NSim, n.burnin=Burnin, n.thin=NThin , debug=FALSE)
detach(HarvestingData)


# Get Simulation from OpenBugs
SArrayMQ22= HarvestingMQ22$sims.array   # Data Arrays
vname=attr(SArrayMQ22,"dimnames")[3][[1]]  # Variable Names
MQ22.coda <- as.mcmc.list(HarvestingMQ22)


# Print Summary statistics of parameters oif Interest
RNames <- rownames(HarvestingMQ22[["summary"]][,c(1:3,5,7)]) # List of parameters
df_Prt <- data.frame(Parameter = RNames)
df_Prt <- cbind(df_Prt, as_tibble(HarvestingMQ22[["summary"]][,c(1:3,5,7)]))
df_Prt %>%
  filter(substr(Parameter,1,2)=="b[" |
           substr(Parameter,1,3)=="tau") %>% 
  kbl(booktabs = TRUE, digits = 4, 
      caption = "OpenBugs Summary - Harvesting Date vs. Yield (Model-2)") %>% 
  kable_styling(latex_options = "striped")

# Plot TracePlots * b[1], b[2], b[3] and tau *

# Sampling 5000 points to generate "thinner" traceplots
set.seed(792)

L <- NSim-Burnin

S <- if (Sz < L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe
D_Workb1<- data.frame(ValCh1=SArrayMQ22[S,1,paste0("b[1]")],
                      ValCh2=SArrayMQ22[S,2,paste0("b[1]")],
                      ValCh3=SArrayMQ22[S,3,paste0("b[1]")])

D_Workb2 <- data.frame(ValCh1=SArrayMQ22[S,1,paste0("b[2]")],
                       ValCh2=SArrayMQ22[S,2,paste0("b[2]")],
                       ValCh3=SArrayMQ22[S,3,paste0("b[2]")])

D_Workb3 <- data.frame(ValCh1=SArrayMQ22[S,1,paste0("b[3]")],
                       ValCh2=SArrayMQ22[S,2,paste0("b[3]")],
                       ValCh3=SArrayMQ22[S,3,paste0("b[3]")])

D_Worktau <- data.frame(ValCh1=SArrayMQ22[S,1,paste0("tau")],
                        ValCh2=SArrayMQ22[S,2,paste0("tau")],
                        ValCh3=SArrayMQ22[S,3,paste0("tau")])

# Plot Traceplots for selected data-frames
P1 <- D_Workb1 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta[1])) +
  theme_bw()
P2 <- D_Workb2 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta[2])) +
  theme_bw()
P3 <- D_Workb3 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta[3])) +
  theme_bw()
P4 <- D_Worktau %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(tau)) +
  theme_bw()

# Plot all Graphs in a same frame
ggpubr::ggarrange(P1+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"), 
                  P2+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"),
                  P3+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"),
                  P4+theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank(),
                           legend.position="none"), ncol = 2, nrow = 2)

# Remove Working Variable to free memory
remove(D_Workb1, D_Workb2, D_Workb3, D_Worktau )

# Plot ACF Plots * b[1], b[2], b[3] and tau *

# Build working dataframe
D_Workb1 <- data.frame(ValCh1=SArrayMQ22[,1,paste0("b[1]")])
D_Workb2 <- data.frame(ValCh1=SArrayMQ22[,1,paste0("b[2]")])
D_Workb3 <- data.frame(ValCh1=SArrayMQ22[,1,paste0("b[3]")])
D_Worktau <- data.frame(ValCh1=SArrayMQ22[,1,paste0("tau")])

# Plot ACF Plots for selected data-frames
P1 <- ggAcf(D_Workb1$ValCh1, lag.max = 50)+
  labs(x = "Lag", y = expression(beta[1])) +
  ggtitle(NULL)+
  theme_bw()
P2 <- ggAcf(D_Workb2$ValCh1, lag.max = 50)+
  labs(x = "Lag", y = expression(beta[2])) +
  ggtitle(NULL)+
  theme_bw()
P3 <- ggAcf(D_Workb3$ValCh1, lag.max = 50)+
  labs(x = "Lag", y = expression(beta[3])) +
  ggtitle(NULL)+
  theme_bw()
P4 <- ggAcf(D_Worktau$ValCh1, lag.max = 50)+
  labs(x = "Lag", y = expression(tau)) +
  ggtitle(NULL)+
  theme_bw()

# Plot all Graphs in a same frame
ggpubr::ggarrange(P1+theme(axis.title.x=element_blank(),
                           legend.position="none"), 
                  P2+theme(axis.title.x=element_blank(),
                           legend.position="none"),
                  P3+theme(axis.title.x=element_blank(),
                           legend.position="none"),
                  P4+theme(axis.title.x=element_blank(),
                           legend.position="none"), ncol = 2, nrow = 2)
# Remove Working Variable to free memory
remove(D_Workb1, D_Workb2, D_Workb3, D_Worktau )

# Routines to Calculate DIC, Deviance of Model - Author: Michael Escobar (Thanks!)

# OBS>:- Flow remodeled to document each step

# Calculates the Deviance by hand
devNormFunc <- function(beta0, beta1, beta2, tau, x, y) {
  mu <- beta0 + beta1*x + beta2*x^2
  return(-2*sum(log(dnorm(y,mu,1/sqrt(tau)))))
}

CalcModelStats <- function (Model, TxtModel) {
  
  # STEP 1 - Calculate the Likelihood and respective Inverse
  M1.pinv <- Model$mean$p.inv       # Get 1/p(x)
  PLogLI <- -1*log(M1.pinv)                  # Get p(x)
  
  # Collects Inv.p(x) Likelihood, P[Log(Likelihood)]
  TbWork1 <- cbind(M1.pinv,1/M1.pinv,PLogLI)
  colnames(TbWork1) <- c("p.inv","p(x)","PLogLI");#TbWork1
  Pseudom2logL= -2*sum(PLogLI);#Pseudom2logL
  #cat("\n-------------------\n\n")
  
  # Collect Residuals (TbWork2)
  TbWork2<-cbind(Model$mean$res,
                 t(apply(Model$sims.list$res.rep,2,
                         function(x){c(quantile(x,probs=c(0.025,.975)),mean(x),sd(x))})))
  colnames(TbWork2)=c("res","2.5%","97.5%","mean","SD");#TbWork2
  #cat("\n-------------------\n\n")
  
  # Collect Std.Residuals (TbWork3)
  TbWork3<-cbind(Model$mean$stdres,
                 t(apply(Model$sims.list$stdres.rep,2,
                         function(x){c(quantile(x,probs=c(0.025,.975)),mean(x),sd(x))})))
  colnames(TbWork3)=c("stdres","2.5%","97.5%","mean","SD");#TbWork3
  #cat("\n-------------------\n\n")
  
  # Collect Deviance (TbWork4)
  M1.devrep<-Model$sims.list$dev.rep
  TbWork4<-c(Model$mean$dev,
             quantile(M1.devrep,probs=c(0.025,.975)),
             mean(M1.devrep),sd(M1.devrep))
  names(TbWork4)=c("Deviance","2.5%","97.5%","mean","SD");#TbWork4
  #cat("\n-------------------\n\n")
  
  # Collect Chi-Deviance (TbWork5)
  M1.chidev2rep<-Model$sims.list$chidev2.rep
  TbWork5<-c(Model$mean$chidev2.obs,
             quantile(M1.chidev2rep,probs=c(0.025,.975)),
             mean(M1.chidev2rep),sd(M1.chidev2rep))
  names(TbWork5)=c("Chi-Dev2","2.5%","97.5%","mean","SD");#TbWork5
  #cat("\n-------------------\n\n")
  
  # Calculates Non-Calibrated ("pval-stats") - REVIEW - 
  #apply(Model$sims.list$p.smaller,2,mean)
  #Model$mean$chidev2.pval
  #Model$mean$dev.pval
  #cat("\n-------------------\n\n")
  
  # Comparing intrinsic and self calculated value for Deviance:
  M1.dev <-Model$sims.list$dev
  M1.deviance <-Model$sims.list$deviance
  TbWork6 <-rbind(
    quantile(M1.dev,probs=c(0.025,.25,.5,.75,.975)),
    quantile(M1.deviance,probs=c(0.025,.25,.5,.75,.975)) )
  rownames(TbWork6)=c("SelfProgramed:","Openbugs Made:");# TbWork6    
  # cat("\n-------------------\n\n")
  
  # c(Model$mean$deviance,Model$mean$dev)    
  Dbar<-mean(M1.dev);#Dbar   # Calculated - Deviance 
  pd2<-0.5*var(M1.dev);#pd2
  
  
  # Setup Parameters for Model
  beta0Bar<- Model$mean$b[1]
  beta1Bar<- Model$mean$b[2]
  beta2Bar<- ifelse (dim(Model$sims.list$b)[2]==3,Model$mean$b[3],0)
  tauBar <- Model$mean$tau
  Dhat <- devNormFunc(beta0Bar, beta1Bar, beta2Bar, tauBar, 
                      HarvestingData$x,HarvestingData$y);#Dhat 
  pd1<-Dbar-Dhat;#pd1
  
  DIC1<-Dbar+pd1; #DIC1
  DIC2<-Dbar+pd2; #DIC2
  
  
  cat("\n------- Goodness of Fit Statistics for",TxtModel," --------")
  cat("\nDeviance: \n");
  print(TbWork6, dig=4)
  cat("\nDbar: ",format(Dbar, digits = 2, nsmall = 4))
  # cat("\nDhat: ",format(Dhat, digits = 2, nsmall = 4))
  cat("\npD:  ",format(pd2, digits = 2, nsmall = 4))
  cat("\nDIC: ",format(DIC2, digits = 2, nsmall = 4))
  cat("\n--------------------------------------------------------\n")
}


CalcModelStats(HarvestingMQ21,"Model-1")
CalcModelStats(HarvestingMQ22,"Model-2")

# Setup for Kuo-Mallick Model Variables - LAST VERSION

# Mean and Variance of Y
mu_y <- mean(HarvestingData$y)
sd_y <- sd(HarvestingData$y)

# Mean and Variance of X
mu_x <- mean(HarvestingData$x)
sd_x <- sd(HarvestingData$x)

# Mean and Variance of X^2
mu_x2 <- mean(HarvestingData$x^2)
sd_x2 <- sd(HarvestingData$x^2)

pp <- 0.5

# OBS:- I tested 'tau0' with several values and all decided in favor of M2 (which makes perfectly sense!) 
#       I ended up choosing a value inside the interval of 1/16<tau0<4
tau0 <- 0.3

# Standardizing Variables X and Y
sx <- (HarvestingData$x-mu_x)/sd_x          # Standardized X
sx2 <- (HarvestingData$x^2-mu_x2)/sd_x2     # Standardized X^2
sy <- (HarvestingData$y-mu_y)/sd_y          # Standardized y

cat("### Kuo-Mallick Model Selection
model{
  for (i in 1:N) {
    sy[i]~dnorm(mu[i],tau)
    mu[i]<- delta[1]*beta[1]*sx[i] + delta[2]*beta[2]*sx2[i]
  }
  
  for (ix in 1:2) {
    beta[ix]~dnorm(0,tau0)
  }
  
  tau~dgamma(.5,.01) 
  
  # Prior Distribution of delta[]
  for(k in 1:2){
    delta[k]~dbern(pp)
  }
  
  # Model Deltas
  for(i1 in 1:2){
    for(i2 in 1:2){
      mod[i1,i2]<-equals((2-i1),delta[1])*equals((2-i2),delta[2])
     }
  }

}", file="Harvest_KM.txt")

dataDiag_KM<-list("sx", "sy", "mu_y", "sd_y", "N",
                  "mu_x", "sd_x", "pp", "tau0")

initsDiag_KM<-function(){ list(beta=rnorm(2), tau=runif(1,.5,1))}

paramDiag_KM<-c("beta", "tau", "delta", "mod")


set.seed(9212)
attach(HarvestingData)

HarvestDiag_KM<-bugs(dataDiag_KM,initsDiag_KM, paramDiag_KM, model.file="Harvest_KM.txt",
                     n.chains=1, n.iter=20000, n.burnin=5000,
                     n.thin=10)
detach(HarvestingData)


# Calculus of abeta[]'s
abeta1 <- HarvestDiag_KM$sims.array[,,paste0("beta[1]")]*sd_x/sd_y
abeta2 <- HarvestDiag_KM$sims.array[,,paste0("beta[2]")]*sd_x/sd_y

t_abeta1 <- c(mean(abeta1),sd(abeta1),quantile(abeta1,c(.025,.5,.975)))
t_abeta2 <- c(mean(abeta2),sd(abeta2),quantile(abeta2,c(.025,.5,.975)))

abeta <- rbind(t_abeta1,t_abeta2)
colnames(abeta) <- c("mean","sd","2.5%","50%","97.5%")
rownames(abeta) <- c("abeta[1]","abeta[2]")

# Print Summary statistics of parameters oif Interest
RNames <- rownames(HarvestDiag_KM[["summary"]][,c(1:3,5,7)]) # List of parameters
df_KM <- data.frame(Parameter = RNames)
df_KM <- cbind(df_KM, as_tibble(HarvestDiag_KM[["summary"]][,c(1:3,5,7)]))
df_KM %>%
  kbl(booktabs = TRUE, digits = 4,
      caption = "OpenBugs Summary - Kuo-Mallik Model Selection Statistics") %>% 
  kable_styling(latex_options = "striped")

# Jeffrey's Criteria for Bayes Factor
JeffCrit <- function(BF) {
  if (BF>=0.0 && BF<=.5) return("Not worth more than a bare mention")
  else{
    if (BF>.5 && BF<=1.0) return("Substantial")
    else {
      if (BF>1.0 && BF<=2.0) return("Strong")
      else return("Decisive")
    }
  }
}

# Kass, Raftery's Criteria for Bayes Factor
KassRafteryCrit <- function(BF) {
  if (BF>=0.0 && BF<=2.0) return("Not worth more than a bare mention")
  else{
    if (BF>2.0 && BF<=6.0) return("Positive")
    else {
      if (BF>6.0 && BF<=10.0) return("Strong")
      else return("Decisive")
    }
  }
}

Kuo_MallickCrit <- function(PModA, PModB, TxtModA, TxtModB) {
  BF_KM <- PModA/PModB
  # Print Kuo-Mallik Bayes Factor Analysis for Model-A over Model-B
  cat("\n-------- Kuo-Mallik Bayes Factor Analysis (",TxtModA,"/",TxtModB,") ---------")
  cat("\n\nPosterior Probability",TxtModA,"| Data: ",format(PModA, digits = 2, nsmall = 4))
  cat("\n\nPosterior Probability",TxtModB,"| Data: ",format(PModB, digits = 2, nsmall = 4))
  cat("\n\n>>> Bayes Factor : ",format(BF_KM, digits = 2, nsmall = 4))
  cat("\n\n>>> Jeffrey's Evidence of",TxtModA,"over",TxtModB,": ",JeffCrit(log10(BF_KM)))
  cat("\n\n>>> Kass/Raftery's Evidence of",TxtModA,"over",TxtModB,": ",KassRafteryCrit(2*log(BF_KM)))
  cat("\n---------------------------------------------------------------\n")
}


# Calculation of BF - It is equivalent to compare
#    mod[1,1] (Model-2) agains mod[1,2] (Model-1)
PM1 <- df_KM[which(df_KM$Parameter=="mod[1,2]"),"mean"]#;PM1  # This is Probability of Model-1
PM2 <- df_KM[which(df_KM$Parameter=="mod[1,1]"),"mean"]#;PM2  # This is Probability of Model-2

Kuo_MallickCrit (PM2, PM1, "M2", "M1")

# Item (1) - Print Residuals
RNames <- rownames(HarvestingMQ21[["summary"]][,c(1:3,5,7)]) # List of parameters
df_Prt <- data.frame(Parameter = RNames)
df_Prt <- cbind(df_Prt, as_tibble(HarvestingMQ21[["summary"]][,c(1:3,5,7)]))
df_Prt %>%
  filter(substr(Parameter,1,4)=="res[") %>% 
  kbl(booktabs = TRUE, digits = 4,
      caption = "OpenBugs Summary - Residuals for Model-1") %>% 
  kable_styling(latex_options = "striped")

# Item (2) - Print Standardized Residuals
RNames <- rownames(HarvestingMQ21[["summary"]][,c(1:3,5,7)]) # List of parameters
df_Prt <- data.frame(Parameter = RNames)
df_Prt <- cbind(df_Prt, as_tibble(HarvestingMQ21[["summary"]][,c(1:3,5,7)]))
df_Prt %>%
  filter(substr(Parameter,1,7)=="stdres[") %>% 
  kbl(booktabs = TRUE, digits = 4,
      caption = "OpenBugs Summary - Standardized Residuals for Model-1") %>% 
  kable_styling(latex_options = "striped")

# Item (3) - Print Chance of Getting Extreme Observations
RNames <- rownames(HarvestingMQ21[["summary"]][,c(1:2)]) # List of parameters
df_Prt <- data.frame(Parameter = RNames)
df_Prt <- cbind(df_Prt, as_tibble(HarvestingMQ21[["summary"]][,c(1:2)]))
df_Prt %>%
  filter(substr(Parameter,1,10)=="p.smaller[") %>% 
  kbl(booktabs = TRUE, digits = 4,
      caption = "OpenBugs Summary - Chance of Extreme Observations for Model-1") %>% 
  kable_styling(latex_options = "striped")

# Get Distribution of Residuals/Std Residuals under Predictive Distribution
SArrayDiagM1= HarvestingMQ21$sims.array      # Data Arrays
vname=attr(SArrayDiagM1,"dimnames")[3][[1]]  # Variable Names


# Plot Sample densities to Compare Distributions of Residuals

# Build working dataframe
D_WorkRes <- as_tibble(data.frame(LvlStd=as.vector(
  SArrayDiagM1[,1,c(paste0("res.rep[", 1:16, "]"))]), Chain = "C1"))
D_WorkstdRes <- as_tibble(data.frame(LvlStd=as.vector(
  SArrayDiagM1[,1,c(paste0("stdres.rep[", 1:16, "]"))]), Chain = "C1"))

# Ploting the Sensities
P1 <- D_WorkRes %>%
  ggplot(aes(x = LvlStd))+
  geom_density(aes(x=LvlStd, group=Chain, colour=Chain), size=0.8)+ 
  labs(x = "Residuals", y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()
P2 <- D_WorkstdRes %>%
  ggplot(aes(x = LvlStd))+
  geom_density(aes(x=LvlStd, group=Chain, colour=Chain), size=0.8)+ 
  labs(x = "Std.Residuals", y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

ggpubr::ggarrange(P1+theme(legend.position="none"), 
                  P2+theme(legend.position="none"), ncol = 1, nrow = 2)
remove(D_WorkRes,D_WorkstdRes)

# Sample points for Residuals
set.seed(792)
LR <-  length(SArrayDiagM1[,1,c(paste0("res.rep[", 1:16, "]"))])
S <- sort(sample(1:LR, 30000, replace = FALSE))

# Build working dataframe
D_WorkRes <- rbind(
  data.frame(
    Iter = c(1:length(S)),
    Res=as.vector(
      SArrayDiagM1[,1,c(paste0("res.rep[", 1:16, "]"))][S]), Model = "M1"),
  data.frame(
    Iter = c(1:length(S)),
    Res=as.vector(
      SArrayDiagM2[,1,c(paste0("res.rep[", 1:16, "]"))][S]), Model = "M2"))

# Ploting Residuals
D_WorkRes %>%
  ggplot(aes(x=Iter))+
  geom_point(aes(y=Res, colour=Model), size=0.6)+ 
  labs(x = "Iteration", y = "Residual") +
  scale_color_manual(values = cbp2)+
  theme_bw()

remove(D_WorkRes)

# Sample points for Residuals
set.seed(792)
LR <-  length(SArrayDiagM1[,1,c(paste0("res.rep[", 1:16, "]"))])
S <- sort(sample(1:LR, 30000, replace = FALSE))

# Build working dataframe
D_WorkRes <- rbind(
  data.frame(
    Iter = c(1:length(S)),
    Res=as.vector(
      SArrayDiagM1[,1,c(paste0("res.rep[", 1:16, "]"))][S]), Model = "M1"),
  data.frame(
    Iter = c(1:length(S)),
    Res=as.vector(
      SArrayDiagM2[,1,c(paste0("res.rep[", 1:16, "]"))][S]), Model = "M2"))

# Calculate Absolute Error for each model
D_WorkAbsErr <- rbind(
  data.frame(
    Iter = c(1:length(S)),
    AbsErr=as.vector(SArrayDiagM1[,1,c(paste0("res.rep[", 1:16, "]"))][S])/LR, Model = "M1"),
  data.frame(
    Iter = c(1:length(S)),
    AbsErr=as.vector(SArrayDiagM2[,1,c(paste0("res.rep[", 1:16, "]"))][S])/LR, Model = "M2"))

# Calculate Relative Error for each model
D_WorkRelErr <- rbind(
  data.frame(
    Iter = c(1:length(S)),
    RelErr=as.vector(SArrayDiagM1[,1,c(paste0("res.rep[", 1:16, "]"))][S]/
                       SArrayDiagM1[,1,c(paste0("mu[", 1:16, "]"))][S]), Model = "M1"),
  data.frame(
    Iter = c(1:length(S)),
    RelErr=as.vector(SArrayDiagM2[,1,c(paste0("res.rep[", 1:16, "]"))][S]/
                       SArrayDiagM2[,1,c(paste0("mu[", 1:16, "]"))][S]), Model = "M2"))

# Ploting Residuals
D_WorkRes %>%
  ggplot(aes(x=Iter))+
  geom_point(aes(y=Res, colour=Model), size=0.6)+ 
  labs(x = "Iteration", y = "Residual") +
  scale_color_manual(values = cbp2)+
  theme_bw()

# Ploting Absolute error
P2 <- D_WorkAbsErr %>%
  ggplot(aes(x=Iter))+
  geom_point(aes(y=AbsErr, colour=Model), size=0.6)+ 
  labs(x = "Iteration", y = "Absolute Error") +
  scale_color_manual(values = cbp2)+
  theme_bw()

# Ploting Relative error
P3 <- D_WorkRelErr %>%
  ggplot(aes(x=Iter))+
  geom_point(aes(y=RelErr, colour=Model), size=0.6)+ 
  labs(x = "Iteration", y = "Relative Error") +
  scale_color_manual(values = cbp2)+
  theme_bw()

ggpubr::ggarrange(P2+theme(axis.title.x=element_blank()), 
                  P3+theme(axis.title.x=element_blank()), ncol = 1, nrow = 2)


remove(D_WorkAbsErr, D_WorkRelErr, D_WorkRes)

