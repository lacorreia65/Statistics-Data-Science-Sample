library(tidyverse)
library(R2OpenBUGS)
library(kableExtra)
library(coda)


library(ggplot2)
library(forecast)

# The palette with black - Used in Graphs with :
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp3 <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", 
          "#C3D7A4", "#52854C", "#4E84C4", "#293352")

# Triangular Distribution
g=function(x){
  (x>0)*(x<1)*((x<=0.5)*4*x+ (x>0.5)*(4-4*x))
}

# Set parameters
N = 1000

set.seed(123)

# Generating two independent samples U1 and U2 from Uniform[0,1]
U1 <- runif(N)
U2 <- runif(N)

# Calculating the exact distribution of Z:= (U1+U2)/2
Z <- (U1+U2)/2

# Estimate E(X) and E(Y) Considerin that Z ~Triang(g)
EZ <- mean(Z)
VarZ <- var(Z)
SE_Z <- sd(Z)/sqrt(N)

# Print Statistics for Approx. Z
cat("\n---- Summary Statistics for approx. Z=(U1+U2)/2 ---------------")
cat("\nE(X) : ",format(EZ, digits = 2, nsmall = 4)) 
cat("\nS.E.: ",format(SE_Z, digits = 5, nsmall = 4))
CR_Z_up <- EZ+2*SE_Z; CR_Z_lo <- EZ-2*SE_Z; 
cat("\n\n95% Credible Region for true-Mean (", 
    format(CR_Z_lo, digits = 6, nsmall = 4), ",",
    format(CR_Z_up, digits = 6, nsmall = 4),")\n")
cat("\nVar(X) : ",format(VarZ, digits = 5, nsmall = 4)) 
cat("\n---------------------------------------------------------------\n")


dfz <- data.frame (ZS = Z)
ggplot(data=dfz, aes(x=ZS))+
  geom_histogram(color="black", fill="white", bins=40)+
  scale_color_manual(values = cbp2)+
  labs(x = "g(.) sample", y = "Histogram of Z")+
  theme_bw()
remove(dfz)

# Item (b) - here h_1(X) = X => h_1(U3)=U3
set.seed(725)
U3 <- runif(N)

# Estimate E(X) via Importance Sampling 
EX_IS <- (1/N)*sum((U3*g(U3))/(1))

# here h_2(X) = X^2 => h_1(U3)=U3^2
EX_IS2 <- (1/N)*sum((U3^2*g(U3))/(1))
VarX_IS <- EX_IS2-EX_IS^2
SE_IS <- sd((U3*g(U3))/(1))/sqrt(N)

# Print Statistics for Importance Sampling
cat("\n---- Summary Statistics for Importance Sampling Algorithm -----")
cat("\nE(X) : ",format(EX_IS, digits = 2, nsmall = 4)) 
cat("\nS.E.: ",format(SE_IS, digits = 5, nsmall = 4))
CR_IS_up <- EX_IS+2*SE_IS; CR_IS_lo <- EX_IS-2*SE_IS; 
cat("\n\n95% Credible Region for true-Mean (", 
    format(CR_IS_lo, digits = 6, nsmall = 4), ",",
    format(CR_IS_up, digits = 6, nsmall = 4),")\n")
cat("\nVar(X) : ",format(VarX_IS, digits = 5, nsmall = 4)) 
cat("\n---------------------------------------------------------------\n")


# Adapted from Robert, Casella - pp.66
xx=g(U3)
mxx <- mean(xx)
estint=cumsum(xx)/(1:N)
esterr=sqrt(cumsum((xx-estint)^2))/(1:N)
dxx <- data.frame(yest = estint,
                  low = estint-2*esterr,
                  up = estint+2*esterr)
ggplot(data=dxx, aes(x=c(1:N)))+
  geom_line(aes(y=yest), size=0.5)+
  geom_line(aes(y=low), color="#E69F00", size=0.5)+
  geom_line(aes(y=up), color="#E69F00", size=0.5)+
  scale_color_manual(values = cbp2)+
  labs(x = "g(.) sample", y = "Density Integration Estimate")+
  ylim(mxx+20*c(-esterr[N],esterr[N]))+
  theme_bw()
remove(xx, mxx, estint, esterr, dxx)

# Initialize Variables
set.seed(923)

xrange <- 1 # Only values in range from 0 to 'xrange' are of interest
M <- 2      # Upper Limit

acc <- rej <- 0
y <- rep(0,N)   # vector of sampled data

while (acc <= N) {
  # Propose a 'x' on support of g
  x <- runif(1, min = 0, max = xrange)
  
  # Generate Accept/Rejection criteria for each fitted value
  u <- runif(1)
  
  # Maximum of value for distribution 'g'
  if (u <= g(x)/(M)) {
    acc <- acc + 1
    y[acc] <- x
  }
  else rej <- rej + 1
}

# Estimate E(X) and Var(X) using the sample of accepted values
EX_AcpRej <- mean(y)
VarX_AcpRej <- var(y)
SE_AcpRej <- sd(y)/sqrt(N)

# Print Statistics for Accept/Reject
cat("\n---- Summary Statistics for Accept-Reject Algorithm -----")
cat("\nE(X) : ",format(EX_AcpRej, digits = 2, nsmall = 4)) 
cat("\nS.E.: ",format(SE_AcpRej, digits = 5, nsmall = 4))
CR_AcpRej_up <- EX_AcpRej+2*SE_AcpRej; CR_AcpRej_lo <- EX_AcpRej-2*SE_AcpRej; 
cat("\n\n95% Credible Region for true-Mean (", 
    format(CR_AcpRej_lo, digits = 6, nsmall = 4), ",",
    format(CR_AcpRej_up, digits = 6, nsmall = 4),")\n")
cat("\nVar(X) : ",format(VarX_AcpRej, digits = 5, nsmall = 4)) 
cat("\n\nRate of Acceptance = ", format(acc/(acc+rej), digits = 2, nsmall = 4))
cat("\n---------------------------------------------------------------\n")


dfy <- data.frame (YS = y)
ggplot(data=dfy, aes(x=YS))+
  geom_histogram(color="black", fill="white", bins=40)+
  scale_color_manual(values = cbp2)+
  labs(x = "g(.) sample", y = "Histogram of Y")+
  theme_bw()
remove(dfy)

#Setup Variables
set.seed(936)

alpha = function(x, y){ 
  min(1, g(x) / g(y))
}
x = rep(0, N)
acc <- 0

# Loop Sampling from the Chain
for(t in 2:N){
  ystar  <- runif(1, min=0, max=xrange)
  T <- runif(1)
  if( T <= alpha(ystar,x[t-1])) 
  {
    x[t] <- ystar
    acc <- acc+1
  } 
  else x[t] <- x[t-1]
}

# Print Rate of Acceptance
EX_MH <- mean(x)
Var_MH <- var(x)
SE_MH <- sd(x)/sqrt(N)

# Print Statistics for Metropolis-Hastings
cat("\n---- Summary Statistics for Metropolis-Hastings Algorithm -----")
cat("\nE(X) : ",format(EX_MH, digits = 2, nsmall = 4)) 
cat("\nS.E.: ",format(SE_MH, digits = 5, nsmall = 4))
CR_MH_up <- EX_MH+2*SE_MH; CR_MH_lo <- EX_MH-2*SE_MH; 
cat("\n\n95% Credible Region for true-Mean (", 
    format(CR_MH_lo, digits = 6, nsmall = 4), ",",
    format(CR_MH_up, digits = 6, nsmall = 4),")\n")
cat("\nVar(X) : ",format(Var_MH, digits = 5, nsmall = 4)) 
cat("\n\nRate of Acceptance = ", format(acc/N, digits = 2, nsmall = 4))
cat("\n---------------------------------------------------------------\n")


dfx <- data.frame (XS = x)
ggplot(data=dfx, aes(x=XS))+
  geom_histogram(color="black", fill="white", bins=40)+
  scale_color_manual(values = cbp2)+
  labs(x = "g(.) sample", y = "Histogram of X")+
  theme_bw()
remove(dfx)

# Prepare Data - Get delta[1-12]
# D3_Work <- as.data.frame(DiabetDrug_M1[["summary"]][4:15,c(1:3,5,7)])

# Generates Studies labels
StudyLst <- rapply(list(c("Z-Estimate", "Imp.Samp.","Accept/Reject","Metr.Hast.")), 
                   sprintf, fmt = "%10s", how = "replace")

# Collects C.R's from Summary report 
CR_Estim <- data.frame(ID = StudyLst[[1]], 
                       lower = c(CR_Z_lo, CR_IS_lo, CR_AcpRej_lo, CR_MH_lo), 
                       estim = c(EZ, EX_IS, EX_AcpRej, EX_MH), 
                       upper = c(CR_Z_up, CR_IS_up, CR_AcpRej_up, CR_MH_up))

# Organize data to generate graphs & analysis
P1 <- CR_Estim %>% 
  ggplot() +
  geom_point(aes(x=ID, y=estim, color=ID), size=2.5) +
  geom_errorbar(aes(x=ID, ymin=lower, ymax=upper, color=ID), width = 1) +
  xlab("Algorithm") +
  ylab(expression("Estimate / C.R.")) +
  scale_color_manual(values = cbp2)+
  theme_bw()

P1 + theme(axis.text.x=element_blank())

# Setup Data-set - Read Data
RatGrowth=read.table("Data/RatData11weeks.csv", header=TRUE, sep = ",")

# Setup Variables (ORIGINAL DATABASE)
N <- nrow(RatGrowth)  # Number of Items in data-set
ni <- length(unique(RatGrowth$IDinDose))  # No. of Rats on each dose level
nj <- ncol(RatGrowth)-2      # No. of Weeks of treatment
nk <- length(unique(RatGrowth$dose))      # No. of dose levels

ProdRun <- TRUE  # Control Variable for Testing/Production run

if (ProdRun) {
  # Production Setup OpenBugs running parameters
  NSim <- 30000    # No. of simulations for production
  NChain <- 3      # No. of chains for production
  NThin <- 8      # n.thin parameter for production
  Burnin <- 10000  # Burn-In parameter for production
  Sz <- 5000       # Size of samples for trace/acf plots
} else {
  # Testing Setup OpenBugs running parameters
  NSim <- 5000    # No. of simulations for production
  NChain <- 3      # No. of chains for production
  NThin <- 5      # n.thin parameter for production
  Burnin <- 1000  # Burn-In parameter for production
  Sz <- 1000       # Size of samples for trace/acf plots
  
}


# Printing the Data-Frame
RatGrowth %>%
  kbl(booktabs = TRUE, digits = 4, longtable = TRUE,
      caption = "Data - Growth in Rats under Treatment") %>% 
  kable_styling(latex_options = "striped")
 

# Model 0 - Standard Model

# Setup variables
y <- as.matrix(RatGrowth[,c(paste0("week",c(1:nj)))])
doselevel <- RatGrowth[,"dose"]
doselevel.c <- RatGrowth[,"dose"]-mean(RatGrowth[,"dose"])
week <- c(1:nj)
week.c <- c(1:nj)-mean(c(1:nj))

# Setup Model in OpenBugs
cat("
model{
  for (i in 1:N) {
    for (j in 1:nj) {
      y[i,j] ~dnorm(mu[i,j], tau.w)
      mu[i,j] <- beta0[i]+beta1[i]*week[j]
    }
    beta0[i]~dnorm(mu0[i], tau.b0)
    beta1[i]~dnorm(mu1[i], tau.b1)
    mu0[i] <- beta00+beta01*doselevel[i]
    mu1[i] <- beta10+beta11*doselevel[i]
  }
  beta00 ~ dnorm(100.0, 0.00001)
  beta01 ~ dnorm(0.0, 0.0001)
  beta10 ~ dnorm(0.0, 0.0001)
  beta11 ~ dnorm(0.0, 0.0001)
  
  s.y ~dunif(0.0, 250.0)
  s.b0 ~dunif(0.0, 250.0)
  s.b1 ~dunif(0.0, 250.0)

  tau.w <- pow(s.y, -2)
  tau.b0 <- pow(s.b0, -2)
  tau.b1 <- pow(s.b1, -2)
}", file="RatsGrowthM0.txt")
  
# Setup Parameters
paramsM0=c("tau.w", "beta0", "beta1",
           "beta00", "beta01", "beta10", "beta11", 
           "tau.b0", "tau.b1")

bugM0.dat=list("y", "week", "doselevel", "N", "nj")  # what variable you need in the model


# Setup Initial Values- Stochastic Components
set.seed(963)
initM0.fun=function(){ list(  
  beta0 = rnorm(N,50.0,10.0), 
  beta1 = rnorm(N,50.0,10.0),
  beta00 = rnorm(1,100.0, 0.00001),
  beta01 = rnorm(1,0.0, 0.0001),
  beta10 = rnorm(1,0.0, 0.0001),
  beta11 = rnorm(1,0.0, 0.0001),
  s.y = runif(1, 0.0, 250.0),  
  s.b0 = runif(1, 0.0, 250.0), 
  s.b1 = runif(1, 0.0, 250.0) 
  ) }

# Run Open Bugs - Parameters according with 'ProdRun' flag
set.seed(2602)
attach(RatGrowth)
RatGrowthM0=bugs(bugM0.dat, initM0.fun, paramsM0, model.file="RatsGrowthM0.txt",
                 n.chains=NChain, n.iter=NSim, n.burnin=Burnin, n.thin=NThin, debug=FALSE)
detach(RatGrowth)

# Get Simulation from OpenBugs
SArrayM0 <-  RatGrowthM0$sims.array   # Data Arrays
vname <- attr(SArrayM0,"dimnames")[3][[1]]  # Variable Names

# Get Summary statistics of parameters of Interest
RNames <- rownames(RatGrowthM0[["summary"]]) # List of parameters
df_SummaryM0 <- data.frame(Parameter = RNames)
df_SummaryM0 <- cbind(df_SummaryM0, as_tibble(RatGrowthM0[["summary"]]))
rownames(df_SummaryM0) <- RNames


# Plot TracePlots * BETA00, BETA01, BETA10, BETA11*

# Sampling 5000 points to generate "thinner" traceplots
set.seed(312)

L <- NSim-Burnin

S <- if (Sz < L) sort(sample(1:L, Sz, replace = FALSE)) else 1:Sz

# Build working dataframe
D_WorkBeta00 <- data.frame(ValCh1=SArrayM0[S,1,paste0("beta00")],
                           ValCh2=SArrayM0[S,2,paste0("beta00")],
                           ValCh3=SArrayM0[S,3,paste0("beta00")])

D_WorkBeta01 <- data.frame(ValCh1=SArrayM0[S,1,paste0("beta01")],
                           ValCh2=SArrayM0[S,2,paste0("beta01")],
                           ValCh3=SArrayM0[S,3,paste0("beta01")])

D_WorkBeta10 <- data.frame(ValCh1=SArrayM0[S,1,paste0("beta10")],
                           ValCh2=SArrayM0[S,2,paste0("beta10")],
                           ValCh3=SArrayM0[S,3,paste0("beta10")])

D_WorkBeta11 <- data.frame(ValCh1=SArrayM0[S,1,paste0("beta11")],
                           ValCh2=SArrayM0[S,2,paste0("beta11")],
                           ValCh3=SArrayM0[S,3,paste0("beta11")])

# Plot Traceplots for selected data-frames
P1 <- D_WorkBeta00 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta['00'])) +
  theme_bw()
P2 <- D_WorkBeta01 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta['01'])) +
  theme_bw()
P3 <- D_WorkBeta10 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta['10'])) +
  theme_bw()
P4 <- D_WorkBeta11 %>%
  ggplot(aes(seq(from=1,to=Sz)))+
  geom_line(aes(y=ValCh1, colour=1), size=0.8)+ 
  geom_line(aes(y=ValCh2, colour=2), size=0.8)+ 
  geom_line(aes(y=ValCh3, colour=3), size=0.8)+ 
  labs(y = expression(beta['11'])) +
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
remove(D_WorkBeta00, D_WorkBeta01, D_WorkBeta10, D_WorkBeta11 )

# Plot ACF Plots * BETA00, BETA01, BETA10, BETA11*

# Build working dataframe
D_WorkBeta00 <- data.frame(ValCh1=SArrayM0[,1,paste0("beta00")])
D_WorkBeta01 <- data.frame(ValCh1=SArrayM0[,1,paste0("beta01")])
D_WorkBeta10 <- data.frame(ValCh1=SArrayM0[,1,paste0("beta10")])
D_WorkBeta11 <- data.frame(ValCh1=SArrayM0[,1,paste0("beta11")])

# Plot ACF Plots for selected data-frames
P1 <- ggAcf(D_WorkBeta00$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(beta['00'])) +
  ggtitle(NULL)+
  theme_bw()
P2 <- ggAcf(D_WorkBeta01$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(beta['01'])) +
  ggtitle(NULL)+
  theme_bw()
P3 <- ggAcf(D_WorkBeta10$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(beta['10'])) +
  ggtitle(NULL)+
  theme_bw()
P4 <- ggAcf(D_WorkBeta11$ValCh1, lag.max = 100)+
  labs(x = "Lag", y = expression(beta['11'])) +
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
remove(D_WorkBeta00, D_WorkBeta01, D_WorkBeta10, D_WorkBeta11 )

# Print statistics of interest - Beta0i and Beta1i
df_Prt <- df_SummaryM0[c(paste0("beta0[",c(1:N),"]"),
                         paste0("beta1[",c(1:N),"]")),c(2:4,6,8)]

df_Prt %>%
  kbl(booktabs = TRUE, digits = 4, longtable = TRUE, 
      caption = "Summary - Beta0i and Beta1i for Rat Growth (Run-0)") %>% 
  kable_styling(latex_options = "striped")

# Print statistics of interest - beta00-beta11
df_Prt <- df_SummaryM0[c("beta00", "beta01", "beta10", "beta11"),c(2:4,6,8)]

df_Prt %>%
  kbl(booktabs = TRUE, digits = 4,
      caption = "Summary - beta00-beta11 for Rat Growth (Run-0)") %>% 
  kable_styling(latex_options = "striped")

# Print statistics of interest - Precisions
df_Prt <- df_SummaryM0[c("tau.w", "tau.b0", "tau.b1"),c(2:4,6,8)]

df_Prt %>%
  kbl(booktabs = TRUE, digits = 4,
      caption = "Summary - Precisions for Rat Growth (Run-0)") %>% 
  kable_styling(latex_options = "striped")

# Create Temporary data-frame
D1_Work <- RatGrowth
colnames(D1_Work) <- c("dose", "IDinDose", "week01","week02","week03","week04","week05",
                         "week06","week07","week08","week09","week10","week11")

P2 <- D1_Work %>% 
  pivot_longer(cols = week01:week11, names_to = "weekStr", values_to = "rweight") %>%
  mutate(weekNo = as.integer(substr(weekStr,5,length(weekStr)-5))) %>% 
  ggplot(aes(group=weekNo))+ 
  geom_boxplot(aes(x=weekStr, y=rweight)) +
  labs(y="Weight (in g)")+
  facet_grid(~dose)+
  theme_bw()
P2 + theme(axis.title.x=element_blank(),
           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
remove(D1_Work)


# Prepare Data-Frame for next steps
D1_Work <- rbind(data.frame(beta00=SArrayM0[,1,paste0("beta00")],
                            beta01=SArrayM0[,1,paste0("beta01")],
                            beta10=SArrayM0[,1,paste0("beta10")],
                            beta11=SArrayM0[,1,paste0("beta11")], Chain=factor(1)),
                 data.frame(beta00=SArrayM0[,2,paste0("beta00")],
                            beta01=SArrayM0[,2,paste0("beta01")],
                            beta10=SArrayM0[,2,paste0("beta10")],
                            beta11=SArrayM0[,2,paste0("beta11")], Chain=factor(2)),
                 data.frame(beta00=SArrayM0[,3,paste0("beta00")],
                            beta01=SArrayM0[,3,paste0("beta01")],
                            beta10=SArrayM0[,3,paste0("beta10")],
                            beta11=SArrayM0[,3,paste0("beta11")], Chain=factor(3)))

# Plot the Graphs of beta00, beta01, beta10, beta11
P1 <- D1_Work %>%
  ggplot(mapping = aes(x = beta00, group = Chain))+
  geom_density(aes(colour=Chain), size=0.5)+
  labs(x = expression(beta['00']), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

P2 <- D1_Work %>%
  ggplot(mapping = aes(x = beta01, group = Chain))+
  geom_density(aes(colour=Chain), size=0.5)+
  labs(x = expression(beta['01']), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

P3 <- D1_Work %>%
  ggplot(mapping = aes(x = beta10, group = Chain))+
  geom_density(aes(colour=Chain), size=0.5)+
  labs(x = expression(beta['10']), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

# Plot the Graphs of Tau
P4 <-D1_Work %>%
  ggplot(mapping = aes(x = beta11, group = Chain))+
  geom_density(aes(colour=Chain), size=0.5)+
  labs(x = expression(beta['11']), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

ggpubr::ggarrange(P1, P2, P3, P4, ncol = 2, nrow = 2)

# Remove Working Variable to free memory
remove(D1_Work)

# Prepare Data-Frame for next steps
D2_Work <- rbind(data.frame(tau.w=SArrayM0[,1,paste0("tau.w")],
                            tau.b0=SArrayM0[,1,paste0("tau.b0")],
                            tau.b1=SArrayM0[,1,paste0("tau.b1")], Chain=factor(1)),
                 data.frame(tau.w=SArrayM0[,2,paste0("tau.w")],
                            tau.b0=SArrayM0[,2,paste0("tau.b0")],
                            tau.b1=SArrayM0[,2,paste0("tau.b1")], Chain=factor(2)),
                 data.frame(tau.w=SArrayM0[,3,paste0("tau.w")],
                            tau.b0=SArrayM0[,3,paste0("tau.b0")],
                            tau.b1=SArrayM0[,3,paste0("tau.b1")], Chain=factor(3)))


# Plot the Graphs of TAU_W, TAU_Beta0, TAU_Beta1
P1 <- D2_Work %>%
  ggplot(mapping = aes(x = tau.w, group = Chain))+
  geom_density(aes(colour=Chain), size=0.5)+
  labs(x = expression(tau['w']), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

P2 <- D2_Work %>%
  ggplot(mapping = aes(x = tau.b0, group = Chain))+
  geom_density(aes(colour=Chain), size=0.5)+
  labs(x = expression(tau['b0']), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

# Plot the Graphs of Tau
P3 <-D2_Work %>%
  ggplot(mapping = aes(x = tau.b1, group = Chain))+
  geom_density(aes(colour=Chain), size=0.5)+
  labs(x = expression(tau['b1']), y = "Density") +
  scale_color_manual(values = cbp2)+
  theme_bw()

ggpubr::ggarrange(P1, P2, P3, ncol = 1, nrow = 3)

# Remove Working Variable to free memory
remove(D2_Work)
