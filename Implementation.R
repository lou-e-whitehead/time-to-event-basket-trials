################################################################################
# Code to implement Bayesian models proposed in manuscript 
# 'Bayesian borrowing in basket trials with time-to-event outcomes 
# Lou E. Whitehead, l.whitehead2@ncl.ac.uk
# 10th July 2024
################################################################################
################################################################################

################################################################################
# Required packages ############################################################
################################################################################
library(R2jags)
library(survival)

################################################################################
# Load veteran dataset (used in Motivating Example) from survival package ######
################################################################################
data(cancer, package="survival") 

################################################################################
# Format the data accordingly ##################################################
################################################################################
# extract relevant columns of the dataset into a new data.frame
veteran1 = data.frame('id'=1:length(veteran$trt), veteran$trt, veteran$celltype, veteran$time, veteran$status) 

# rename columns just for context / clarity
colnames(veteran1)<-c('id', 'trt', 'basket', 'time','status') 

# look at the first few rows of the data!
head(veteran1)

# set up data to pass to the models (using veteran1 dataframe)
t<-veteran1$time # time to death or time to censoring
is.na(t)<-veteran1$status==0 # assign NA to any censored survival times (in this dataset, censored times have status==0)
is.censored<-1-veteran1$status # is.censored = 1 if censored, 0 otherwise (required for model implementation via 'dinterval' distribution in JAGS [see section 9.2.4 of JAGS version 4.3.0 user manual for more details])
t.cen<-veteran1$time+veteran1$status  # Element i of t.cen is the censoring time for case i when this case
                                      # was censored. If the event was observed for case i then the “censoring
                                      # time” (t.cen) must be greater than the observed event time.
Trt = veteran1$trt # treatment assignment indicator; in the veteran dataset 1=std trt, 2=test trt)
Module <- as.integer(veteran1$basket) # hypothetical basket assignment based on cell type; squamous = 1, small cell = 2, adeno = 3, large = 4
k=length(unique(Module)) # number of baskets

nk <-vector() # number of patients in each basket
for (i in 1:k){
  nk[i] = length(Trt[Module == i])
}
totaln=sum(nk) # total number of patients in the veteran study

################################################################################
# Set working directory to where model files are stored ########################
################################################################################
setwd("C:/Users/yourfilepathhere")

################################################################################
# Pooled Bayesian analysis (estimate treatment efficacy in the entire cohort)
################################################################################
# Set prior hyperparameters
prior.beta0<-c(0,1e2) 
prior.beta1<-c(0,1e2)
a<-1.1
b<-1.1

# Set cutoff.theta. Posterior probability, pCat = P(theta<cutoff.theta)
cutoff.theta = 1

# Parameters to return - 'theta' = hazard ratio, pCat = Pr(theta < 1|Data)
parameters<- c('theta', 'pCat')

# run model
set.seed(1234) # set.seed for reproducibility
for (i in 1:k){
  pooled_data<-list(n=totaln,t=t,is.censored=is.censored,t.cen=t.cen,Trt=Trt, 
                       prior.beta0=prior.beta0, prior.beta1=prior.beta1,a=a,b=b,
                       cutoff.theta=cutoff.theta)
  
  pooled_jags <- jags(data=pooled_data, 
                          parameters.to.save=parameters, 
                          model.file="Strat-W.txt", # use stratified
                          n.chains = 4,
                          n.burnin = 2000, 
                          n.iter = 4000)
}

################################################################################
# Stratified Bayesian analysis (Strat-W)
# (estimate treatment efficacy in each basket with no borrowing between baskets)
################################################################################
strat_data = list()
strat_jags = list()
set.seed(1234)
for (i in 1:k){
  strat_data[[i]]<-list(n=nk[i],t=t[Module == i],is.censored=is.censored[Module == i],t.cen=t.cen[Module == i],Trt=Trt[Module == i], 
                       prior.beta0=prior.beta0, prior.beta1=prior.beta1,a=a,b=b,
                       cutoff.theta=cutoff.theta) # 
  
  strat_jags[[i]] <- jags(data=strat_data[[i]], 
                          parameters.to.save=parameters, 
                          model.file="Strat-W.txt",
                          n.chains = 4,
                          n.burnin = 2000, 
                          n.iter = 4000)
}
################################################################################
# Bayesian hierarchical model (BHM-W) 
# (estimate treatment efficacy in each basket with borrowing between baskets
# using standard Bayesian hierarchical model)
################################################################################
# Set prior hyperparameters for BHM
prior.B0<-c(0,1e2)  
prior.B1<-c(0,1e2)
a<-1.1
b<-1.1
prior.sig.HN0<-1
prior.sig.HN1<-1

# Set cutoff.theta. Posterior probability, pCat = P(theta<cutoff.theta)
cutoff.theta = 1

# Parameters to return
parameters<- c('theta', 'pCat')

bhm_data<-list(totaln=totaln, Module=Module, k=k, t=t,is.censored=is.censored,t.cen=t.cen,Trt=Trt, 
              prior.B0=prior.B0, prior.B1=prior.B1,a=a,b=b,prior.sig.HN0=prior.sig.HN0, prior.sig.HN1=prior.sig.HN1,
              cutoff.theta=cutoff.theta) 
set.seed(1234)
bhm_jags <- jags(data=bhm_data, 
                  parameters.to.save=parameters, 
                  model.file="BHM-W.txt",
                  n.chains = 4,
                  n.iter = 4000,
                  n.burnin = 2000)

################################################################################
# Exchangeable-Nonexchangeable model (EXNEX-W) 
# (estimate treatment efficacy in each basket with borrowing between baskets
# using EXNEX model)
################################################################################
# Set prior hyperparameters for EXNEX
prior.B0<-c(0,100)
prior.B1<-c(0,50)
a<-rep(1.1,k)
b<-rep(1.1,k)
prior.sig.HN0<-1
prior.sig.HN1<-1 
pMix = c(0.5, 0.5)
nex.beta0 <- rep(0, k)
nex.sig.beta0 <- rep(1e2, k) 
nex.beta1 <- rep(0, k) 
nex.sig.beta1 <- rep(50, k) 

# Set cutoff.theta. Posterior probability, pCat = P(theta<cutoff.theta)
cutoff.theta = 1

# Parameters to return
parameters<- c('theta', 'pCat')

set.seed(1234)
exnex_data<-list(totaln=totaln, Module=Module, k=k, t=t,is.censored=is.censored,t.cen=t.cen,Trt=Trt, 
                prior.B1=prior.B1, 
                a=a,b=b, 
                prior.sig.HN1=prior.sig.HN1, 
                nex.beta0=nex.beta0, nex.sig.beta0=nex.sig.beta0,
                nex.beta1 = nex.beta1, nex.sig.beta1 = nex.sig.beta1, 
                pMix = pMix, 
                cutoff.theta=cutoff.theta)

parameters<- c('theta', 'pCat')
exnex_jags <- jags(data=exnex_data, 
                    parameters.to.save=parameters, 
                    model.file="EXNEX-W.txt",
                    n.chains = 4,
                    n.iter = 4000,
                    n.burnin = 2000)
