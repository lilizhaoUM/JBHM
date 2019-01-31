rm(list=ls())
require(rjags)

#----------------------------------------------------------------------------------------------
# simulate survival times from two weibull distributions for two groups (DIPG and HGG in the paper)
# Weibull distribution is parametrized by rweibull(n, a, b) with S(x) = exp(- (x/b)^a)

# input:
#   N is the total sample size
#   N1 is the sample size at stage one
#   prop is the proportion subjects in group 2
#   wd1 and wd2 are the two landmarking time point for the two groups, respectively
#   rate_accrual is the patient enrollment rate;
#   a is the shape parameter in the weibull distribution, and 
#   b1 and b2 are scale parameters for the two groups, respectively

# output: 
#   y1 is the vector of failure indicators and w1 is vector of uniform weights for subjects at stage one;
#   y is the vector of failure indicators for all subjects at the end of the study
#----------------------------------------------------------------------------------------------

Simdata<-function(N=20,N1=10,prop=0.5,rate_accrual=1/2.4,wd1=8,wd2=12,a=1,b1=b1,b2=b2){
  
  group=rbinom(N,1,prop)+1   # group indicator for DIPG and HGG
  b=c(b1,b2)
  wd=c(wd1,wd2)
  
  #--simulate survival data from weibull
  yt=sapply(1:N,function(j) rweibull(1,a,b[group[j]])) # S(x) = exp(- (x/b)^a) and rweibull(n, a, b)
  arrival_t=cumsum(rexp(N, rate_accrual))
  # the first patient start at time 0
  arrival=arrival_t-arrival_t[1]
  
  # the time at first stage (time enroll the N1+1 patient)
  T1=arrival[N1+1]
  # Kd is follow up for patients at stage 1
  Kd=T1-arrival[arrival<T1]
  WD=sapply(1:N, function(j) wd[group[j]])
  WD1=WD[arrival<T1]
  minfollow=pmin(WD1, Kd)
  
  y1=ifelse(yt[1:N1]<=minfollow,1,0)
  w1=ifelse((Kd<WD1 & yt[1:N1]>Kd), Kd/WD1,1)  # uniform weights
  
  # At the final stage (last enrolled patient + 12 months)
  y=ifelse(yt<=WD,1,0)
  list(y1=y1,w1=w1,y=y) 
}

#-----------------------------------------------------------------------
# jags files for for K strata
#---------------------------------------------------------------------------
modelmultiple="
model{
for(i in 1:N){
y[i]  ~  dbin(p[i]*w[i], 1)
logit(p[i]) <- alpha[g[i]] + beta[g[i]]
}

for(k in 1:K){
beta[k] ~ dnorm(mu.beta,tau)
}

mu.beta ~ dnorm(0,prior.beta)

# uniform prior for standard deviation
tau <- pow(sigma, -2)
sigma ~ dunif(0, prior.sigma)

# Or the inverse gamma prior for precision
#tau ~ dgamma(1, prior.sigma) 

for(k in 1:K){ 
pe[k]=1/(1+exp(alpha[k]+beta[k]))  # success prob
}
}"

#----------------------------------------------------------------
# jags files with a single stratum
#------------------------------------------------------------------
modelsingle="
model{
for(i in 1:N){
y[i]  ~  dbin(p[i]*w[i], 1)
logit(p[i]) <- alpha + beta
}
beta ~ dnorm(0,prior.beta)
pe=1/(1+exp(alpha+beta))
}"


#----------------------------------------------------------------------------------------------
# Function for running one trial, given data at first stage (ys1, ws1,gs1) 
#     and data at second stage (ys2, gs2)
# H0: null success rates;  H1: alternative success rates
# pL: the threshold for stopping
# ns:  the number of MCMC interations in jags
# beta ~ N(0,prior.beta) & mu.beta ~ N(0,prior.beta)
# sigma ~ U(0, prior.sigma)
#----------------------------------------------------------------------------------------------

JBHM <- function(ys1,ws1,gs1,ys2,ws2,gs2,H0=H0,H1=H1,pL=0.1,prior.beta=prior.beta,prior.sigma=prior.sigma,ns=niter){  
  
  ngrp=length(unique(gs1)) # number of strata
  # uninteresting (NULL) rate are transferred to logodds
  alpha=H0
  for(k in 1:ngrp){
    alpha[k]=log((1-H0[k])/H0[k])  
  }
  
  #------------------------------------------------------------------
  #  stage one
  #------------------------------------------------------------------
  data=list(y=ys1,g=gs1,N=length(ys1),K=ngrp,w=ws1,alpha=alpha,prior.sigma=prior.sigma,prior.beta=prior.beta)
  init=list(beta=c(0,0),mu.beta=0,sigma=1)
  model=jags.model(textConnection(modelmultiple),data=data,inits=init, quiet=TRUE)
  update(model,n.iter=1000,progress.bar="none")
  output=coda.samples(model=model,variable.names=c("pe"), n.iter=ns,thin=1)
  Mat = as.matrix(output)
  prob_s1=rep(0,ngrp)
  for(k in 1:ngrp){
    prob_s1[k]=sum(ifelse(Mat[,k]>H1[k],1,0))/ns  #stop for stratum 1
  }
  
  stop1=ifelse(prob_s1[1]< pL,1,0) # 1 indicates stopping for stratum 1
  stop2=ifelse(prob_s1[2]< pL,1,0) # 1 indicates stopping for stratum 2
  
  #------------------------------------------------------------------
  # final stage
  #------------------------------------------------------------------
  prob_s2=rep(0,ngrp)
  
  if(stop1==0 & stop2==0){
    
    init=list(beta=c(0,0),mu.beta=0,sigma=1)
    data=list(y=ys2,g=gs2,N=length(ys2),K=ngrp,w=ws2,alpha=alpha,prior.sigma=prior.sigma,prior.beta=prior.beta)
    model=jags.model(textConnection(modelmultiple),data=data,inits=init, quiet=TRUE)
    update(model,n.iter=1000,progress.bar="none")
    output=coda.samples(model=model,variable.names=c("pe"), n.iter=ns,thin=1)
    Mat = as.matrix(output)
    for(k in 1:ngrp){
      prob_s2[k]=sum(ifelse(Mat[,k]>H0[k],1,0))/ns  
    }
    
  }else if(stop1==1 & stop2==0){ # only the second group with the survival outcome in stage 2
    
    init=list(beta=0)
    data=list(y=ys2[gs2==2],N=sum(gs2==2),w=ws2,alpha=alpha[2],prior.beta=prior.beta)
    model=jags.model(textConnection(modelsingle),data=data,inits=init, quiet=TRUE)
    update(model,n.iter=1000,progress.bar="none")
    output=coda.samples(model=model,variable.names=c("pe"),n.iter=ns,thin=1)
    Mat = as.matrix(output)
    prob_s2[2]=sum(ifelse(Mat>H0[2],1,0))/ns
    
  }else if(stop1==0 & stop2==1){ # only the first group with the survival outcome in stage 2
    
    init=list(beta=0)
    data=list(y=ys2[gs2==1],N=sum(gs2==1),w=ws2,alpha=alpha[1],prior.beta=prior.beta)
    model=jags.model(textConnection(modelsingle),data=data,inits=init, quiet=TRUE)
    update(model,n.iter=1000,progress.bar="none")
    output=coda.samples(model=model,variable.names=c("pe"), n.iter=ns,thin=1)
    Mat = as.matrix(output)
    prob_s2[1]=sum(ifelse(Mat>H0[1],1,0))/ns
    
  }
  
  list(prob_s1=prob_s1,prob_s2=prob_s2,stop=c(stop1,stop2))
}

######################################
# Run simulations
#########################################
n=c(20,20)
Ng1=10; Ng2=10  # sample size at stage one
N=sum(n)
H0=c(0.1,0.5) # null 
H1=c(0.3,0.75) # alternative

# for PFS outcome
wd1=8   # PFS at 8 months
wd2=12  # PFS at 1 year
a.shape=1   # exponential distribution
prior.sigma=1  
prior.beta=0.01

BB=2000
  
  #-----------------------------------------------------------
  # obtain type I error simulated under H0
  #-----------------------------------------------------------
Prob_s2=Stop=matrix(0,BB,2)
  for (ntrial in 1:BB){
    
    set.seed(ntrial*10+10+454*ntrial)
    # simulate all subjects for the stratum with a survival outcome (this stratum has two subgroups)
    simd=Simdata(N=n[2],N1=Ng2,prop=0.5,rate_accrual=1/2.4,wd1=8,wd2=12,
                 a=a.shape,b1=wd1/((-log(H0[2]))^(1/a.shape)),b2=wd2/((-log(H0[2]))^(1/a.shape)))
    
    # simulate all subjects for the stratum with a binary outcome 
    yg1=1-rbinom(n[1],1,H0[1]) 
    
    # data at stage one: combine two strata
    ys1=c(yg1[1:Ng1],simd$y1) # failure indicators for both strata in stage one
    ws1=c(rep(1,Ng1),simd$w1) # weights for both strata in stage one (for the group with a binary outcome, all weights are one)
    gs1=rep(c(1,2),c(Ng1,Ng2))
    
    # data at stage two (no weight as all reached the study endpoint)
    ys2=c(yg1,simd$y) # failure indicators for both strata in stage one
    ws2=rep(1,N) # all weights are one as all subjects have complete information
    gs2=rep(c(1,2),n)
    
    res=JBHM(ys1,ws1,gs1,ys2,ws2,gs2,H0=H0,H1=H1,pL=0.1,
             prior.beta=prior.beta,prior.sigma=prior.sigma,ns=2000)
    
    Stop[ntrial,]=res$stop
    Prob_s2[ntrial,]=res$prob_s2
  } 

  Reject=colMeans(ifelse(Prob_s2>0.9,1,0))
  PET=colMeans(Stop)
  Es=PET*c(Ng1,Ng2)+(1-PET)*n
  res_H0=cbind.data.frame(type1=Reject,Es_H0=Es,PET_H0=PET)

  
  #-----------------------------------------------------------
  # obtain power simulated under H1
  #-----------------------------------------------------------
  Stop=Prob_s2=matrix(0,BB,2)
  for (ntrial in 1:BB){
    
    set.seed(ntrial*10+10+454*ntrial)
    # simulate all subjects for the stratum with a survival outcome (this stratum has two subgroups)
    simd=Simdata(N=n[2],N1=Ng2,prop=0.5,rate_accrual=1/2.4,wd1=8,wd2=12,
                 a=a.shape,b1=wd1/((-log(H1[2]))^(1/a.shape)),b2=wd2/((-log(H1[2]))^(1/a.shape)))
    
    # simulate all subjects for the stratum with a binary outcome 
    yg1=1-rbinom(n[1],1,H1[1]) 
    
    # data at stage one: combine two strata
    ys1=c(yg1[1:Ng1],simd$y1) # failure indicators for both strata in stage one
    ws1=c(rep(1,Ng1),simd$w1) # weights for both strata in stage one (for the group with a binary outcome, all weights are one)
    gs1=rep(c(1,2),c(Ng1,Ng2))
    
    # data at stage two (no weight as all reached the study endpoint)
    ys2=c(yg1,simd$y) # failure indicators for both strata in stage one
    ws2=rep(1,N) # all weights are one as all subjects have complete information
    gs2=rep(c(1,2),n)
    
    res=JBHM(ys1,ws1,gs1,ys2,ws2,gs2,H0=H0,H1=H1,pL=0.1,
             prior.beta=prior.beta,prior.sigma=prior.sigma,ns=2000)
    
    Stop[ntrial,]=res$stop
    Prob_s2[ntrial,]=res$prob_s2
  } 
  
  Reject=colMeans(ifelse(Prob_s2>0.9,1,0))
  PET=colMeans(Stop)
  Es=PET*c(Ng1,Ng2)+(1-PET)*n
  res_H1=cbind.data.frame(power=Reject,Es_H1=Es,PET_H1=PET)
  
  res_H0;res_H1;
  
