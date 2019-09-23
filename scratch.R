library(plyr)
library(lme4)
library(MASS) ## for glmmPQL


m<-10
k1<-10
k2<-8
p1<-0.18
p2<-0.03
rho<-0.01
alpha<-.05
nreps<-200

set.seed(101)

int <- expand.grid(iid = c(1:m), site = c(1:k1), ttt = 0)
ctr <- expand.grid(iid = c(1:m), site = c((k1+1):(k1+k2)), ttt = 1)
    
expdat <- rbind(int,ctr)
    
    
#expdat$site <- factor(ifelse(expdat$ttt == 0,expdat$site, expdat$site + k))
    
 
o2<-p2/(1-p2)
o1<- p1/(1-p1)

beta <- c(log(o2), log(o1/o2))
names(beta)<-c("(Intercept)","ttt")


rho <- 0.01
sigma2 <-(pi ^ 2) / 3
theta <-sqrt((rho*sigma2)/(1-rho))


names(theta)<-c("site.(Intercept)")


ss <- simulate(~ttt  + (1 | site), nsim = nreps, family = binomial, 
               newdata = expdat, newparams = list(theta = theta,   beta = beta))


expdat$resp <- ss[, 1]
fit1 <- MASS::glmmPQL(resp ~ ttt , random = ~ 1 | site, family = binomial, data=expdat)

fitsim <- function(i) {
  expdat$resp <- ss[, i]
  return <- tryCatch({
    res<-coef(summary(MASS::glmmPQL(resp ~ ttt , random = ~ 1 | site, family = binomial, data=expdat)))["ttt", ]
    names(res)<-c("est","se","df","t","p")
    return(res)
  },
  warning =function(e) {
    #message(e)  # print error message
    return(c(est=NA, se=NA, df=NA , t=NA, p=NA))
  },
  
  error=function(e) {
    #message(e)  # print error message
    return(c(est=NA, se=NA, df=NA , t=NA, p=NA))
  })
  
  return(return)
}
   
fitAll <- ldply(seq(nreps), function(i) fitsim(i))
   
pow<-with(fitAll, mean( p < 0.05, na.rm = TRUE))
   
pow



pow_simclusterRCT(m=5,k1=8,k2=10,p1=.11,p2=.03,rho=.01,alpha=.05,nreps = 500)



pow_simclusterRCT(m=15,k1=8,k2=10,p1=.13,p2=.03,rho=.01,alpha=.05,nreps = 500)



pow_simclusterRCT(m=5,k1=8,k2=10,p1=.20,p2=.03,rho=.01,alpha=.05,nreps = 500)

source("utils.R")

sapply(
    seq(from=.11, to=.21, by=0.02), 
    sim_cRCT, 
    p2=.03,k1=8,k2=10,rho=.01,alpha=.05,nreps=100,method=2,m=15)


sim_cRCT(p1=0.11,p2=.03,k1=8,k2=10,rho=.01,alpha=.05,nreps=100,method=2,m=15)
sim_cRCT(0.15,p2=.03,k1=8,k2=10,rho=.01,alpha=.05,nreps=100,method=2,m=15)

