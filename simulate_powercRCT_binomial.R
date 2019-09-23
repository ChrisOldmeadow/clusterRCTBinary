library(plyr)
library(lme4)
library(MASS) ## for glmmPQL





pow_simclusterRCT<-function(m,k1,k2,p1,p2,rho,alpha,nreps){

    int <- expand.grid(iid = c(1:m), site = c(1:k1), ttt = 0)
    ctr <- expand.grid(iid = c(1:m), site = c((k1+1):(k1+k2)), ttt = 1)
    
    expdat <- rbind(int,ctr)
    
    
    #expdat$site <- factor(ifelse(expdat$ttt == 0,expdat$site, expdat$site + k))
    
    
    set.seed(101)
    nsim <- 200
    
    o2<-p2/(1-p2)
    o1<- p1/(1-p1)
    
    beta <- c(log(o2), log(o1/o2))
    names(beta)<-c("(Intercept)","ttt")
    
    
    rho <- 0.01
    sigma2 <-(pi ^ 2) / 3
    theta <-sqrt((rho*sigma2)/(1-rho))
    
    
    names(theta)<-c("site.(Intercept)")
    
    
    ss <- simulate(~ttt  + (1 | site), nsim = nsim, family = binomial, 
                   newdata = expdat, newparams = list(theta = theta,   beta = beta))
    
    
    expdat$resp <- ss[, 1]
    fit1 <- lme4::glmer(resp ~ ttt + (1 | site), family = binomial, data=expdat, control=glmerControl(optimizer="bobyqa", check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    
   # fit1 <- MASS::glmmPQL(resp ~ ttt , random = 1 | site, family = binomial, data=expdat)
    
    fitsim <- function(i) {
      
      return <- tryCatch({
        res<-coef(summary(lme4::refit(fit1, ss[[i]])))["ttt", ]
        names(res)<-c("est","se","z","p")
        return(res)
      },
      warning =function(e) {
        #message(e)  # print error message
        return(c(est=NA, se=NA, z=NA, p=NA))
      },
      
      error=function(e) {
        #message(e)  # print error message
        return(c(est=NA, se=NA, z=NA, p=NA))
      })
      
      return(return)
    }
    
    fitAll <- ldply(seq(nsim), function(i) fitsim(i))
    
    pow<-with(fitAll, mean( p < 0.05, na.rm = TRUE))
    
    
    return(pow)

}

pow_simclusterRCT(m=5,k1=8,k2=10,p1=.11,p2=.03,rho=.01,alpha=.05,nreps = 500)



pow_simclusterRCT(m=15,k1=8,k2=10,p1=.13,p2=.03,rho=.01,alpha=.05,nreps = 500)



pow_simclusterRCT(m=5,k1=8,k2=10,p1=.20,p2=.03,rho=.01,alpha=.05,nreps = 500)


