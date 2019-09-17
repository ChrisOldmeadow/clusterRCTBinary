
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





sim_cRCT <-function(m,p1,p2,k1,k2,rho,alpha,nreps){
  
 
  
  # specify the parameters the simulate function:   
      # "theta" - in the case of single variance terms, that's just the standard deviation (on the logit scale) of each random effect: e.g. theta=c(1,0.5) (among-individual variation is 4x among-observation variance). 
      # "beta" is the fixed-effects parameters, in this case (intercept,treat) – also all on the logit scale.
 
  
  expdat <- expand.grid(iid = c(1:m), site = c(1:k), ttt = c(0,1))
  
  expdat$site <- factor(ifelse(expdat$ttt == 0,expdat$site, expdat$site + k))
  
  
  set.seed(101)
  
  
  o2<-p2/(1-p2)
  o1<- p1/(1-p1)
  
  beta <- c(log(o2), log(o1/o2))
  names(beta)<-c("(Intercept)","ttt")
  
  
  
  sigma2 <-(pi ^ 2) / 3
  theta <-sqrt((rho*sigma2)/(1-rho))
  
  
  names(theta)<-c("site.(Intercept)")
  
  
  ss <- simulate(~ttt  + (1 | site), nsim = nreps, family = binomial, 
                 newdata = expdat, newparams = list(theta = theta,   beta = beta))
  
  
  expdat$resp <- ss[, 1]
  fit1 <- lme4::glmer(resp ~ ttt + (1 | site), family = binomial, data=expdat, control=glmerControl(optimizer="bobyqa", check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
  
  fitAll <- ldply(seq(nreps), function(i) fitsim(i))
  
  pow <- with(fitAll, mean( p < alpha, na.rm = TRUE))
  return(pow)

}


powersims <-function(p1,p2,
                     k1,k2,
                     rho,
                     alpha, 
                     nreps, 
                     start, 
                     end, 
                     by){
  
  out <- data.frame("Sample Size" = seq(start,end,by),  "Power" = sapply(
    seq(start, end, by), 
    sim_cRCT, 
    p1,p2,k1,k2,rho,alpha,nreps
  ))
  return(out)
}
