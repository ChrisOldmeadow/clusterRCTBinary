library(plyr)

m<- 100 # average particpants per site
k <- 10 #sites per arm
p1<-0.1
p2<-0.18
expdat <- expand.grid(iid = c(1:m), site = c(1:k), ttt = c(0,1))

expdat$site <- factor(ifelse(expdat$ttt == 0,expdat$site, expdat$site + k))


set.seed(101)
nsim <- 200


beta <- c(qlogis(p2), qlogis(p2-p1))
names(beta)<-c("(Intercept)","ttt")


rho <- 0.05
sigma2 <-(pi ^ 2) / 3
theta <-sqrt((rho*sigma2)/(1-rho))


names(theta)<-c("site.(Intercept)")


ss <- simulate(~ttt  + (1 | site), nsim = nsim, family = binomial, 
               newdata = expdat, newparams = list(theta = theta,   beta = beta))


expdat$resp <- ss[, 1]
fit1 <- lme4::glmer(resp ~ ttt + (1 | site), family = binomial, data=expdat, control=glmerControl(optimizer="bobyqa", check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))

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

with(fitAll, mean( p < 0.05, na.rm = TRUE))





