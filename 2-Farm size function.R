#########################################################################
######### Estimate farm size function parameters ########################

############ C1 ###############################

llC1<-function(k0,r0,alpha,Cm) { 
  
  
  datesc$kernel<-(k0/(1+((datesc$distance/r0)^alpha)))*((datesc$NIav*datesc$NSav)^Cm)
  datinf$kernel<-(k0/(1+((datinf$distance/r0)^alpha)))*((datinf$NIav*datinf$NSav)^Cm)
  datend$kernel<-(k0/(1+((datend$distance/r0)^alpha)))*((datend$NIav*datend$NSav)^Cm)
  datreesc$kernel<-(k0/(1+((datreesc$distance/r0)^alpha)))*((datreesc$NIav*datreesc$NSav)^Cm)
  
  lambdesc<-datesc[, sum(kernel),by=list(D,IDsus)]
  lambdinf<-datinf[, sum(kernel),by=list(D,IDsus)]
  lambdend<-datend[, sum(kernel),by=list(D,IDsus)]
  lambdreesc<-datreesc[, sum(kernel),by=list(D,IDsus)]
  
  logPesc<-sum(lambdesc[,V1]) #calculate log P escape for each susceptible farms from t=1 until t max-1
  logPinf<-sum(log(1-(exp(-1*lambdinf[,V1])))) 
  logPend<-sum(lambdend[,V1]) #calculate log P escape to the end for each susceptible farms from t=1 until t max
  logPreesc<-sum(lambdreesc[,V1]) #calculate log P escape for recover farm
  return ((-1)*(logPinf-logPesc-logPend-logPreesc)) # use- instead of + in formula because minimization negative loglikelihood = maximize likelihood but it is simpler
}



fit.LPllC1 <-mle2( llC1,start = list(k0=0.0061  , r0 = 0.18, alpha =	1.57, Cm = 0.48), skip.hessian=FALSE,  optimizer =  "optimx")
Kfit.LPllC1 <-profile(fit.LPllC1 , alpha = 0.05)
plot(Kfit.LPllC1)
confint(Kfit.LPllC1, level = 0.95)# to get profile CI


############ C2 ###############################

llC2<-function( k0, r0, alpha, Cs,Ci) {
  
  datesc$kernel<-(k0/(1+((datesc$distance/r0)^alpha)))*((datesc$NSav)^Cs)*((datesc$NIav)^Ci)
  datinf$kernel<-(k0/(1+((datinf$distance/r0)^alpha)))*((datinf$NSav)^Cs)*((datinf$NIav)^Ci)
  datend$kernel<-(k0/(1+((datend$distance/r0)^alpha)))*((datend$NSav)^Cs)*((datend$NIav)^Ci)
  datreesc$kernel<-(k0/(1+((datreesc$distance/r0)^alpha)))*((datreesc$NSav)^Cs)*((datreesc$NIav)^Ci)
  
  lambdesc<-datesc[, sum(kernel),by=list(D,IDsus)]
  lambdinf<-datinf[, sum(kernel),by=list(D,IDsus)]
  lambdend<-datend[, sum(kernel),by=list(D,IDsus)]
  lambdreesc<-datreesc[, sum(kernel),by=list(D,IDsus)]
  
  logPesc<-sum(lambdesc[,V1]) #calculate log P escape for each susceptible farms from t=1 until t max-1
  logPinf<-sum(log(1-(exp(-1*lambdinf[,V1])))) 
  logPend<-sum(lambdend[,V1]) #calculate log P escape to the end for each susceptible farms from t=1 until t max
  logPreesc<-sum(lambdreesc[,V1]) #calculate log P escape for recover farm
  return ((-1)*(logPinf-logPesc-logPend-logPreesc)) # use- instead of + in formula because minimization negative loglikelihood = maximize likelihood but it is simpler
}

fit.LPllC2<-mle2( llC2,start = list(k0=0.0059  , r0 = 0.17, alpha =	1.55, Cs = 0.45, Ci = 0.76), skip.hessian=FALSE, method = "nlminb" , optimizer="optimx",lower = c(k0=0  , r0 = 0, alpha =	0, Cs = 0.1, Ci = 0.1), upper = c(k0=0.1  , r0 = 1, alpha =	3, Cs = 1, Ci = 3))
Kfit.LPllC2<-profile(Kfit.LPllC2, alpha = 0.05)
plot(Kfit.LPllC2)
confint(Kfit.LPllC2, level = 0.95)

############# C3 ###############################

llC3<-function( k0 , r0 ,alpha,Am ) { 
  
  
  datesc$kernel<-(k0/(1+((datesc$distance/r0)^alpha)))*(1+(Am*((datesc$NSav)-1)*((datesc$NIav)-1)))
  datinf$kernel<-(k0/(1+((datinf$distance/r0)^alpha)))*(1+(Am*((datinf$NSav)-1)*((datinf$NIav)-1)))
  datend$kernel<-(k0/(1+((datend$distance/r0)^alpha)))*(1+(Am*((datend$NSav)-1)*((datend$NIav)-1)))
  datreesc$kernel<-(k0/(1+((datreesc$distance/r0)^alpha)))*(1+(Am*((datreesc$NSav)-1)*((datreesc$NIav)-1)))
  
  lambdesc<-datesc[, sum(kernel),by=list(D,IDsus)]
  lambdinf<-datinf[, sum(kernel),by=list(D,IDsus)]
  lambdend<-datend[, sum(kernel),by=list(D,IDsus)]
  lambdreesc<-datreesc[, sum(kernel),by=list(D,IDsus)]
  
  logPesc<-sum(lambdesc[,V1]) #calculate log P escape for each susceptible farms from t=1 until t max-1
  logPinf<-sum(log(1-(exp(-1*lambdinf[,V1])))) 
  logPend<-sum(lambdend[,V1]) #calculate log P escape to the end for each susceptible farms from t=1 until t max
  logPreesc<-sum(lambdreesc[,V1]) #calculate log P escape for recover farm
  return ((-1)*(logPinf-logPesc-logPend-logPreesc)) # use- instead of + in formula because minimization negative loglikelihood = maximize likelihood but it is simpler
}


fit.LPllC3 <-mle2( llC3, start = list(k0=0.0059  , r0 = 0.17, alpha =	1.55, Am = 0.5), skip.hessian=FALSE, method = "nlminb" , optimizer="optimx",lower = c(k0=0  , r0 = 0, alpha =	0, Am=0), upper = c(k0=0.1  , r0 = 1, alpha =	3, Am=2))
Kfit.LPllC3 <-profile(Kfit.LPllC3, alpha = 0.05)
plot(Kfit.LPllC3)
confint(Kfit.LPllC3, level = 0.95)

############# C4 ###############################

llc4<-function(k0 , r0 ,alpha,As,Ai ){
  datesc$kernel<-(k0/(1+((datesc$distance/r0)^alpha)))*(1+(As*(datesc$NSav-1)))*(1+(Ai*(datesc$NIav-1)))
  datinf$kernel<-(k0/(1+((datinf$distance/r0)^alpha)))*(1+(As*(datinf$NSav-1)))*(1+(Ai*(datinf$NIav-1)))
  datend$kernel<-(k0/(1+((datend$distance/r0)^alpha)))*(1+(As*(datend$NSav-1)))*(1+(Ai*(datend$NIav-1)))
  datreesc$kernel<-(k0/(1+((datreesc$distance/r0)^alpha)))*(1+(As*(datreesc$NSav-1)))*(1+(Ai*(datreesc$NIav-1)))
  
  lambdesc<-datesc[, sum(kernel),by=list(D,IDsus)]
  lambdinf<-datinf[, sum(kernel),by=list(D,IDsus)]
  lambdend<-datend[, sum(kernel),by=list(D,IDsus)]
  lambdreesc<-datreesc[, sum(kernel),by=list(D,IDsus)]
  
  logPesc<-sum(lambdesc[,V1]) #calculate log P escape for each susceptible farms from t=1 until t max-1
  logPinf<-sum(log(1-(exp(-1*lambdinf[,V1])))) 
  logPend<-sum(lambdend[,V1]) #calculate log P escape to the end for each susceptible farms from t=1 until t max
  logPreesc<-sum(lambdreesc[,V1]) #calculate log P escape for recover farm
  return ((-1)*(logPinf-logPesc-logPend-logPreesc)) # use- instead of + in formula because minimization negative loglikelihood = maximize likelihood but it is simpler
}


fit.LPllc4 <-mle2( llc4, start = list(k0=0.0059  , r0 = 0.17, alpha =	1.55, As = 0.5, Ai = 0.7), skip.hessian=FALSE, method = "nlminb" , optimizer="optimx",lower = c(k0=0  , r0 = 0, alpha =	0, As=0, Ai = 0), upper = c(k0=0.1  , r0 = 1, alpha =	3, As=3, Ai=3))
Kfit.LPllc4<-profile(fit.LPllc4, alpha = 0.05)
plot(Kfit.LPllc4)
confint(Kfit.LPllc4, level = 0.95)

############# C5 ###############################
llc5<-function(k0,r0,alpha,d){
  
  datesc$kernel<-(k0/(1+((datesc$distance/r0)^alpha)))*((1-exp(-1*(datesc$NSav/d)))*(1-exp(-1*(datesc$NIav/d))))
  datinf$kernel<-(k0/(1+((datinf$distance/r0)^alpha)))*((1-exp(-1*(datinf$NSav/d)))*(1-exp(-1*(datinf$NIav/d))))
  datend$kernel<-(k0/(1+((datend$distance/r0)^alpha)))*((1-exp(-1*(datend$NSav/d)))*(1-exp(-1*(datend$NIav/d))))
  datreesc$kernel<-(k0/(1+((datreesc$distance/r0)^alpha)))*((1-exp(-1*(datreesc$NSav/d)))*(1-exp(-1*(datreesc$NIav/d))))
  
  lambdesc<-datesc[, sum(kernel),by=list(D,IDsus)]
  lambdinf<-datinf[, sum(kernel),by=list(D,IDsus)]
  lambdend<-datend[, sum(kernel),by=list(D,IDsus)]
  lambdreesc<-datreesc[, sum(kernel),by=list(D,IDsus)]
  
  logPesc<-sum(lambdesc[,V1]) #calculate log P escape for each susceptible farms from t=1 until t max-1
  logPinf<-sum(log(1-(exp(-1*lambdinf[,V1])))) 
  logPend<-sum(lambdend[,V1]) #calculate log P escape to the end for each susceptible farms from t=1 until t max
  logPreesc<-sum(lambdreesc[,V1]) #calculate log P escape for recover farm
  return ((-1)*(logPinf-logPesc-logPend-logPreesc)) # use- instead of + in formula because minimization negative loglikelihood = maximize likelihood but it is simpler
}


fit.LPllc5 <-mle2( llc5, start = list(k0=0.0059  , r0 = 0.17, alpha =	1.55, d =0.67), skip.hessian=FALSE, method = "nlminb" , optimizer="optimx",lower = c(k0=0  , r0 = 0, alpha =	0, d=0.1), upper = c(k0=0.1  , r0 = 1, alpha =	3, d = 1.5))
Kfit.LPllc5 <-profile(fit.LPllc5 , alpha = 0.05)
plot(Kfit.LPllc5)
confint(Kfit.LPllc5, level = 0.95)

############# C6 ###############################

llc6<-function(k0,r0,alpha,ds,di){ 
  
  datesc$kernel<-(k0/(1+((datesc$distance/r0)^alpha)))*((1-exp(-1*(datesc$NSav/ds)))*(1-exp(-1*(datesc$NIav/di))))
  datinf$kernel<-(k0/(1+((datinf$distance/r0)^alpha)))*((1-exp(-1*(datinf$NSav/ds)))*(1-exp(-1*(datinf$NIav/di))))
  datend$kernel<-(k0/(1+((datend$distance/r0)^alpha)))*((1-exp(-1*(datend$NSav/ds)))*(1-exp(-1*(datend$NIav/di))))
  datreesc$kernel<-(k0/(1+((datreesc$distance/r0)^alpha)))*((1-exp(-1*(datreesc$NSav/ds)))*(1-exp(-1*(datreesc$NIav/di))))
  
  lambdesc<-datesc[, sum(kernel),by=list(D,IDsus)]
  lambdinf<-datinf[, sum(kernel),by=list(D,IDsus)]
  lambdend<-datend[, sum(kernel),by=list(D,IDsus)]
  lambdreesc<-datreesc[, sum(kernel),by=list(D,IDsus)]
  
  logPesc<-sum(lambdesc[,V1]) #calculate log P escape for each susceptible farms from t=1 until t max-1
  logPinf<-sum(log(1-(exp(-1*lambdinf[,V1])))) 
  logPend<-sum(lambdend[,V1]) #calculate log P escape to the end for each susceptible farms from t=1 until t max
  logPreesc<-sum(lambdreesc[,V1]) #calculate log P escape for recover farm
  return ((-1)*(logPinf-logPesc-logPend-logPreesc)) # use- instead of + in formula because minimization negative loglikelihood = maximize likelihood but it is simpler
}


fit.LPllc6<-mle2( llc6, start = list(k0=0.0005  , r0 = 0.17, alpha =	1.55, ds =0.67,di =0.7), skip.hessian=FALSE, method = "nlminb" , optimizer="optimx",lower = c(k0=0  , r0 = 0, alpha =	0, d=0.1, di =0.1), upper = c(k0=0.03  , r0 = 1, alpha =	3, ds=1.17, di =1.17))
Kfit.LPllc6<-profile(fit.LPllc6, alpha = 0.05)
plot(Kfit.LPllc6)
confint(Kfit.LPllc6, level = 0.95)