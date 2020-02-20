#################################################
#### kernel with h and delta  ###########
################################################

llOA1.1CI<-function(k0,r0,alpha,h, delta){ 
  
  datesc$kernel<-(k0/(1+((datesc$distance/r0)^alpha)))+delta
  datinf$kernel<-(k0/(1+((datinf$distance/r0)^alpha)))+delta
  datend$kernel<-(k0/(1+((datend$distance/r0)^alpha)))+delta
  datreesc$kernel<-(k0/(1+((datreesc$distance/r0)^alpha)))+delta
  
  lambdesc<-datesc[, sum(kernel),by=list(D,IDsus)]+h
  lambdinf<-datinf[, sum(kernel),by=list(D,IDsus)]+h
  lambdend<-datend[, sum(kernel),by=list(D,IDsus)]+h
  lambdreesc<-datreesc[, sum(kernel),by=list(D,IDsus)]+h
  
  
  logPesc<-sum(lambdesc[,V1]) #calculate log P escape for each susceptible farms from t=1 until t max-1
  logPinf<-sum(log(1-(exp(-1*lambdinf[,V1])))) 
  logPend<-sum(lambdend[,V1]) #calculate log P escape to the end for each susceptible farms from t=1 until t max
  logPreesc<-sum(lambdreesc[,V1]) #calculate log P escape for recover farm
  
  return ((-1)*(logPinf-logPesc-logPend-logPreesc)) # use- instead of + in formula because minimization negative loglikelihood = maximize likelihood but it is simpler
}



fit.OA1.1CI <-mle2(llOA1.1CI, start = list(k0 = 0.004268534 ,r0 = 0.2926315,alpha= 1.938818,h =  0.0017, delta = 0.002), lower = c(k0 = 0, r0 = 0, alpha = 0, h = 0, delta = 0) )
Kfit.OA1.1CI<-profile(fit.OA1.1CI, alpha = 0.05)
plot(Kfit.OA1.1CI)
confint(Kfit.OA1.1CI, level = 0.95)
