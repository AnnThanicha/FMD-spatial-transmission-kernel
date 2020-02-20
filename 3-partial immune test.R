############### Randomly assign immune status to susceotible farms ####################
#########################################################################

#### function for randomly assign immune status to susceotible farms ####

random_immmune<- function (stepimmune, immunepercent){ 
  for(i in 1:round(335/stepimmune)){
    set.seed(111)
    sampleimmune<-FMDdata$ID[(sample( which(is.na(FMDdata$stopinfectiousdate)), round(immunepercent*length(which(is.na(FMDdata$stopinfectiousdate))))))]
    FMDdata$stopinfectiousdate[FMDdata$ID %in%  sampleimmune ]<-stepimmune*i
  }
  FMDdata$waningdate<-FMDdata$stopinfectiousdate+225
  FMDdata$waningdate[FMDdata$waningdate>335] <-NA #waning date more than 335 is NA
  return(FMDdata)
}

##################### fucntion for kernel estimation ###########################################
ll1CI <-function(k0,r0,alpha) { 
  
  datesc$kernel<-(k0/(1+((datesc$distance/r0)^alpha))) 
  datinf$kernel<-(k0/(1+((datinf$distance/r0)^alpha)))
  datend$kernel<-(k0/(1+((datend$distance/r0)^alpha)))
  datreesc$kernel<-(k0/(1+((datreesc$distance/r0)^alpha)))
  
  lambdesc<-datesc[, sum(kernel),by=list(D,IDsus)]
  lambdinf<-datinf[, sum(kernel),by=list(D,IDsus)]
  lambdend<-datend[, sum(kernel),by=list(D,IDsus)]
  lambdreesc<-datreesc[, sum(kernel),by=list(D,IDsus)]
  
  
  logPesc<-sum(lambdesc[,V1]) #calculate log P escape for each susceptible farms from t=1 until t max-1
  logPinf<-sum(log(1-(exp(-1*lambdinf[,V1])))) 
  logPend<-sum(lambdend[,V1]) #calculate log P escape to the end for each susceptible farms from t=1 until t max
  logPreesc<-sum(lambdreesc[,V1]) #calculate log P escape for recover farm
  
  return ((-1)*(logPinf-logPesc-logPend-logPreesc)) # because mle2 will find the parameters that minimized negative loglikelihood (i.e., maximize likelihood) so the return value have to be negative loglikelihood 
}

############5% immune every 30 days##################################################
rand30day_0.05<-random_immune(stepimmune = 30, immunepercent=0.05)
summary(rand30day_0.05)

datesc_rand30day_0.05<-prep_datesc(x=rand30day_0.05,n=max(rand30day_0.05$stopinfectiousdate,na.rm=TRUE)-1, c=LPECdistkm) #x= outbreak data, c=distance matrix time max-1 because it is a group of escape until t infection
datinf_rand30day_0.05<-prep_datinf(x=rand30day_0.05,n=max(rand30day_0.05$stopinfectiousdate,na.rm=TRUE),c=LPECdistkm)
datend_rand30day_0.05<-prep_datend(x=rand30day_0.05,n=max(rand30day_0.05$stopinfectiousdate,na.rm=TRUE),c=LPECdistkm)
datreesc_rand30day_0.05<-prep_datreesc(x=rand30day_0.05,n=max(rand30day_0.05$stopinfectiousdate,na.rm=TRUE),c=LPECdistkm)

datesc<-datesc_rand30day_0.05
datinf<-datinf_rand30day_0.05
datend<-datend_rand30day_0.05
datreesc<-datreesc_rand30day_0.05


fit.rand30day_0.05_1CI <-mle2( ll1CI,  start = list(k0=0.0153237 , r0 = 1.162419, alpha =	0.7598452 ), skip.hessian=FALSE, optimizer =  "optimx")
rand30day_0.05_1CI<-profile(fit.rand30day_0.05_1CI, alpha = 0.05)
plot(rand30day_0.05_1CI)
confint(rand30day_0.05_1CI, level = 0.95)# to get profile CI


######10% immune every 30days #######################################################
rand30day_0.1<-random_immune(stepimmune = 30, immunepercent=0.1)
summary(rand30day_0.1)

datesc_rand30day_0.1<-prep_datesc(x=rand30day_0.1,n=max(rand30day_0.1$stopinfectiousdate,na.rm=TRUE)-1, c=LPECdistkm) #x= outbreak data, c=distance matrix time max-1 because it is a group of escape until t infection
datinf_rand30day_0.1<-prep_datinf(x=rand30day_0.1,n=max(rand30day_0.1$stopinfectiousdate,na.rm=TRUE),c=LPECdistkm)
datend_rand30day_0.1<-prep_datend(x=rand30day_0.1,n=max(rand30day_0.1$stopinfectiousdate,na.rm=TRUE),c=LPECdistkm)
datreesc_rand30day_0.1<-prep_datreesc(x=rand30day_0.1,n=max(rand30day_0.1$stopinfectiousdate,na.rm=TRUE),c=LPECdistkm)

datesc<-datesc_rand30day_0.1
datinf<-datinf_rand30day_0.1
datend<-datend_rand30day_0.1
datreesc<-datreesc_rand30day_0.1

fit.rand30day_0.1_1CI <-mle2(ll1CI,  start = list(k0=0.01708773, r0 = 1.121095, alpha =	0.7554839), skip.hessian=FALSE, optimizer =  "optimx")
rand30day_0.1_1CI<-profile(fit.rand30day_0.1_1CI, alpha = 0.05)
plot(rand30day_0.1_1CI)
confint(rand30day_0.1_1CI, level = 0.95)# to get profile CI



####################################################################################################
######  Assumption that immune farms did not became susceptible ##########################

ll1wCI<-function(k0,r0,alpha) { # this function is for
  
  datesc$kernel<-(k0/(1+((datesc$distance/r0)^alpha)))
  datinf$kernel<-(k0/(1+((datinf$distance/r0)^alpha)))
  datend$kernel<-(k0/(1+((datend$distance/r0)^alpha)))
  
  lambdesc<-datesc[, sum(kernel),by=list(D,IDsus)]
  lambdinf<-datinf[, sum(kernel),by=list(D,IDsus)]
  lambdend<-datend[, sum(kernel),by=list(D,IDsus)]
  
  
  logPesc<-sum(lambdesc[,V1]) #calculate log P escape for each susceptible farms from t=1 until t max-1
  logPinf<-sum(log(1-(exp(-1*lambdinf[,V1])))) 
  logPend<-sum(lambdend[,V1]) #calculate log P escape to the end for each susceptible farms from t=1 until t max
  
  return ((-1)*(logPinf-logPesc-logPend)) # use- instead of + in formula because minimization negative loglikelihood = maximize likelihood but it is simpler
}


fitll1wCI <-mle2( ll1wCI,start = list(k0=0.005376565 , r0 = 0.188094559, alpha =	1.563751854), skip.hessian=FALSE,  optimizer =  "nlminb")#trace>0 to get output
Kfitll1wCI<-profile(fitll1wCI, alpha = 0.05)
plot(Kfitll1wCI)
confint(Kfitll1wCI, level = 0.95)# to get profile CI

#######################################################################################################
########################### Plot partial immune kernels comparison  ##########################################
library (ggplot2)
C0 <-function(k0,r0,alpha, distance){ k0/(1+((distance/r0)^alpha))}
t<-seq(0,3,0.01)
kernelC0<-C0(distance=t, k0=0.005376565, r0=0.188094559, alpha = 1.563751854 )  
kernelrand30day_0.05<-C0(distance=t, k0=0.005416304, r0= 0.1916735, alpha =  1.538266) 
kernelrand30day_0.1<-C0(distance=t, k0=0.005704608 , r0=0.2018732, alpha = 1.552981)   
kernelrand30day_0.2<-C0(distance=t, k0=0.005956393, r0=0.2103479, alpha = 1.514261 )  

kernelcompare<-data.frame(t,kernelC0,kernelrand30day_0.05,kernelrand30day_0.1,kernelrand30day_0.2)

labelx<-seq(0,3,1)

plotimmune<-ggplot(data = kernelcompare, aes(x = t)) + 
  geom_line(aes(x=t,y = kernelC0, color="0%"), size=1)+
  geom_line(aes(x=t,y = kernelrand30day_0.05, color="5%"), size=1)+
  geom_line(aes(x=t,y = kernelrand30day_0.1, color="10%"), size=1)+
  geom_line(aes(x=t,y = kernelrand30day_0.2, color="20%"), size=1)+
  scale_color_manual(name = "country", values = c("0%"= "black", "5%" = "red", "10%"="green","20%"="blue"))


plotimmune<-plotimmune+xlab("distance(km)") + ylab("transmission rate (per day)") + scale_x_continuous(breaks = labelx ,limits = c(0,3),expand = c(0, 0))+ scale_y_continuous(limits = c(0,0.008),expand = c(0, 0)) 
plotimmune

