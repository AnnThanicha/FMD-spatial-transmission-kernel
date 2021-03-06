##################################################################################
################ prepare distance data in long format ############################
library(reshape2)

distLPEClong<-melt(data= LPECdist) ## reshape matrix to long format data
distLPEClong<- setNames(distLPEClong, c("IDsus","IDinf","distance"))
summary(distLPEClong)

######## Find which distribution are fitted with the infectious duration from outbreak data ########

library(fitdistrplus)
FMDdata$duration<-FMDdata$stopinfectiousdate-FMDdata$startinfectiousdate # calculate infectious period
duration<-FMDdata$duration[!is.na(FMDdata$duration)]+1#remove NA value 
hist(duration)
summary(duration)
descdist(duration, discrete = TRUE, boot=1000)
gamma_duration <- fitdist(duration, "gamma") # You can test by trying other distribution and check which one is the best fit.
plot(gamma_duration) # I decide to use gamma distribution for stochatically generating infectious duration


######### Calculate the Rh ###########################
Rh<-function(dat,k0,r0,alpha,duration){
  dat$kernel<-(k0/(1+((dat$distance/r0)^alpha)))
  dat$pinf<-1-(exp(-1*dat$kernel*duration))
  dat[which(dat[,"IDsus"] == dat[,"IDinf"]), "pinf"]<-0  # probability of infection in the same farm = 0
  dat1<-aggregate(pinf~IDinf, sum, data=dat) # aggregate probability of infection generated by each infectious farm
  return( dat1)
}
set.seed(01)
baselineRh<-Rh(dat = distLPEClong, k0=3.372114e-03,r0= 4.316170e-01, alpha = 2.803463e+00, duration = rgamma(nrow(distLPEClong),shape =  1.86604707,rate = 0.09852991 ))
#infectious duration is randomized from gamma distribution based on infectious duration of outbreak data

write.table(baselineRh,  "D:/Phd thesis/R code/R-lumpayaklang2/Rh.txt", sep="\t")# export reproduction number
#I will use this csv to make a risk map in Qgis.
# check this link for how to make a map in Qgis : https://www.polarmicrobes.org/tutorial-on-qgis-how-to-make-a-map/
#Another options using sf package to making a map.

