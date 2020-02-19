########## import data ### Please check data pattern in readme file ####### 
library(readxl)
FMDdata <- read_excel("FMDdata.xlsx")
FMDdist <- read_excel("FMDdist.xlsx")

summary(FMDdata)
avpop<-mean(FMDdata$pop)#calculate average farm size

######### calculate Euclidean distance ############
library(sp)
library(sf)
library(raster)

# function for calculate Euclidean distance matrix 
cal.dist<- function (lon,lat)
{crdref <- CRS('+proj=longlat +datum=WGS84') # make sure that your latitude and longitude is in WGS84 project system
i<- cbind(c(lon),c(lat))
pointmap <- SpatialPoints(i, proj4string=crdref)
dst <- pointDistance(pointmap, lonlat= TRUE)
dst <- as.dist(dst)
D<-as.matrix(dst)
return(D)
}

LPECdist<- cal.dist(lon=FMDdist$longitude, lat = FMDdist$latitude)
LPECdist<-LPECdist/1000 # to change meters to kilometers
# Now we got the matrix of Euclidean distance between farms in km unit

########## Function to preparing data ############################
#Why did I decide to prepare data first?
#Because I try to avoid using for-loop during optimization preocess, which is not computational efficiency. 
#So I prepare data into the format that I can calculate a kernel for each pair of farm on each day. 
#This way we can vectorized function and reduce calculation time.  

####### kernel for farms that escape until the time Tinfection-1 ########
prep_datesc<-function (x,n,c) {
  datesc<-data.frame()
  for (i in 1:n) {
    a<-x$ID[which((x$startinfectiousdate<=i ) & (x$stopinfectiousdate>=i))] #a = IDinf will be infectious until stopinfect date
    b<-x$ID[which(x$infectiondate>i)] #b = IDsus
    D<-i
    IDinf<-rep(a,each=length(b)) #to replicate it
    IDsus<-rep(b,length(a))
    df<-data.frame(cbind(D,IDinf,IDsus)) #bind it together
    if (length(df) == 3) {datesc<-rbind(datesc,df)} else {next}
  }
  distance<-rep(0,nrow(datesc)) # blank vector for distance
  for (i in 1:nrow(datesc)) {distance[i]<-c[datesc$IDinf[i],datesc$IDsus[i]]} #find distance between 2 farms
  datesc$distance<-distance
  
  popinf<-rep(0,nrow(datesc))
  popsus<-rep(0,nrow(datesc))
  for (i in 1:nrow(datesc)) {
    popinf[i]<-x$pop[x$ID==datesc$IDinf[i]] # find pop in each farm
    popsus[i]<-x$pop[x$ID==datesc$IDsus[i]]
  }
  datesc<-cbind(datesc,popinf,popsus)
  datesc$NIav<-datesc$popinf/avpop # to create column with the number of animal in infectious farm/average population
  datesc$NSav<-datesc$popsus/avpop # to create column with the number of animal in susceptible farm/average population
 
  return (datesc)
}


####### kernel for farms on the Tinfection ##########
prep_datinf<-function (x,n,c){
  datinf<-data.frame()
  for (i in 1:n) {
    a<-x$ID[which((x$startinfectiousdate<=i ) & (x$stopinfectiousdate>=i))] #a = IDinf will be infectious until stopinfect date
    b<-x$ID[which(x$infectiondate==i)] #b = IDsus (that become infection this day)
    D<-i 
    IDinf<-rep(a,each=length(b)) #to replicate it
    IDsus<-rep(b,length(a))
    df<-data.frame(cbind(D,IDinf,IDsus))
    if (length(df) == 3) {datinf<-rbind(datinf,df)} else {next}
  }
  distance<-rep(0,nrow(datinf)) # blank vector for distance
  for (i in 1:nrow(datinf)) {distance[i]<-c[datinf$IDinf[i],datinf$IDsus[i]]} #find distance between 2 farms
  datinf$distance<-distance
  
  popinf<-rep(0,nrow(datinf))
  popsus<-rep(0,nrow(datinf))
  for (i in 1:nrow(datinf)) {
    popinf[i]<-x$pop[x$ID==datinf$IDinf[i]] # find pop in each farm
    popsus[i]<-x$pop[x$ID==datinf$IDsus[i]]
  }
  datinf<-cbind(datinf,popinf,popsus)
  datinf$NIav<-datinf$popinf/avpop
  datinf$NSav<-datinf$popsus/avpop
 
  return(datinf)
}


######### kerkernel for farms that escape until the end of outbreak ########
prep_datend<-function (x,n,c){
  datend<-data.frame()
  for (i in 1:n) {
    a<-x$ID[which((x$startinfectiousdate<=i ) & (x$stopinfectiousdate>=i))] #a = IDinf will be infectious until stopinfect date
    b<-x$ID[(is.na(x$infectiondate))&(is.na(x$waningdate))] #ID sus which have never be infected or recovered from the start till the end 
    D<-i
    IDinf<-rep(a,each=length(b)) #to replicate it
    IDsus<-rep(b,length(a))
    df<-data.frame(cbind(D,IDinf,IDsus)) #bind it together
    if (length(df) == 3) {datend<-rbind(datend,df)} else {next}
  }
  distance<-rep(0,nrow(datend))# blank vector for distance
  for (i in 1:nrow(datend)) {distance[i]<-c[datend$IDinf[i],datend$IDsus[i]]} #find distance between 2 farms
  datend$distance<-distance
  
  popinf<-rep(0,nrow(datend))
  popsus<-rep(0,nrow(datend))
  for (i in 1:nrow(datend)) {
    popinf[i]<-x$pop[x$ID==datend$IDinf[i]] # find pop in each farm
    popsus[i]<-x$pop[x$ID==datend$IDsus[i]]
  }
  datend<-cbind(datend,popinf,popsus)
  datend$NIav<-datend$popinf/avpop
  datend$NSav<-datend$popsus/avpop
  
  return(datend)
}


######### kernel for recovered farms that became susceptible and escaping until the the end of outbreak  ########
prep_datreesc<-function (x,n,c){
  datreesc<-data.frame()
  for (i in 1:n) {
    a<-x$ID[which((x$startinfectiousdate<=i ) & (x$stopinfectiousdate>=i))]#a = IDinf will be infectious until stopinfect date
    b<-x$ID[which(x$waningdate<=i)] #b = IDsus (that immunity waning on that day)
    D<-i
    IDinf<-rep(a,each=length(b)) #to replicate it
    IDsus<-rep(b,length(a))
    df<-data.frame(cbind(D,IDinf,IDsus))
    if (length(df) == 3) {datreesc<-rbind(datreesc,df)} else {next}
  }
  distance<-rep(0,nrow(datreesc)) # blank vector for distance
  for (i in 1:nrow(datreesc)) {distance[i]<-c[datreesc$IDinf[i],datreesc$IDsus[i]]} #find distance between 2 farms
  datreesc$distance<-distance
  
  popinf<-rep(0,nrow(datreesc))
  popsus<-rep(0,nrow(datreesc))
  for (i in 1:nrow(datreesc)) {
    popinf[i]<-x$pop[x$ID==datreesc$IDinf[i]] # find pop in each farm
    popsus[i]<-x$pop[x$ID==datreesc$IDsus[i]]
  }
  datreesc<-cbind(datreesc,popinf,popsus)
  datreesc$NIav<-datreesc$popinf/avpop
  datreesc$NSav<-datreesc$popsus/avpop
  
  return(datreesc)
}

################  Run prepare data function ##########################

datesc_LP<-prep_datesc(x=FMDdata,n=max(FMDdata$stopinfectiousdate,na.rm=TRUE)-1, c=LPECdist) #x= outbreak data, c=distance matrix time max-1 because it is a group of escape until t infection
datinf_LP<-prep_datinf(x=FMDdata,n=max(FMDdata$stopinfectiousdate,na.rm=TRUE),c=LPECdist)
datend_LP<-prep_datend(x=FMDdata,n=max(FMDdata$stopinfectiousdate,na.rm=TRUE),c=LPECdist)
datreesc_LP<-prep_datreesc(x=FMDdata,n=max(FMDdata$stopinfectiousdate,na.rm=TRUE),c=LPECdist)

#D = day, IDinf = ID infectious farm, IDsus = IDsusceptble farm, distance = distance between farms (km), 
#popinf = population in the infectious farm, popsus = population in the susceptible farm, 
#NIav = popinf/average farm size, NSav= popsus/average farm size

#to quicken optimization time, I decide to change data.frame to data.table which process faster.
library(data.table)
datesc<-setDT(datesc_LP)
datinf<-setDT(datinf_LP)
datend<-setDT(datend_LP)
datreesc<-setDT(datreesc_LP)

################### Optimization process ################################
library(optimx)
library(bbmle)

#################################################
#### baseline kernel  ###########
################################################

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


fitll1CI <-mle2(minuslogl = ll1CI, start = list(k0 = 0.005, r0 = 0.2, alpha = 1.5), skip.hessian = FALSE)
# Here, I used optimizer "optimx" but there are other optimizers that you can use check ?mle2 for details.
Kfitll1CI<-profile(fitll1CI, alpha = 0.05) # profile likelihood
plot(Kfitll1CI)
confint(Kfitll1CI, level = 0.95)# to get profile CI


