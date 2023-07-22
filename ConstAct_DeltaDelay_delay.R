#Antibodies affect the time delay between a successful activation and the time it takes for 
#vl to reach detection threshold. 

rm(list = ls())

library(cubature)
library(pracma)
library(Bhat)
library(readr)

### Extract data ###
Data1 <- read.csv("TimeToVR/Branching Process/LastUndetFirstDet_B2.csv", header=TRUE)
Data1 <- read.csv("TimeToVR/Branching Process/LastUndetFirstDet_C5.csv", header=TRUE)
Data3 <- read.csv("TimeToVR/Branching Process/LastUndetFirstDet_C4.csv", header=TRUE)


# FirstDet <- c(Data1$FirstDet,Data2$FirstDet,Data3$FirstDet)
# LastUndet <- c(Data1$LastUndet,Data2$LastUndet,Data3$LastUndet)
# Neut_Percent <- c(Data1$Neut_Percent,Data2$Neut_Percent,Data3$Neut_Percent)/100

FirstDet <- Data1$FirstDet
LastUndet <- Data1$LastUndet
Neut_Percent <- Data1$Neut_Percent

animals<-length(FirstDet)

## Keep data up to 60 days
## For monkeys with FirstDet>60, set LastUndet=60 and FirstDet=NA
LastUndet[which(LastUndet>60)]=60
FirstDet[which(FirstDet>60)]=NA

f_det=44/100


AvgPnorbd <- function(f,a,delay1,delay2,t){#no rebound function
  #delay=delay1+ (delay2-delay1)*f#^2
  delay=delay1+ delay2*f
  #delay=delay1+ (delay1+delay2)*f
  return(exp(-a*(t-delay))*1/f_det)
}

### Negative Log-likelihood for an individual
### derived from cumulative probability of rebound in [t-1,t]
nloglik <- function(params,FirstDet,LastUndet, Neut_Percent) {
  # for an individual
  a = params[1]
  delay1 = params[2]
  delay2=params[3]
  
  t1 = FirstDet
  t2 = LastUndet
  f = Neut_Percent
  
  #Compute pvr1
  if (is.na(t1)) {
    pvr1 = 1
  } else {
    if(is.na(f)){
      pvr1 <- 1 - integrate(AvgPnorbd, lower=0, upper=f_det, a=params[1],
                            delay1=params[2],delay2=params[3],t=t1)$value
    } else {
      #delay=delay1+ (delay2-delay1)*f#^2
      delay=delay1+ delay2*f
      #delay=delay1+ (delay1+delay2)*f
      pvr1 <- 1 - exp(-a*(t1-delay)) 
    }
  }
  if (pvr1<0) {pvr1=0}
  
  #Compute pvr2
  if(is.na(f)){
    pvr2 <- 1 - integrate(AvgPnorbd, lower=0, upper=f_det, a=params[1],
                          delay1=params[2],delay2=params[3],t=t2)$value
  } else {
    #delay=delay1+ (delay2-delay1)*f#^2
    delay=delay1+ delay2*f
    #delay=delay1+ (delay1+delay2)*f
    pvr2 <- 1 - exp(-a*(t2-delay))  
  }
  if (pvr2<0) {pvr2=0}
  
  # careful about the possibility of numerical junk when pvr1, pvr2 very near 1.
  if (pvr1<=pvr2) {
    lpvr <- 1e16  
  } else {
    lpvr <- -log(pvr1-pvr2)
  }
  return(lpvr)
}

nloglikall <- function(params) {
  ### negative log likelihood summed over all data 
  nlogfsum <- 0
  for (jj in 1:length(FirstDet)) {
    nlogfsum <- nlogfsum + nloglik(params,FirstDet[jj],LastUndet[jj], Neut_Percent[jj])
  }
  return(nlogfsum)
}

### optimization time
x0de <- list(label=c("a", "delay1", "delay2"),est=c(0.1,2,3),
             low=c(1e-14, 1e-14, 1e-14),upp=c(5, 20, 40))
### Run optimization dfp--Davidon–Fletcher–Powell algorithm``
mlede=dfp(x0de, nloglikall)

############# log-likelihood & AICs/BIC
nlogL <- nloglikall(mlede$est)
AIC = 2*length(mlede$est) - 2*(-nlogL) 
BIC = -2*(-nlogL) + length(mlede$est)*log(length(FirstDet)) 

AIC
BIC


#Plot
#grey area-- empirical pvr
#black solid line-- median pvr
f<-Neut_Percent
a<-mlede$est[1]
delay1<-mlede$est[2]
delay2<-mlede$est[3]

firstdet_alldata <- NULL
lastundet_alldata <- NULL
# model
ts<-seq(0,60,by=1)

pvr<-as.data.frame(matrix(data=NA, nrow=length(ts), ncol=10)) # NULL
for (jj in c(1:length(ts))) {
  # empirical
  firstdet_alldata[jj] <- length(which(FirstDet<=ts[jj]))/length(FirstDet)
  lastundet_alldata[jj] <- length(which(LastUndet<=ts[jj]))/length(LastUndet)
  # model - MAKE MULTIPLE CURVES, COMPUTE MEDIAN FROM THIS.
  
  for(i in 1:animals){
    if(is.na(f[i])){ #no neutralization
      pvr[jj, i] <- 1 - integrate(AvgPnorbd, lower=0, upper=f_det, a=a,
                                  delay1=delay1,delay2=delay2,t=ts[jj])$value
      if (pvr[jj, i]<0) {pvr[jj, i]=0}# if t<delay, pvr=0 for that monkey
    }
    else{ #there is neutralization
      #delay=delay1+(delay2-delay1)*f[i]#^2
      delay=delay1+delay2*f[i]
      #delay=delay1+(delay1+delay2)*f[i]
      
      pvr[jj, i] <- 1 - exp(-a*(ts[jj]-delay))
      if (ts[jj]<=delay) {pvr[jj,i]=0} # if t<delay, pvr=0 for that monkey
    }
  }
}

median_vector<-NULL  #the median value of pvr at time ts

for (i in 1:nrow(pvr)){
  x<-as.numeric(pvr[i, ])
  median_vector[i]=median(x)
}

x<-which(median_vector<0)
median_vector[x]=0

plot(ts, firstdet_alldata, type="l", xlim=c(0, 60), ylim=c(0, 1),
     xlab="Time since ATI (days)", ylab="Cum. Prob. Rbd")
par(new=TRUE)
plot(ts, lastundet_alldata, type="l", xlim=c(0, 60), ylim=c(0, 1),
     xlab="", ylab="")
polygon(c(ts, rev(ts)), c(lastundet_alldata, rev(firstdet_alldata)),
        col = "grey", border = NA)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ts, median_vector, type="l", col="red", xlim=c(0, 60), ylim=c(0, 1), xlab="", ylab="")
legend("topleft", legend=c("Empirical", "Model"), col=c("grey", "red"), lty=1, cex=0.8 )