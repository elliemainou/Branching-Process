#Antibodies affect the time delay between a successful activation and the time it takes for 
#vl to reach detection threshold. 

rm(list = ls())

library(cubature)
library(pracma)
library(Bhat)
library(readr)

### Extract data ###
#Data1 <- read.csv("TimeToVR/Branching Process/LastUndetFirstDet_B2.csv", header=TRUE)
Data1 <- read.csv("TimeToVR/Branching Process/LastUndetFirstDet_C5.csv", header=TRUE)
Data1 <- read.csv("TimeToVR/Branching Process/LastUndetFirstDet_C4.csv", header=TRUE)

# FirstDet <- c(Data1$FirstDet,Data2$FirstDet,Data3$FirstDet)
# LastUndet <- c(Data1$LastUndet,Data2$LastUndet,Data3$LastUndet)
# IPDA_Data<-c(Data1$IPDA,Data2$IPDA,Data3$IPDA)
# Det_IPDA_Data= c(Data1$Det_IPDA,Data2$Det_IPDA,Data3$Det_IPDA)
# CD4_Data<-c(Data1$CD4_Wk0_ul,Data2$CD4_Wk0_ul,Data3$CD4_Wk0_ul)

FirstDet <- Data1$FirstDet
LastUndet <- Data1$LastUndet
IPDA_Data<-Data1$IPDA
Det_IPDA_Data= Data1$Det_IPDA
CD4_Data<-Data1$CD4_Wk0_ul

mu=19.33
sigma=0.53
CD4_Data<-CD4_Data*1e3*162 #over the entire monkey 

animals<-length(FirstDet)

## Keep data up to 60 days
## For monkeys with FirstDet>60, set LastUndet=60 and FirstDet=NA
LastUndet[which(LastUndet>60)]=60
FirstDet[which(FirstDet>60)]=NA

AvgPnorbd_CD4 <- function(IPDA,lCD4,a,delay,t){# we have IPDA data but no CD4
  LR=log10(exp(lCD4)*IPDA/1e6)
  pdf_CD4=(1/(sigma*sqrt(2*pi)))*exp(-(((lCD4)-mu)^2)/(2*sigma^2))
  return(exp(-a*LR*(t-delay))*pdf_CD4) #a=(1-q)*a
}

AvgPnorbd_CD4_IPDA <- function(x,Det_IPDA,a,delay,t){
  IPDA=x[1]
  lCD4=x[2]
  LR=log10(exp(lCD4)*IPDA/1e6)
  pdf_CD4=(1/(sigma*sqrt(2*pi)))*exp(-(((lCD4)-mu)^2)/(2*sigma^2))
  return(exp(-a*LR*(t-delay))*pdf_CD4*(1/Det_IPDA)) #a=(1-q)*a
}

### Negative Log-likelihood for an individual
### derived from cumulative probability of rebound in [t-1,t]
nloglik <- function(params,FirstDet,LastUndet,IPDA_Data,CD4_Data,Det_IPDA_Data) {
  # for an individual
  a = params[1]
  delay = params[2]
  
  t1 = FirstDet
  t2 = LastUndet
  IPDA = IPDA_Data
  CD4 = CD4_Data
  Det_IPDA=Det_IPDA_Data
  
  
  if (is.na(t1)) {
    pvr1 = 1
  } else if (t1-delay<0) {
    pvr1 <- 0
  } else {
    if(is.na(CD4) & !is.na(IPDA)){
      pvr1 <- 1 - integrate(AvgPnorbd_CD4, lower=15, upper=25, a=params[1],
                            delay=params[2], IPDA=IPDA_Data,t=t1)$value  
    } else if(is.na(CD4) & is.na(IPDA)){
      pvr1 <- 1-adaptIntegrate(AvgPnorbd_CD4_IPDA, lowerLimit = c(0.01, 15), 
                               upperLimit = c(Det_IPDA_Data,25), a=params[1],
                               delay=params[2], Det_IPDA=Det_IPDA_Data,t=t1)$integral
    }else{
      LR=log10(IPDA*CD4/1e6)
      pvr1 <- 1 - exp(-a*LR*(t1-delay)) 
    }
  }
  if (t2-delay<0) {
    pvr2 = 0
  } else {
    if(is.na(CD4_Data) & !is.na(IPDA_Data)){
      pvr2 <- 1 - integrate(AvgPnorbd_CD4, lower=15, upper=25, a=params[1],
                            delay=params[2], IPDA=IPDA_Data,t=t2)$value
    } else if(is.na(CD4_Data) & is.na(IPDA_Data)){
      pvr2 <- 1-adaptIntegrate(AvgPnorbd_CD4_IPDA, lowerLimit = c(0.01, 15), 
                               upperLimit = c(Det_IPDA_Data,25), a=params[1],
                               delay=params[2], Det_IPDA=Det_IPDA_Data,t=t2)$integral
    } else {
      LR=log10(IPDA*CD4/1e6)
      pvr2 <- 1 - exp(-a*LR*(t2-delay))  
    }
  }  
  
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
    nlogfsum <- nlogfsum + nloglik(params,FirstDet[jj],LastUndet[jj], 
                                   IPDA_Data[jj],CD4_Data[jj],Det_IPDA_Data[jj])
  }
  return(nlogfsum)
}

### optimization time
x0de <- list(label=c("a","delay"),est=c(0.1,10),
             low=c(1e-14, 1e-14),upp=c(10, 40))
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
IPDA<-IPDA_Data
Det_IPDA<-Det_IPDA_Data
CD4<-CD4_Data

a<-mlede$est[1]
delay<-mlede$est[2]

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
    if(is.na(CD4_Data[i]) & !is.na(IPDA_Data[i])){ #no CD4 and IPDA
      pvr[jj, i] <- 1 - integrate(AvgPnorbd_CD4, lower=15, upper=25, a=a,
                                  delay=delay, IPDA=IPDA_Data[i],t=ts[jj])$value
      if (ts[jj]<=delay) {pvr[jj, i]<-0}# if t<delay, pvr=0 for that monkey
    }else if (is.na(CD4_Data[i]) & is.na(IPDA_Data[i])){
      pvr[jj, i] <-1-adaptIntegrate(AvgPnorbd_CD4_IPDA, lowerLimit = c(0.01, 15), 
                                    upperLimit = c(Det_IPDA_Data[i],25), a=a,
                                    delay=delay, Det_IPDA=Det_IPDA_Data[i],t=ts[jj])$integral
    }
    else{ #there are both CD4 and IPDA 
      LR=log10(IPDA[i]*CD4[i]/1e6)
      pvr[jj, i] <- 1 - exp(-a*LR*(ts[jj]-delay))
      if (ts[jj]<=delay) {pvr[jj,i]<-0} # if t<delay, pvr=0 for that monkey
    }
  }
}


median_vector<-NULL  #the median value of pvr at time ts

for (i in 1:nrow(pvr)){
  x<-as.numeric(pvr[i, ])
  median_vector[i]=median(x)
}


median_vector[1:(delay+1)]=0

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
