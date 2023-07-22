rm(list = ls())

library(cubature)
library(pracma)
library(Bhat)
library(readr)

### Extract data ###
Data1 <- read.csv("TimeToVR/Branching Process/LastUndetFirstDet_B2.csv", header=TRUE)
Data2 <- read.csv("TimeToVR/Branching Process/LastUndetFirstDet_C5.csv", header=TRUE)
Data3 <- read.csv("TimeToVR/Branching Process/LastUndetFirstDet_C4.csv", header=TRUE)

FirstDet <- c(Data1$FirstDet, Data2$FirstDet, Data3$FirstDet)
LastUndet <- c(Data1$LastUndet, Data2$LastUndet, Data3$LastUndet)
Neut_Percent <- c(Data1$Neut_Percent, Data2$Neut_Percent, Data3$Neut_Percent)/100
CD4_Data<-c(Data1$CD4_Wk0_ul,Data2$CD4_Wk0_ul,Data3$CD4_Wk0_ul)
IPDA_Data<-c(Data1$IPDA, Data2$IPDA, Data3$IPDA)
Det_IPDA_Data=c(Data1$Det_IPDA, Data2$Det_IPDA, Data3$Det_IPDA)

# FirstDet <- Data1$FirstDet
# LastUndet <- Data1$LastUndet
# Neut_Percent <- Data1$Neut_Percent/100
# CD4_Data<-Data1$CD4_Wk0_ul
# IPDA_Data<-Data1$IPDA
# Det_IPDA_Data= Data1$Det_IPDA

f_det=44/100  #antibody neutralization detection threshold
#parameters of the lognormal distribution of CD4 T cells 
mu=(19.33)
sigma=(0.53)
CD4_Data<-CD4_Data*1e3*162 #over the entire monkey 


animals=length(FirstDet)

## Keep data up to 60 days
## For monkeys with FirstDet>60, set LastUndet=60 and FirstDet=NA
LastUndet[which(LastUndet>60)]=60 
FirstDet[which(FirstDet>60)]=NA

n=1
AvgPnorbd_CD4_f <- function(x,IPDA,q0,delay1,delay2,k,a,t){# we have IPDA data but no CD4
  f=x[1]
  lCD4= x[2]
  LR=log10(exp(lCD4)*IPDA/1e6)  #LR size over the entire monkey, take the log10
  pdf_CD4=(1/(sigma*sqrt(2*pi)))*exp(-(((lCD4)-mu)^2)/(2*sigma^2)) #pdf of lognormal for CD4s
  q=min(q0*f,1)
  delay=delay1+delay2*f^n/(k^n+f^n)
  return(exp(-a*(1-q)*LR*(t-delay))*pdf_CD4*(1/f_det))
}

AvgPnorbd_CD4_IPDA_f <- function(x,Det_IPDA,q0,delay1,delay2,k,a,t){ 
  f=x[1]
  IPDA=x[2]
  lCD4=x[3]
  LR=log10(exp(lCD4)*IPDA/1e6)  
  pdf_CD4=(1/(sigma*sqrt(2*pi)))*exp(-(((lCD4)-mu)^2)/(2*sigma^2)) 
  q=min(q0*f,1)
  delay=delay1+delay2*f^n/(k^n+f^n)
  return(exp(-a*(1-q)*LR*(t-delay))*pdf_CD4*(1/Det_IPDA)*(1/f_det))
}

### Negative Log-likelihood for an individual
### derived from cumulative probability of rebound in [t-1,t]
nloglik <- function(params,FirstDet,LastUndet, IPDA_Data, CD4_Data, 
                    Neut_Percent,Det_IPDA_Data) {
  # for an individual
  q0 = params[1]
  delay1 = params[2]
  delay2 = params[3]
  k = params[4]
  a = params[5]
  
  t1 = FirstDet
  t2 = LastUndet
  IPDA = IPDA_Data
  CD4 = CD4_Data
  f = Neut_Percent
  Det_IPDA = Det_IPDA_Data
  
  if (is.na(t1)) {
    pvr1 = 1
  } else {
    if(is.na(CD4) & !is.na(IPDA) & is.na(f)){ #no CD4 Data, yes detectable IPDA, no neut Abs-- double integral
      pvr1 <-  1- adaptIntegrate(AvgPnorbd_CD4_f, lowerLimit=c(0,15),
                                 upperLimit=c(f_det,25),IPDA=IPDA,
                                 q0=params[1],delay1=params[2],delay2=params[3],
                                 k=params[4],a=params[5],t=t1)$integral
    }else if(is.na(CD4) & is.na(IPDA) & is.na(f)){#no CD4 Data, detectable IPDA or neut Abs-- triple integral
      pvr1 <- 1- adaptIntegrate(AvgPnorbd_CD4_IPDA_f, lowerLimit=c(0,1e-2,15),
                                upperLimit=c(f_det,Det_IPDA,25),Det_IPDA=Det_IPDA,
                                q0=params[1],delay1=params[2],delay2=params[3],
                                k=params[4],a=params[5],t=t1)$integral
    }else{ #yes CD4 Data, detectable IPDA, and neut Abs-- no integral
      LR=log10(CD4*IPDA/1e6)
      q=min(q0*f,1)
      delay=delay1+delay2*f^n/(k^n+f^n)
      pvr1 <- 1 - exp(-LR*a*(1-q)*(t1-delay))
    }
  }
  if (pvr1<0){pvr1=0}
  
  
  if(is.na(CD4) & !is.na(IPDA) & is.na(f)){ #IPDA, no CD4 or f-- double integral 
    pvr2 <- 1-adaptIntegrate(AvgPnorbd_CD4_f, lowerLimit=c(0,15),
                             upperLimit=c(f_det,25),IPDA=IPDA,
                             q0=params[1],delay1=params[2],delay2=params[3],
                             k=params[4],a=params[5],t=t2)$integral
  }else if(is.na(CD4) & is.na(IPDA) & is.na(f)){#no IPDA, CD4 or f-- triple integral 
    pvr2 <- 1- adaptIntegrate(AvgPnorbd_CD4_IPDA_f, lowerLimit=c(0,1e-2,15),
                              upperLimit=c(f_det,Det_IPDA,25),Det_IPDA=Det_IPDA,
                              q0=params[1],delay1=params[2],delay2=params[3],
                              k=params[4],a=params[5],t=t2)$integral
  }else{
    LR=log10(CD4*IPDA/1e6)
    q=min(q0*f,1)
    delay=delay1+delay2*f^n/(k^n+f^n) 
    pvr2 <- 1 - exp(-LR*(1-q)*a*(t2-delay))
  }
  if (pvr2<0){pvr2=0}
  
  
  #there are no cases where we have f but either IPDA or CD4 are missing 
  
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
  #  print(params)
  for (jj in 1:length(FirstDet)) {
    nlogfsum <- nlogfsum + nloglik(params,FirstDet[jj],LastUndet[jj], IPDA_Data[jj],
                                   CD4_Data[jj],Neut_Percent[jj],Det_IPDA_Data[jj])
    #  print(jj)    
  }
  return(nlogfsum)
}

### optimization time
x0de <- list(label=c("q0","delay1","delay2","k","a"),est=c(0.1,5,2,0.01,0.01),
             low=c(1e-14,1e-14,1e-14,1e-14,1e-14),upp=c(2,300,300,100,20))

### Run optimization dfp
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
CD4<-CD4_Data
IPDA<-IPDA_Data
Det_IPDA<-Det_IPDA_Data
LR<-log10(CD4*IPDA/1e6)
q0<-mlede$est[1]
delay1<-mlede$est[2]
delay2<-mlede$est[3]
k<-mlede$est[4]
a<-mlede$est[5]

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
    if(is.na(CD4[i]) & !is.na(IPDA[i]) & is.na(f[i])){ #IPDA, no CD4 or f-- double integral 
      pvr[jj, i] <- 1-adaptIntegrate(AvgPnorbd_CD4_f, lowerLimit=c(0,15),
                                     upperLimit=c(f_det,25),IPDA=IPDA[i],
                                     q0=q0,delay1=delay1,delay2=delay2,k=k,a=a,t=ts[jj])$integral
      if (pvr[jj,i]<0) {pvr[jj,i]=0}
    }else if(is.na(CD4[i]) & is.na(IPDA[i]) & is.na(f[i])){#no IPDA, CD4 or f-- triple integral 
      pvr[jj, i] <- 1- adaptIntegrate(AvgPnorbd_CD4_IPDA_f, lowerLimit=c(0,1e-2,15),
                                      upperLimit=c(f_det,Det_IPDA[i],25),Det_IPDA=Det_IPDA[i],
                                      q0=q0,delay1=delay1,delay2=delay2,k=k,a=a,t=ts[jj])$integral
      if (pvr[jj,i]<0) {pvr[jj,i]=0}
    }else{ #there is neutralization
      q=min(q0*f[i],1)
      delay=delay1+delay2*f[i]^n/(k^n+f[i]^n) 
      pvr[jj, i] <- 1 - exp(-LR[i]*a*(1-q)*(ts[jj]-delay))
      if (ts[jj]<=delay) {pvr[jj,i]<-0} # if t<delay, pvr=0 for that monkey
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
