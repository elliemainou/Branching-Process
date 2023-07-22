#Fit baseline prob of unsuccessful activitation 
library(cubature)
library(pracma)
library(Bhat)
library(readr)

### Extract data ###
# Data1 <- read.csv("TimeToVR/Branching Process/LastUndetFirstDet_B2.csv", header=TRUE)
# Data2 <- read.csv("TimeToVR/Branching Process/LastUndetFirstDet_C5.csv", header=TRUE)
Data1 <- read.csv("TimeToVR/Branching Process/LastUndetFirstDet_C4.csv", header=TRUE)

#All groups
# FirstDet <- c(Data1$FirstDet,Data2$FirstDet,Data3$FirstDet)
# LastUndet <- c(Data1$LastUndet,Data2$LastUndet,Data3$LastUndet)

# Each group
FirstDet <- Data1$FirstDet
LastUndet <- Data1$LastUndet

animals<-length(FirstDet)

## Keep data up to 60 days
## For monkeys with FirstDet>60, set LastUndet=60 and FirstDet=NA
LastUndet[which(LastUndet>60)]=60
FirstDet[which(FirstDet>60)]<-NA

### Negative Log-likelihood for an individual
### derived from cumulative probability of rebound in [t-1,t]
nloglik <- function(params,FirstDet,LastUndet) {
  # for an individual
  a = params[1] 
  delay = params[2]
  t1 = FirstDet
  t2 = LastUndet

  if (is.na(t1)) {
    pvr1 = 1
  } else if (t1-delay<0) {
    pvr1 <- 0
  } else {
    pvr1 <- 1 - exp(-a*(t1-delay)) #
  }
  if (t2-delay<0) {
    pvr2 = 0
  } else {
    pvr2 <- 1 - exp(-a*(t2-delay)) #
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
    nlogfsum <- nlogfsum + nloglik(params,FirstDet[jj],LastUndet[jj])
  }
  return(nlogfsum)
}

### optimization time
#x0de <- list(label=c("a","delay"),est=c(0.02,5),low=c(1e-14,1e-14),upp=c(5,30))

#for C4
x0de <- list(label=c("a","delay"),est=c(0.02,20),low=c(1e-14,1e-14),upp=c(5,100))

### Run optimization dfp
mlede=dfp(x0de, nloglikall)

# ############# log-likelihood & AICs/BIC
nlogL <- nloglikall(mlede$est)
AIC = 2*length(mlede$est) - 2*(-nlogL)
BIC = -2*(-nlogL) + length(mlede$est)*log(length(FirstDet)) # Bayes Information Criterion

AIC
BIC

#Plot
#grey area-- empirical pvr
#black solid line-- median pvr
a<-mlede$est[1]
delay<-mlede$est[2]

firstdet_alldata <- NULL
lastundet_alldata <- NULL
# model
n=60
ts<-seq(0,n,by=1)


pvr<-as.data.frame(matrix(data=NA, nrow=length(ts), ncol=10)) # NULL 
for (jj in c(1:length(ts))) {
  # empirical
  firstdet_alldata[jj] <- length(which(FirstDet<=ts[jj]))/length(FirstDet)
  lastundet_alldata[jj] <- length(which(LastUndet<=ts[jj]))/length(LastUndet)
  # model - MAKE MULTIPLE CURVES, COMPUTE MEDIAN FROM THIS.
  for(i in 1:animals){
    pvr[jj, i] <- 1 - exp(-a*(ts[jj]-delay))
    }
  
  if (ts[jj]<=mlede$est[2]) {pvr[jj]<-0}
}

median_vector<-NULL  #the median value of pvr at time ts

for (i in 1:nrow(pvr)){
  x<-as.numeric(pvr[i, ])
  median_vector[i]=median(x)
}

#for ts<delay, pvr=0 
median_vector[1:(delay+1)]=0

plot(ts, firstdet_alldata, type="l", xlim=c(0, n), ylim=c(0, 1), 
     xlab="Time since ATI (days)", ylab="Cum. Prob. Rbd")
par(new=TRUE)
plot(ts, lastundet_alldata, type="l", xlim=c(0, n), ylim=c(0, 1),
     xlab="", ylab="")
polygon(c(ts, rev(ts)), c(lastundet_alldata, rev(firstdet_alldata)),
        col = "grey", border = NA)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ts, median_vector, type="l", col="red", xlim=c(0, n), ylim=c(0, 1), xlab="", ylab="")
legend("topleft", legend=c("Empirical", "Model"), col=c("grey", "red"), lty=1, cex=0.8 )

