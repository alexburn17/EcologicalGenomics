# P. Alexander Burnham
# first 35 lines written by Lauren Ash
# Feb 13, 2017


source("https://bioconductor.org/biocLite.R")
biocLite("RNASeqPower")
library(RNASeqPower)


rep2<-rnapower(depth=10, n=2, cv=0.2, effect=c(1.5, 2, 2.5, 4), alpha=.05)

plot(rep2, type='b', col="blue", xlab="Effect size (fold change)",ylab="Power",ylim=c(0,1), lty=2)

rep3<-rnapower(depth=10, n=3, cv=0.2, effect=c(1.5, 2, 2.5, 4), alpha=.05)
points(rep3, type='b', col="purple", lty=2)

rep5<-rnapower(depth=10, n=5, cv=0.2, effect=c(1.5, 2, 2.5, 4), alpha=.05)
points(rep5, type='b', col="green", lty=2)

rep10<-rnapower(depth=10, n=10, cv=0.2, effect=c(1.5, 2, 2.5, 4), alpha=.05)
points(rep10, type='b', col="red", lty=2)

rep2V<-rnapower(depth=10, n=2, cv=0.6, effect=c(1.5, 2, 2.5, 4), alpha=.05)
points(rep2V, type='b', col="blue")

rep3V<-rnapower(depth=10, n=3, cv=0.6, effect=c(1.5, 2, 2.5, 4), alpha=.05)
points(rep3V, type='b', col="purple")

rep5V<-rnapower(depth=10, n=5, cv=0.6, effect=c(1.5, 2, 2.5, 4), alpha=.05)
points(rep5V, type='b', col="green")

rep10V<-rnapower(depth=10, n=10, cv=0.6, effect=c(1.5, 2, 2.5, 4), alpha=.05)
points(rep10V, type='b', col="red")

#--------------------------------------------------------------------------------

Alexrep1 <- rnapower(depth=10, n=20, cv=0.5, effect=c(1.5, 2, 2.5, 4), alpha=.05)

plot(Alexrep1, type='b', col="blue", xlab="Effect size (fold change)",ylab="Power",ylim=c(0,1), lty=2)






