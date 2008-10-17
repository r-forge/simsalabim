################################################################################
## 19.04. 2007
## Lukas Gudmundsson
##
## Implementation ofthe SSA gapfilling procedure presented by Kondrashov and
## Ghil (Nonlin. Processes Geophys., 13, 151–159, 2006)
##
## modified
## 10.10.2008: various bugs and inconsistencies fixed
## 15.10.2008: updated cross validation.
################################################################################

kgSSAM <- function(x,L,toeplitz=FALSE,K,max.inner=100,eps=0.05,
                   getFreq = TRUE,verbose=FALSE){
   mpos <- is.na(x)
   x[!mpos]<-x[!mpos]-mean(x,na.rm=TRUE)
   x[mpos]<- 0
   nrmse <- rep(NA,K)
   names(nrmse) <- 1:K
   xRC <- matrix(NA,length(x),K)
   colnames(xRC) <- paste("RC",1:K,sep="")
   for(k in 1:K){
       for(i in 1:max.inner){
           x.ssa <- decompSSA(x=x,L=L,toeplitz=toeplitz,getFreq=FALSE)
           xRC[,k] <- reconSSA(dSSA=x.ssa,x=x,groups=list(1:k))
           nrmse[k] <- sqrt(mean((x[mpos]-xRC[mpos,k])^2))/sqrt(mean(x[mpos]^2))
           x[mpos] <- xRC[mpos,k]
           if(nrmse[k]<=eps)
               break
       }
       if(verbose)
           cat(date(),"comp",k,"of",K,"iter:",i,"\n")
   }
   x.ssa <- decompSSA(x=x,L=L,toeplitz=toeplitz,getFreq=getFreq)
   x.ssa$RC <- xRC
   x.ssa$nrmse <- nrmse
   x.ssa$seriesName <- deparse(substitute(x))
   x.ssa$call <- match.call()
   class(x.ssa) <- "kgSSAM"
   return(x.ssa)
}

cvSSA.kgSSAM <- function(x,L,K=10,cv.crit="bestApprox",cv.gap=0.05,
                         runs=200,verboseCV=TRUE,...){
### Parameter Estimation based on cross validation
### values vor cv.crit: "betsApprox" ; "stabSignal"
   nrmse <- matrix(NA,K,runs)
   nmpos <- !is.na(x)
   nmpos <- (1:length(x))[nmpos]
   cvn <- round(cv.gap*length(nmpos))
   origKgSSAM <- kgSSAM(x,L=L,K=K,...)
   bench <- switch(cv.crit,
                   "stabSignal"=origKgSSAM$RC,
                   "bestApprox"=matrix(x, ncol=1, nrow=length(x)),
                   stop("No valid performance criterion")
                   )
   for(r in 1:runs){
       if(verboseCV)
           cat(date(),"loop",r,"of",runs,"\n")
       sam <- sample(nmpos,cvn)
       xx <- x
       xx[sam] <- NA
       RC.cv <- kgSSAM(xx,L=L,K=K,...)$RC
       nrmse[,r] <- sqrt(colMeans((bench[sam,]-RC.cv[sam,])^2))
   }
   colnames(nrmse) <- paste("run",1:runs,sep="")
   rownames(nrmse) <- paste("RC",1:K,sep="")
   origKgSSAM$CVerror <- nrmse
   origKgSSAM$CVcall <- match.call()
   origKgSSAM$class <- c("kgSSAM","cvSSA")
   return(origKgSSAM)
}



