################################################################################
## 03.05. 2007
## Lukas Gudmundsson
##
## simple implementation of Monte Carlo SSA MCSSA
################################################################################


.MCSSA <- function(U,lambda,x,n,conf=0.95,toeplitz=FALSE,keepSurr=FALSE,
                   ar.method="mle"){
### U: eigenvectors
### lambda: eigenvalues
### other parameters documented in package documentation
###
### basic, non visible MCSSA function
    if(ncol(U)!=nrow(U))
        stop("Matrix U is not square!")
    if(length(conf)==1){
        conf <- c((1-conf)/2,1-(1-conf)/2)
    } else if(length(conf)!=2 & sum(conf)!=1){
        stop("conf has more than two enties or the sum of the enties is not equal 1")
    }
    conf <- sort(conf)
    N <- length(x)
    M <- ncol(U)
    K <- N-M+1
    ar1<- ar(x,aic=FALSE,order.max=1,demean=TRUE,method=ar.method)
    rm(x)
    XS <- replicate(n,arima.sim(list(ar=ar1$ar),n=N,sd=sqrt(ar1$var.pred)))
    XS <- XS + ar1$x.mean
    if(toeplitz){
        .intFun <- function(x){
            S <- rep(NA,M)
            for(i in 1:M)
                S[i] <- x[1:(N-i+1)]%*%x[i:N]
            S <- S / (N-(1:M))
            S <- toeplitz(S)
            diag(t(U)%*%S%*%U)
        }
    } else {
        .intFun <- function(x){
            XX <- matrix(0,nrow=M,ncol=K)
            XX <- matrix(x[row(XX)+col(XX)-1],nrow=M,ncol=K)
            S <- XX %*% t(XX)
            diag(t(U)%*%S%*%U)
        }
    }
    XS <- apply(X=XS,MARGIN=2,.intFun)
    CL <- t(apply(X=XS,MARGIN=1,quantile,probs=conf,names=FALSE))
    freq <- apply(U,2,.findFreq)
    ord <- order(freq)
    .mct <- list(
                 lambda=lambda[ord],
                 freq=freq[ord],
                 rank=(1:M)[ord],
                 ar1=ar1,
                 upper=CL[,2][ord],
                 lower=CL[,1][ord],
                 conf=conf,
                 N=N,
                 M=M,
                 call=match.call()
                 )
    if(keepSurr)
        .mct$surLambda <- XS[ord,]
    class(.mct) <- "MCSSA"
    return(.mct)
}


MCSSA.decompSSA <- function(dSSA,x,n,conf=0.95,keepSurr=FALSE,
                            ar.method="mle",...){
### method for decompSSA
    .mct <- .MCSSA(U=dSSA$U,lambda=dSSA$lambda,
                   x=x,n=n,conf=conf,toeplitz=dSSA$toeplitz,keepSurr=keepSurr)
    .mct$seriesName <- dSSA$seriesName
    return(.mct)
}

MCSSA <- function(dSSA,...)
    UseMethod("MCSSA")
