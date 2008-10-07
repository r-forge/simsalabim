################################################################################
## 19.08.2007
## Lukas Gudmundsson
##
## Performes an Itterative significance test on the Eigenspectrum of an single
## channes SSA
################################################################################

.sdTest <- function(x,y,conf,L,N,...){
### x : eigenvalues corresponding
### y: dominant frequency
### N : length of the initial ts
### L : embedding dimension
### conf : confidence levels
### y : dominant frequencies of EOFs
###
### performes Sun and Duffy test
    if(length(conf)==1){
        conf <- c((1-conf)/2,1-(1-conf)/2)
    } else if(length(conf)!=2 & sum(conf)!=1){
        stop("conf has more than two enties or the sum of the enties is not equal 1")
    }
    conf <- sort(conf)
    lnorm <- x/sum(x)*100
    dd <- data.frame(y,lnorm)
    ## remove elements with negative eigenvalues
    ## purpose: avoid NaNs in nls and lm!
    ## this assumes that negative eigenvalues are close to zero and due to
    ## numerical limitations.
    dd<-dd[dd$lnorm>0,]
    dd <- dd[order(dd$y),]
    rednoise <- log(lnorm)~log(A/(B+(2*pi*y)^2))
    nu <- round(3*N/L)
    conf <- sort(conf)
    qchi2 <- qchisq(conf,nu)
    for(i in 1:L){
        ## WARNING: estimate of starting values is very crude!!!
        spec.lm <- lm(lnorm*(2*pi*y)^2~lnorm,data=dd)
        sv <- abs(coef(spec.lm))
        names(sv)<-c("A","B")
        nmodel <- nls(rednoise,data=dd,start=sv)
        lfit <- exp(predict(nmodel))
        lfit <- lfit*sum(x)/100
        upper <- nu*lfit/qchi2[1]
        lower <- nu*lfit/qchi2[2]
        sig <- upper<dd$lnorm
        if(sum(sig)==0){
            break
        } else if(i==1){
            sigf <- dd$y[sig]
        } else {
            sigf <- c(sigf,dd$y[sig])
        }
        dd <- dd[!sig,]
    }
    lfit <- (exp(predict(nmodel,newdata=list(y=y,lnorm=lnorm)))*sum(x)/100)
    ord <- order(y)
    .sdt <- list(
                 lambda=x[ord],
                 freq=y[ord],
                 rank=(1:L)[ord],
                 fittedNoise=lfit[ord],
                 upper=(nu*lfit/qchi2[1])[ord],
                 lower=(nu*lfit/qchi2[2])[ord],
                 conf=conf,
                 N=N,
                 L=L,
                 call=match.call()
                 )
    class(.sdt) <- "sdTest"
    return(.sdt)
}


################################################################################
# decompSSA method for computing sdTest
################################################################################
sdTest.decompSSA <- function(dSSA,conf=0.95,...){
# dSSA : object of class SSA, output of SSA()
# conf : the confidence level

  if(is.null(dSSA$freq))
    dSSA$freq <-  apply(dSSA$U,2,.findFreq)

  #.sdt <- sdTest.default(x=dSSA$lambda,y=dSSA$freq,conf=conf,L=dSSA$L,N=dSSA$N)
  .sdt <- .sdTest(x=dSSA$lambda,y=dSSA$freq,conf=conf,L=dSSA$L,N=dSSA$N)
  .sdt$seriesName <- dSSA$seriesName
  return(.sdt)
}

### deactivated 22.08.2008
## sdTest.decompSSAM <- function(dSSA,conf=0.95,...){
##     if(is.null(dSSA$freq))
##         dSSA$freq <-  apply(dSSA$U,2,.findFreq)
##     .sdt <- .sdTest(x=dSSA$lambda,y=dSSA$freq,conf=conf,L=dSSA$L,N=dSSA$N)
##     .sdt$seriesName <- dSSA$seriesName
##     return(.sdt)
## }
sdTest <- function(dSSA,...)
    UseMethod("sdTest")

