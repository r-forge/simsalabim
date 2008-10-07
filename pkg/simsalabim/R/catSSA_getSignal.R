###################################################################
## Implementation of an algorithm for detecting trends and
## periodic components in SSA decompositions
##
## Lukas Gudmundsson
## 22.08.2008
###################################################################

##################################################################
## internal, not user visible functions
##################################################################
.periodogram <- function(e){
    L <- length(e)
    f <- seq.int(from=0,by=1/L,length.out=floor(L/2)+1)
    Y <- fft(e)
    Pyy <- Re(Y * Conj(Y))/L
    Pyy <- Pyy[1:(floor(L/2)+1)]
    return(list(pgrm=Pyy,freq=f))
}


##################################################################
## the main work  function's
##################################################################
getTrend <- function(U,omega0=1/nrow(U),C0){
    L <- nrow(U)
    U.pgrm <- apply(U,2,.periodogram)
    CC <- sapply(U.pgrm,function(pgrm){
        hh <- pgrm$freq<=omega0
        cc <- sum(pgrm$pgrm[hh])
        cc <- cc/sum(pgrm$pgrm)
        return(cc)
    })
    sel <- (1:L)[CC>=C0]
    output <- data.frame(k=sel,
                         smoothness=CC[CC>=C0]
                         )
    return(output)
}

getPeriod <- function(U,r0){
    L <- nrow(U)
    U.pgrm <- apply(U,2,.periodogram)
    domFreq <- sapply(U.pgrm,function(pgrm)
                      pgrm$freq[pgrm$pgrm==max(pgrm$pgrm)])
    deltaFreq <- abs(diff(domFreq))
    sel2 <- (1:(L-1))[deltaFreq<=1/L^2]
    sel2 <- as.list(sel2)
    RR <- sapply(sel2,function(k){
        ppk <- U.pgrm[[k]]$pgrm
        ppk1 <- U.pgrm[[k+1]]$pgrm
        ppk <- ppk/sum(ppk) ## amount of variance explained
        ppk1 <- ppk1/sum(ppk1)
        gamma <- ppk + ppk1
        rr <- gamma[1:((L/2)-1)]+gamma[2:(L/2)]
        rr <- max(rr/2)
        return(rr)
    })
    sel <- sel2[RR>=r0]
    sel <- as.numeric(sel)
    freq <- domFreq[sel]
    output <- data.frame(
                         k1=sel,
                         k2=sel+1,
                         freq=freq,
                         smootheness=RR[RR>=r0])
    return(output)
}


#####################################################################
## A wraper function for getTrend and getPeriod
#####################################################################
getSignal.decompSSA <- function(dSSA,omega0=1/dSSA$L,C0,r0,...){
    trends <- getTrend(dSSA$U,omega0,C0)
    periods <- getPeriod(dSSA$U,r0)
    isTrend <- trends[,1]%in%periods[,1] | trends[,1]%in%periods[,2]
    trendAndPeriod <- trends[isTrend,1]
    trends <- trends[!isTrend,]
    lambda <- dSSA$lambda
    EV <- lambda/sum(lambda)*100
    if(nrow(periods)==0){
        periods$explainedVariance <- numeric()
    } else {
        periods$explainedVariance <- sapply(1:nrow(periods),function(x)
                                            sum(EV[as.numeric(periods[x,1:2])]))
    }
    trends$explainedVariance <- EV[trends[,1]]
    output <- list(
                   periodic=periods,
                   trend=trends,
                   trendAndPeriod=trendAndPeriod,
                   call=match.call()
                   )
    class(output) <- "SSAsignal"
    return(output)
}

getSignal <- function(dSSA,...)
    UseMethod("getSignal")
