#########################################################################
## Lukas Gudmundsson
## 31.8.2007
##
## Implementation of the basic SSA Algorithm
#########################################################################

decompSSA <- function(x,L,toeplitz=FALSE,getFreq=TRUE){
### SSA - decomposition
    N <- length(x)
    if(L>N/2) stop("L > N/2")
    K <- N-L+1

    if(toeplitz){
        S <- rep(NA,L)
        for(i in 1:L)
            S[i] <- x[1:(N-i+1)]%*%x[i:N]
        S <- S / (N-(1:L))
        S <- toeplitz(S)
    } else {
        S <- matrix(0,nrow=L,ncol=K)
        S <- matrix(x[row(S)+col(S)-1],nrow=L,ncol=K)
        S <- S%*%t(S)
    }

    ssa.dec <- eigen(S)
    names(ssa.dec) <- c("lambda","U")
    ssa.dec$rank <- 1:L
    if(getFreq)
        ssa.dec$freq <- apply(ssa.dec$U,2,.findFreq)
    ssa.dec$L <- L
    ssa.dec$N <- N
    ssa.dec$toeplitz=toeplitz
    ssa.dec$seriesName <- deparse(substitute(x))
    ssa.dec$call <- match.call()
    class(ssa.dec) <- "decompSSA"

    return(ssa.dec)
}

.groupSSA <- function(x,Us){
###  ds, Us : subset of a eigentripel
###
### SSA - grouping, see Golyandina et al
    if(!is.matrix(Us) & class(Us)=="numeric")
        Us <- matrix(Us)
    N <- length(x)
    L <- nrow(Us)
    K <- N-L+1
    XX <- matrix(0,nrow=L,ncol=K)
    XX <- matrix(x[row(XX)+col(XX)-1],nrow=L,ncol=K)
    ZZ <- matrix(0,nrow=L,ncol=K)
    for(ii in 1:ncol(Us))
        ZZ <- ZZ + Us[,ii]%*%t(t(XX)%*%Us[,ii])
    return(ZZ)
}

reconSSA <- function(dSSA,x,groups){
### dSSA : output of catSSA
### groups : list containing vectors determining the index of the eigentrippels
###          to be grouped together in one RC
### SSA - reconstruction
    N <- dSSA$N
    L <- dSSA$L
    K <- N-L+1
    RC <- matrix(0,nrow=N,ncol=length(groups))
    IND <- matrix(nrow=L,ncol=K)
    IND <- row(IND)+col(IND)-1
    for(ii in seq(along=groups)){
        XX <- .groupSSA(x,dSSA$U[,groups[[ii]]])
        ##Diagonal Averaging
        .intFun <- function(i,x,ind)
            mean(x[ind==i])
        RC[,ii] <- sapply(1:N,.intFun,x=XX,ind=IND)
    }
    colnames(RC) <- paste("RC",1:length(groups),sep="")
    return(RC)
}
