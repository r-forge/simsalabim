################################################################################
## 15.04. 2007
## Lukas Gudmundsson
##
## omega - correlations
################################################################################

w.cor <- function(R,L){
    if(missing(L)){
        L <- dim(R)[2]
        warning("L not specified assuming complete set of RCs setting L to ",L)
    }
    N <- dim(R)[1]
    K <- N - L +1
    Ls <- min(L,K)
    Ks <- min(L,K)
    LL <- dim(R)[2]
    WC <- matrix(NA,nrow=LL,ncol=LL)
    ww <- rep.int(NA,N)
    if(Ls!=Ks){
        ww[1:(Ls-1)] <- 1:(Ls-1)
        ww[Ls:(Ks-1)] <- Ls
        ww[Ks:N] <- N-(Ks:N)+1
    } else {
        ww[1:(Ls-1)] <- 1:(Ls-1)
        ww[Ls:Ks] <- Ls
        ww[(Ks+1):N] <- N-((Ks+1):N)+1
    }
    .intFun<-function(k,R,ww)
        (ww*R[,k])%*%R
    WC <- sapply(1:LL,.intFun,R=R,ww=ww)
    dd <- sqrt(diag(WC))
    WC <- WC/dd%o%dd
    dimnames(WC) <- list(colnames(R),colnames(R))
    return(WC)
}
