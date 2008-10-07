################################################################################
## Lukas Gudmdundsson
## 1.9.2007
##
## Implementation of:
## The “Caterpillar”-SSA method for analysis oftime series with missing values
## by Golyandina and Osipov
################################################################################

.star <- function(a,b){
### A formal substitution for the inner product for vectors
### a & b containing missing values
  xx <- cbind(a,b)
  n <- dim(xx)[1]
  xx <- na.omit(xx)
  return((xx[,1]%*%xx[,2])*(n/(nrow(xx))))
}


################################################################################
# Stage1: Decomposition
################################################################################

decompSSAM <- function(x,L,tau=0,toeplitz=FALSE,getFreq=TRUE){
### SSA - decomposition for time series with missing values
    N <- length(x)
    if(L>N/2) stop("L > N/2")
    K <- N-L+1
    if(tau<=0 & toeplitz==FALSE ){
        XX <- matrix(0,nrow=L,ncol=K)
        XX <- matrix(x[row(XX)+col(XX)-1],nrow=L,ncol=K)
        CC <- !is.na(colSums(XX))
        S <- XX[,CC]%*%t(XX[,CC])
    } else if(tau>0 & toeplitz==FALSE ){
        XX <- matrix(0,nrow=L,ncol=K)
        XX <- matrix(x[row(XX)+col(XX)-1],nrow=L,ncol=K)
        CC <- colSums(is.na(XX))/L
        CC <- CC<=tau
        XX <- XX[,CC]
        S <- matrix(0,L,L)
        for(ii in 1:L){
            for(jj in 2:L)
                S[c(ii,jj),c(jj,ii)] <- .star(XX[ii,],XX[jj,])
        }
        for(ii in 1:L)
            S[ii,ii] <- .star(XX[ii,],XX[ii,])
    } else {
        S <- rep(NA,L)
        for(i in 1:L)
            S[i] <- .star(x[1:(N-i+1)],x[i:N])
        S <- S/(1+N-(1:L))
        S <- toeplitz(S)
    }
    UU <- eigen(S)
    names(UU) <- c("lambda","U")
    ## crude way to account for possibly negative/zero eigenvalues
    ## assuming that these are due to numerical limitations
    UU$lambda[UU$lambda<=0] <- min(UU$lambda[UU$lambda>0])
    UU$rank <- 1:L
    if(getFreq)
        UU$freq <- apply(UU$U,2,.findFreq)
    UU$L <- L
    UU$N <- N
    UU$numMisssing <- sum(is.na(x))
    UU$toeplitz=toeplitz
    UU$tau=tau
    UU$seriesName <- deparse(substitute(x))
    UU$call <- match.call()
    ### changed 22.08.2008
    ## class(UU) <- "decompSSAM"
    class(UU) <- c("decompSSAM","decompSSA")
    return(UU)
}

.groupSSAM <- function(x,Us,method){
### Us : selection of eigenvectors
### x : the input series
### method : specifying the method to estimate incompleate lagged vectors
###
### SSA - grouping, see Golyandina et al
    if(!is.matrix(Us) & class(Us)=="numeric")
        Us <- matrix(Us)
    N <- length(x)
    L <- nrow(Us)
    K <- N-L+1
    kk<-1:ncol(Us)
    XX <- matrix(0,nrow=L,ncol=K)
    XX <- matrix(x[row(XX)+col(XX)-1],nrow=L,ncol=K)
    XR <- matrix(NA,nrow=L,ncol=K)
    CC <- (1:K)[!is.na(colSums(XX))]
    nCC <- (1:K)[is.na(colSums(XX))]
    ## STEP 3a: projection of the complete lagged vectors
    for(ii in seq(along.with=CC)){
        expr <- paste("(XX[,CC[ii]]%*%Us[,",kk,"])%*%Us[,",kk,"]",sep="",collapse="+")
        XR[,CC[ii]] <- eval(parse(text=expr))
    }
    ## STEP 3b: Projection of the incomplete lagged vectors
    ## substep ALPHA: nonmissing entries
    if(is.list(method)){
        if(!all(names(method)==c("alpha","beta")))
            stop("If method is a list it must contain alpha and beta")
        if(method$alpha=="PI"){ ## Use the PI - projector
            for(ii in seq(along.with=nCC)){
                ip <- !is.na(XX[,nCC[ii]])
                np <- sum(!ip)
                V <-  Us[ip,]
                W <- Us[!ip,]
                if(np==1){ ## necessary as a vector is interpreted as a column matrix
                    V <- matrix(V,ncol=length(kk))
                    W <- matrix(W,ncol=length(kk))
                }
                PI <- V%*%t(V) + V%*%t(W) %*% solve(diag(1,np,np) - W%*%t(W)) %*% W%*%t(V)
                XR[ip,nCC[ii]] <- PI %*% XX[ip,nCC[ii]]
            }
            if(length(nCC)!=0) ## there may be no columns without missings
                rm(V,W,ip,PI) ## to avoid memory overflow. V, W or PI might be very large!
        } else stop("Other method for alpha than PI is currently not implemented.")
        ## substep BETA: missing entries
        if(method$beta=="simultaneous"){
            for(ii in seq(along.with=nCC)){
                p <- is.na(XX[,nCC[ii]])
                np <- sum(p)
                Up <- Us[p,]
                Uip <- Us[!p,]
                if(np==1){ ## necessary as a vector is interpreted as a column matrix
                    Up <- matrix(Up,ncol=length(kk))
                    Uip <- matrix(Uip,ncol=length(kk))
                }
                XR[p,nCC[ii]] <- solve(diag(1,np,np) - Up %*% t(Up)) %*% Up%*%t(Uip) %*% XR[!p,nCC[ii]]
            }
        } else stop("Other method for alpha than simultaneous is currently not implemented.")
    } else if(method == "PC"){
        ## Projection by means of principal components
        for(ii in 1:ncol(XR)){
            expr <- paste(".star(XX[,ii],Us[,",kk,"]) %*% Us[,",kk,"]",collapse="+")
            XR[,ii] <- eval(parse(text=expr))
        }
    }
    return(XR)
}


reconSSAM <- function(dSSAM,x,groups,method=list(alpha="PI",beta="simultaneous")){
### SSA - Grouping and Diagonal averaging
  U <- dSSAM$U;rm(dSSAM)
  N <- length(x)
  L <- nrow(U)
  K <- N-L+1
  RC <- matrix(0,nrow=N,ncol=length(groups))
  IND <- matrix(nrow=L,ncol=K)
  IND <- row(IND)+col(IND)-1
  for(ii in seq(along=groups)){
    XX <- .groupSSAM(x,U[,groups[[ii]]],method=method)
    # Diagonal Averaging
    .intFun <- function(i,x,ind)
      mean(x[ind==i])
    RC[,ii] <- sapply(1:N,.intFun,x=XX,ind=IND)
  }
  colnames(RC) <- paste("RC",1:length(groups),sep="")
  return(RC)
}
