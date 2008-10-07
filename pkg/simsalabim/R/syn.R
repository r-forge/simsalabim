################################################################################
## 20.08.2007
## Lukas Gudmundsson
##
## Determination of the phase synchronization of two time series
################################################################################


.hilbert <- function(x){
### Computes the dicrete hilbert transformation of x
  n <- length(x)
  xh <- fft(x)
  hh <- rep(0,n)
  if(2*floor(n/2)==n){ # eaven
    hh[c(1,n/2+1)] <- 1
    hh[2:(n/2)] <- 2
  } else { # odd
    hh[1] <- 1
    hh[2:((n+1)/2)] <- 2
  }
  xh <- fft(xh*hh,inverse = TRUE)/n
  return(xh)
}

.extractPhase <- function(z){
### Extract phase of the complex input vector z
  phase <- rep(NA,length(z))
  phase[Im(z)>=0] <- Arg(z[Im(z)>=0]) - pi
  phase[Im(z)<0] <- Arg(z[Im(z)<0]) + pi
  return(phase)
}

.cumulatePhase <- function(p){
### Cumulates phases
  n <- length(p)
  ii <- sign(p)
  ii <- ii[1:(n-1)]==1 & ii[2:n]==-1
  ii <- c(FALSE,ii)
  ii[ii==TRUE] <- 2*pi
  ii <- cumsum(ii)
  return(p+ii)
}


phaSyn <- function(x1,x2,method="MRL",M=100){
### Computes phase Synchronisation betwean two series x1 & x2
    z1 <- .hilbert(x1)
    z2 <- .hilbert(x2)
    phase1 <- .extractPhase(z1)
    phase2 <- .extractPhase(z2)
    phase1 <- .cumulatePhase(phase1)
    phase2 <- .cumulatePhase(phase2)
    phi <- (phase1-phase2)%%(2*pi) - pi
    phi[phi > 0] <- -phi[phi > 0]%%pi
    phi[phi == 0] <- pi
    phi[phi < 0] <- -(phi[phi < 0]%%pi)
    if(method=="SH"){
        bins <- seq(from=-pi,to=pi,length.out=M+1) # results in M classes in hist
        hh <- hist(phi,breaks=bins,plot=FALSE)
        p <- hh$counts
        hh <- hh$mids
        p[p!=0] <- p[p!=0]/length(phi)
        S_max <- log(M)
        S <- - sum(p[p!=0]*log(p[p!=0]))
        rho <- (S_max-S)/S_max
    } else if(method=="MRL"){
        rho <- sqrt(mean(cos(phi))^2 + mean(sin(phi))^2)
    } else stop("wrong value for method")
    scrn <- list(
                 rho=rho,
                 phi=phi,
                 call=match.call(),
                 method=method,
                 name_x1=deparse(substitute(x1)),
                 name_x2=deparse(substitute(x2))
                 )
    class(scrn) <- "pSyn"
    return(scrn)
}

phaSynMat <- function(X,method="MRL",M=100,verbose=TRUE){
### Phase synchronisation Matrix
### compute pairwise phase synchronisation of the columns of X.
    RHO <- matrix(0,ncol(X),ncol(X))    # index
    RHO[col(RHO)==row(RHO)] <- 1
    dimnames(RHO) <- list(colnames(X),colnames(X))
    if(any(is.na(X)))
        warning("Some coluns of X contain NA \n phase Synchronisation is set to NA.")
    if(verbose)
        cat(date(),": compute synchronisation:\n")
    for(ii in 1:(ncol(X)-1)){
        if((ii/15==round(ii/15)|ii==1) & verbose)
            cat("\t",date(),"Outer loop :",ii,"of",ncol(X),"\n")
        for(jj in (ii+1):ncol(X)){
            if(!any(is.na(X[,ii]),is.na(X[,jj]))){
                RHO[ii,jj] <- phaSyn(X[,ii],X[,jj],method=method,M=M)$rho
            }else{
                RHO[ii,jj] <- NA
            }
        }
    }
    RHO[col(RHO)<row(RHO)] <- t(RHO)[col(RHO)<row(RHO)]
    return(RHO)
}
