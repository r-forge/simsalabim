################################################################################
## Lukas Gudmundsson
## 30.08.2008
##
##  Simple bootstrapping test for phase synchronizations
#
# 10.09.2007 15:30:31: Single Side added
################################################################################

synBoot <- function(x1,x2,singleSide=FALSE,method="MRL",M=100,R=100){
    .dSide <- function(x1,x2,method,M){
        rr <- phaSyn(sample(x1),sample(x2),method=method,M=M)$rho
    }
    .sSide <- function(x1,x2,method,M){
        rr <- phaSyn(sample(x1),x2,method=method,M)$rho
    }
    bRho <- if(singleSide){
        replicate(R,.dSide(x1,x2,method=method,M=M))
    } else {
        replicate(R,.sSide(x1,x2,method=method,M=M))
    }
    return(bRho)
}

