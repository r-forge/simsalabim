################################################################################
## 19.04. 2007
## Lukas Gudmundsson
##
## Implementation ofthe SSA gapfilling procedure presented by Kondrashov and
## Ghil (Nonlin. Processes Geophys., 13, 151–159, 2006)
################################################################################

kgSSAM <- function(x,L,toeplitz=FALSE,K,max.inner=10,eps=0.98,getFreq = TRUE){
### Kondrashov and Ghil SSA for time series with missings
    if(length(K)==1)
        K <- 1:K
    mpos <- is.na(x)
    x[!mpos]<-x[!mpos]-mean(x,na.rm=TRUE)
    x[mpos]<- 0
    for(k in seq(along=K)){
        for(i in 1:max.inner){
            x.ssa <- decompSSA(x=x,L=L,toeplitz=toeplitz,getFreq=FALSE)
            x.RCsum <- reconSSA(dSSA=x.ssa,x=x,groups=list(K[1:k]))
            if(i==1){
                R2 <- cor(x.RCsum[mpos],x[mpos])^2
                x[mpos]<-x.RCsum[mpos]
                next
            } else {
                R2_old <- R2
                R2 <- cor(x.RCsum[mpos],x[mpos],method = "pearson")^2
                x[mpos]<-x.RCsum[mpos]
            }
            if(R2>eps)
                break
        }
    }
    x.ssa <- decompSSA(x=x,L=L,toeplitz=toeplitz,getFreq=getFreq)
    x[mpos]<-NA
    x.ssa$x <- x
    x.ssa$RC <- x.RCsum
    x.ssa$R2 <- R2
    x.ssa$seriesName <- deparse(substitute(x))
    x.ssa$call <- match.call()
    class(x.ssa) <- "kgSSAM"
    return(x.ssa)
}

CV.kgSSAM <- function(x,L,K_min=1,K_max=10,cv.gap=0.05,runs=30,...){
### Parameter Estimation based on cross validation
    K <- seq.int(K_min,K_max)
    rms <- matrix(NA,length(K),runs)
    nmpos <- !is.na(x)
    nmpos <- (1:length(x))[nmpos]
    cvn <- round(cv.gap*length(nmpos))
    for(r in 1:runs){
        cat(date(),"loop",r,"of",runs,"\n")
        for(k in seq_along(K)){
            sam <- sample(nmpos,cvn)
            xx <- x
            xx[sam] <- NA
            x.ssa <- kgSSAM(xx,L=L,K=K[k],...)
            rms[k,r] <- mean((x.ssa$RC[sam]-x[sam])^2,na.rm=TRUE)
        }
    }
    colnames(rms) <- paste("run",1:runs,sep="")
    rownames(rms) <- paste("RC1_",K,sep="")
    return(rms)
}



