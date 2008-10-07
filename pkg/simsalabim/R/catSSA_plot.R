################################################################################
## 18.08.2007
## Lukas Gudmundsson
##
## Plotting utilities for SSA
################################################################################

plot.decompSSA <-   function(x,by="freq",normalize=FALSE,asFreq=TRUE,ann=TRUE,log="xy",...){
    if(is.null(x$freq))
        x$freq <-  apply(x$U,2,.findFreq)
    if(by == "freq"){
        if(asFreq){
            xx <- x$freq
            xLab <- "frequency"
        } else {
            xx <- 1/x$freq
            xLab <- "period"
        }
    } else if(by == "rank"){
        xx <- x$rank
        xLab <- "rank"
    } else {
        stop("not a valid value for argumet by!")
    }

    if(normalize){
        yLab <- "fraction of explained variance"
        yy <- x$lambda/sum(x$lambda)
    } else {
        yLab <- "Eigenvalues"
        yy <- x$lambda
    }
    plot(xx,yy,ann=FALSE,log=log,...)
    if(ann)
        title(xlab=xLab,ylab=yLab,main=paste("SSA Spectrum of",x$seriesName))
}

### deactivated: 22.08.2008, "decompSSAM" inherits now "decompSSA"
## plot.decompSSAM <-   function(x,by="freq",normalize=FALSE,
##                               asFreq=TRUE,ann=TRUE,log="xy",...){
##     plot.decompSSA(x,by=by,normalize=normalize,asFreq=asFreq,ann=ann,log=log,...)
## }

plot.sdTest <-  function(x,asFreq=TRUE,normalize=TRUE,lam.pch=19,
                         lam.col="black",lam.cex=1,sig.pch=19,sig.col="red",
                         sig.cex=1,conf.col="lightgray",conf.border="lightgray",
                         noise.col="black",noise.lwd=1,ann=TRUE,
                         legend=TRUE,axes=TRUE,...){
    FoP <- x$freq
    oo <- order(FoP)
    FoP <- FoP[oo]
    if(!asFreq)
        FoP <- 1/FoP
    norm <- 1
    if(normalize)
        norm <- sum(x$lambda)
    upper <- x$upper[oo]/norm
    lower <- x$lower[oo]/norm
    nfit <- x$fittedNoise[oo]/norm
    EV <- x$lambda[oo]/norm
    xRange <- range(FoP)
    yRange <- range(EV[EV>0],upper[upper>0],lower[lower>0])
    plot.new()
    plot.window(xRange,yRange,log="xy")
    xx<-c(FoP,FoP[length(FoP):1])
    yy <- c(lower,upper[length(upper):1])
    polygon(xx,yy,col=conf.col,border=conf.border,...)
    points(FoP[EV<=upper],EV[EV<=upper],pch=lam.pch,col=lam.col,cex=lam.cex,...)
    points(FoP[EV>upper],EV[EV>upper],pch=sig.pch,col=sig.col,cex=sig.cex,...)
    lines(FoP,nfit,col=noise.col,lwd=noise.lwd,...)
    if(axes){
        axis(1,...)
        axis(2,...)
    }
    box(...)
    if(ann){
        if(!asFreq){
            xLabel <- "period"
        } else xLabel <- "frequency"
        if(normalize){
            yLabel <- "fraction of explained variance"
        } else  yLabel <- "explained variance"

        conf <- x$conf[x$conf==max(x$conf)]-x$conf[x$conf==min(x$conf)]
        title(xlab=xLabel,ylab=yLabel,main=paste("sdTest of",x$seriesName,
                                      "\np <",conf))
    }
    if(legend){
        if(!asFreq){
            legPos <- "bottomright"
        } else legPos <- "bottomleft"
        leg.txt <- c("eigenvalues", "significant eigenvalues","confidence bounds",
                     "noise model")
        leg.col <- c(lam.col,sig.col,conf.border,noise.col)
        leg.pch <-c(lam.pch,sig.pch,22,-1)
        leg.lty <- c(-1,-1,-1,1)
        leg.bg <- c(1,1,conf.col,1)
        legend(legPos,legend=leg.txt,col=leg.col,pch=leg.pch,
               lty=leg.lty,pt.bg=leg.bg,box.lty=0)
    }
}

plot.MCSSA <- function(x,by="freq",normalize=FALSE,asFreq=TRUE,lam.pch=19,
                       lam.col="black",lam.cex=1,sig.col="red",sig.pch=19,sig.cex=1,
                       conf.col="darkgray",
                       log="xy",ann=TRUE,legend=TRUE,axes=TRUE,...){
    if(by == "freq"){
        if(asFreq){
            xx <- x$freq
            xLab <- "frequency"
        } else {
            xx <- 1/x$freq
            xLab <- "period"
        }
    } else if(by == "rank"){
        xx <- x$rank
        xLab <- "rank"
    } else {
        stop("not a valid value for argumet by!")
    }
    if(normalize){
        yLab <- "fraction of explained variance"
        yy <- x$lambda/sum(x$lambda)
        lower <- x$lower/sum(x$lambda)*100
        upper <- x$upper/sum(x$lambda)*100
    } else {
        yLab <- "Eigenvalues"
        yy <- x$lambda
        lower <- x$lower
        upper <- x$upper
    }
    xRange <- range(xx)
    yRange <- range(yy[yy>0],upper[upper>0],lower[lower>0])
    plot.new()
    plot.window(xRange,yRange,log=log)
    arrows(xx,upper,xx,lower,code=3,
           angle=90,length=0.025,col=conf.col,...)
    points(xx[yy<=upper],yy[yy<=upper],pch=lam.pch,col=lam.col,cex=lam.cex,...)
    points(xx[yy>upper],yy[yy>upper],pch=sig.pch,col=sig.col,cex=sig.cex,...)
    if(axes){
        axis(1,...)
        axis(2,...)
    }
    box(...)
    if(ann){
        if(by=="freq"){
            if(asFreq){
                xLabel <- "frequency"
            } else xLabel <- "period"
        }
        if(normalize){
            yLabel <- "fraction of explained variance"
        } else  yLabel <- "explained variance"

        conf <- x$conf[x$conf==max(x$conf)]-x$conf[x$conf==min(x$conf)]
        title(xlab=xLabel,ylab=yLabel,main=paste("sdTest of",x$seriesName,
                                      "\np <",conf))
    }
    if(legend){
        if(!asFreq){
            legPos <- "bottomright"
        } else legPos <- "bottomleft"
        leg.txt <- c("eigenvalues", "significant eigenvalues","confidence bounds")
        leg.col <- c(lam.col,sig.col,conf.col)
        leg.cex <- c(lam.cex,sig.cex,1)
        leg.pch <-c(lam.pch,sig.pch,-1)
        leg.lty <- c(-1,-1,1)
        legend(legPos,legend=leg.txt,col=leg.col,pch=leg.pch,lty=leg.lty,box.lty=0)
    }
}
