\name{kgSSAM}
\alias{kgSSAM}
\alias{cvSSA.kgSSAM}
\title{Itterative filling in of missing values in time series, based on SSA}
\description{
Implementation of the iterative SSA algorithm for time series with missing values introduced by Kondrashov \& Ghil (2006).
}
\usage{
kgSSAM(x, L, toeplitz=FALSE, K, max.inner=100, eps=0.05,
       getFreq=TRUE, verbose=FALSE)

cvSSA.kgSSAM(x,L,K=10,cv.crit="bestApprox",cv.gap=0.05,
             runs=200,verboseCV=TRUE,...)
}
\arguments{
  \item{x}{A vector representing the time series.}
  \item{L}{ Embedding dimension or window length. }
  \item{toeplitz}{ Use toeplitz variant of SSA for stationary time series.}
  \item{getFreq}{ Whether dominant frequencies of the eigenelements shall
      be determined. }
  \item{K}{ Vector with index number (rank) of the eigenvectors to be
             used the reconstruction. If it is a scalar the first
             \code{K}
             eigenvectors will be used.}
  \item{max.inner}{ Maximum number of iterations for the estimation of
             each component represented by \code{K}. This is only used
             the error of reconstruction is not converging. }
  \item{eps}{Convergence criterion for the iteration which is estimating
             the reconstructed components specified by
             \code{K}. Currently this is implemented as the minimal
             normalized root mean square value between the
             corresponding reconstructed components between two
             iteration steps.}
  \item{cv.gap}{ Proportion of artificially introduced gaps for the
    cross validation procedure.}
  \item{cv.crit}{Optimization criterion for cross
             validation. \code{"bestApprox"}} for the lowest
             root-mean-square error (rms) between the cumulative reconstructed
             components and the original time
             series. \code{"stabSignal"} The rms between the cumulative
             reconstructed (RC) of the original series components and
             the cumulative RC of the original time series.
  \item{runs}{ Number of runs for the cross validation for each
    parameter combination.}
   \item{verbose}{Print information about progress}
   \item{verboseCV}{Print information about progress cross validation
     progress}
  \item{\dots}{ arguments to be passed to \code{kgSSAM}.}
}

\details{
  For time series in many applications, a few leading SSA modes
  correspond to the record's dominant oscillatory and/or trend  modes
  ("smooth" signal), while the rest is noise. SSA gap-filling algorithm
  estimates such "smooth" component from the incomplete times series and
  uses it to fill-in the gaps.  The SSA gap-filling algorithm consists
  of two iteration loops. The inner-loop iteration is started by
  computing the leading SSA eigenvector E1 of the centered time series
  with the unbiased value of the mean, and zero-padded in place of the
  missing points. Then the SSA algorithm is performed again on the new
  time series, in which the RC R1 corresponding to the E1 mode alone was
  used to obtain nonzero values in place of the missing points and
  correct the record's mean, the covariance matrix and the mode E1
  itself. The reconstruction of the missing data is repeated with a new
  estimate of R1 and tested against the previous one, until a
  convergence test has been satisfied. Next, outer loop of iterations is
  performed by adding a second SSA eigenvector E2 for reconstruction,
  starting from the solution with data filled-in by R1.

  The optimum number of SSA modes for gap-filling are found from a set
  of cross-validation experiments; for each such experiment, a fixed
  fraction of available data is flagged as missing, and the
  root-mean-square (rms) error in reconstruction is computed as a
  function of the number  of SSA modes retained. The global minimum in
  error, averaged over all experiments, corresponds to the required
  optimum, and provides an estimate of the actual error in the
  reconstructed dataset. Alternatively the cross-validation criterion
  can be set to asses the similarity of the cumulative reconstructed
  components of the original case and the case with artificially
  introduced missing values. Optimum SSA window size can be found in a
  similar way as well.
}


\value{
  The output of \code{kgSSAM} is an object of class
  \emph{kgSSAM}. The output of \code{cvSSA.kgSSAM} is an object of
  class \emph{cvSSA.kgSSAM} that inherits \emph{kgSSAM}. They
  containin following items:
  \item{lambda }{The eigenvalues, ordered by decreasing value.}
  \item{U }{The eigenvectors (columns), ordered by decreasing eigenvalues.}
  \item{freq }{Dominant frequency of the eigenvectors,
    ordered by decreasing eigenvalues.}
  \item{rank }{Rank of the eigenvalues, ordered by decreasing eigenvalues.}
  \item{N }{Length of the input series.}
  \item{L }{The embedding dimension.}
  \item{toeplitz}{Logical, indicates if toeplitz variant has been used.}
  \item{x}{The input series.}
  \item{RC}{The reconstructed components. Corresponding to \code{K}}
  \item{nrmse}{Normalized root mean square error for the last iteration
    step.}
  \item{CVerror}{Root mean square cross validation error (\eqn{runs * K})
             matrix.}
  \item{seriesName}{Name of input series.}
  \item{call}{Call of kgSSAM}
  \item{CVcall}{Call of }
}

\references{
Kondrashov, D. & Ghil, M.
Spatio-temporal filling of missing points in geophysical data sets.
Nonlinear Processes in Geophysics, 2006, 13, 151-159
\url{http://www.nonlin-processes-geophys.net/13/151/2006/npg-13-151-2006.pdf}
}
\author{ Lukas Gudmundsson and Dmitri Kondrashov}

\section{Warning }{Extreme computational demands.}



\seealso{ \code{\link{decompSSAM}} and \code{\link{decompSSA}}}

\examples{
## generate "benchmark"
x <- sin(seq(0,10*pi,len=200))
x <- x + rnorm(x)/2
xna<-x
xna[100:120] <- NA
x.dc <- decompSSA(x,L=40)
x.rc <- reconSSA(x.dc,x,list(1:2))

## try out different cross - validation criteria
cv.bApp <- cvSSA.kgSSAM(xna, 40,K = 10, cv.gap = 0.15,
                   runs = 5,max.inner = 200,
                   eps = 0.005,
                   cv.crit="bestApprox")

cv.stSig <- cvSSA.kgSSAM(xna, 40,K = 10, cv.gap = 0.15,
                   runs = 5,max.inner = 200,
                   eps = 0.005,
                   cv.crit="stabSignal")

## plot results
layout(rbind(1:2,c(3,3)))
matplot(cv.bApp$CVerror,t="l",col="gray",lty=1,
        main="Best approximation of series",ylab="rms")
lines(meanCV.bApp <- rowMeans(cv.bApp$CVerror),lwd=2,col="red")
points(meanCV.bApp,pch=19)

matplot(cv.stSig$CVerror,t="l",col="gray",lty=1,
        main="Stability of RC",ylab="rms")
lines(meanCV.stSig <- rowMeans(cv.stSig$CVerror),lwd=2,col="red")
points(meanCV.stSig,pch=19)

plot(x,t="l",col="gray")
lines(xna,col="black")
lines(x.rc,col="red",lwd=3)
lines(cv.bApp$RC[,meanCV.bApp==min(meanCV.bApp)],col="blue",lwd=2)
lines(cv.stSig$RC[,meanCV.stSig==min(meanCV.stSig)],col="green",lwd=2,lty=2)

legend("bottom",legend=c("orig","gap (best approx.)","gap (signal stab.)"),
      fill=c("red","blue","green"),horiz=TRUE,
      bg="white")


}

\keyword{ ts }

