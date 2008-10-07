\name{phaSyn}
\alias{phaSyn}
\alias{phaSynMat}
\title{ Phase Synchronisation Analysis}
\description{
  Computes the phase synchronisation index between to time series.
}
\usage{
phaSyn(x1,x2,method="MRL",M=100)
phaSynMat(X,method="MRL",M=100,verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1, x2}{ Numeric vectors of the same length, interpreted as time series.}
  \item{method}{How to determine the phase synchronisation index.
            \code{"MRL"}: mean resultant length. \code{"SH"}: a measure based
            on the Shannon - Entropy.}
  \item{M}{ Number of bins in the histogram of the cyclic phase difference.
            This is only used if \code{method="SH"} and is ignored otherwise.}
  \item{X}{A matrix containing time series in its columns.}
  \item{verbose}{If \code{TRUE} informations about the progress of the
    computation are printed.}
}
\details{
  Phase synchronisation is quantified in two steps: First, the phase of both series is
  extracted by means of the Hilbert -   transformation. Second, the distribuition
  of the cyclic phase difference of both series is checked for uniformity. In
  case of uniformity the series are not synchronous in respect to their phase.

  Uniformity of the cyclic phase difference can be quantified in two ways. If
  \code{method="MRL"}, the dispersion of the barycenter of the phase difference
  on the unit circle is calculated. In the case of \code{method="SH"}, a
  histogram of the phase difference is built and its uniformity quantified by
  means of the Shannon - entropy.

  The indices of both methods are scaled in such a manner that 1 means perfect
  synchronisation and 0 stands for no synchronisation.

}
\value{
  The output of \code{phaSyn} is an object of class \emph{pSyn} containing:
  \item{rho}{ Phase synchronisation index.}
  \item{phi}{ Cyclic phase difference. }
  \item{call}{The call of the generating function.}
  \item{method}{The method used to quantify peakedness of dsitribution.}
  \item{name_x1}{ Name of series \code{x1}.}
  \item{name_x2}{ Name of series \code{x2}.}

  The output of \code{phaSynMat} is a symmetric matrix containing the pairwise
  phase synchronisation index.
}
\references{
Allefeld, C. & Kurths, J. Testing for phase synchronization.
International Journal of Bifurcation and Chaos, 2004, 14, 405-416

Rybski, D.; Havlin, S. & Bunde, A.
Phase synchronization in temperature and precipitation records.
Physica A: Statistical Mechanics and its Applications, Elsevier, 2003, 320, 601-610

}
\author{ Lukas Gudmundsson}

\seealso{ \code{\link{synBoot}}}
\examples{
# phase synchronisation
x1 <- sin(1:100)
x2 <- cos(1:100)
phaSyn(x1,x2)

# phase synchronisation matrix
xx<-matrix(seq(0,6*pi,len=100),ncol=30,nrow=100)
colnames(xx) <- 1:30
xx[,1:10] <- sin(xx[,1:10])
xx[,11:20] <- sin(xx[,11:20]*2*pi)
xx[,21:30] <- sin(xx[,21:30]*4*pi)
cxx <- phaSynMat(xx,verbose=FALSE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
