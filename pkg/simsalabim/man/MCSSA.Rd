\name{MCSSA}
\alias{MCSSA}
%\alias{MCSSA.default}
\alias{MCSSA.decompSSA}

\title{ Monte Carlo SSA (MCSSA)}
\description{
Simple implementation of MCSSA
}
\usage{
MCSSA(dSSA, ...)
%\method{MCSSA}{default}(U, lambda, x, n, conf = 0.95, toeplitz = FALSE, keepSurr = FALSE)
\method{MCSSA}{decompSSA}(dSSA, x, n, conf = 0.95,
      keepSurr = FALSE,ar.method="mle", ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dSSA}{Object of class \code{\link{decompSSA}}.}
  \item{x}{ A numeric vector, interpreted as time series.}
  \item{n}{ Number of surrogates to be computed. }
  \item{conf}{Confidence limit.}
  \item{keepSurr}{ Whether to keep the surrogate estimates of lambda or
    not.}
  \item{ar.method}{Method to estimate AR1 parameters, passed on to
    \code{\link{ar}}. One of \code{"yule-walker"}, \code{"burg"},
    \code{"ols"}, \code{"mle"}, \code{"yw"}. }
  \item{\dots}{Arguments passed to methods.}
}
\details{
A simple implementation of the MCSSA, testing against the AR1 hypothesis.
The parameters of the AR1 process are estimated via \code{\link{ar}} .

The confidence level \code{conf} can be specified either as number < 1 or
as vector of length two with \code{sum(conf)==1}.
}
\value{
  An object of class \emph{MCSSA}
  \item{lambda }{The eigenvalues, ordered by decreasing frequency.}
  \item{freq }{Dominant frequency of the eigenvectors, ordered by
    decreasing value.}
  \item{rank }{Rank of the eigenvalues, ordered by decreasing
    frequency.}
  \item{ar1 }{The noise model, output of \code{\link{ar.mle}}}
  \item{upper }{Upper confidence limit, ordered by decreasing
    frequency.}
  \item{lower }{Lower confidence limit, ordered by decreasing
    frequency.}
  \item{conf }{ Vector of length two, specifying the upper and lower
    confidence bounds of the test.}
  \item{N }{Length of the initial series.}
  \item{L }{The embedding dimension.}
  \item{call }{ Call of the generating function.}
  \item{surLambda}{A matrix containing the estimates of the surrogate
      eigenvalues, ordered by the dominant frequency in its columns.
      Only available if \code{keepSurr=TRUE}.}

}
\references{
Allen, M. & Smith, L. Monte Carlo SSA: Detecting irregular
oscillations in the Presence of Colored Noise. Journal of Climate,
1996, 9, 3373-3404

}
\author{ Lukas Gudmundsson}
\section{Warning }{Large computational demands}

\seealso{ \code{\link{sdTest}}, \code{\link{decompSSA}},
  \code{\link{plot.MCSSA}}
}
\examples{
n <- arima.sim(n=200,alist(ar=0.4))
n <- as.vector(scale(n))
x <- sin(seq(0,15*pi,len=200))
xn <- x+n/2
xn.dc <- decompSSA(xn,L=40)
xn.mct <- MCSSA(xn.dc,xn,n=50)
plot(xn.mct)
}

\keyword{ ts }
\keyword{ htest }
