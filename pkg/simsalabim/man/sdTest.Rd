\name{sdTest}
\alias{sdTest}
%\alias{sdTest.default}
\alias{sdTest.decompSSA}
\alias{sdTest.decompSSAM}
\title{ Test SSA eigenspectrum against red noise }
\description{
Performs a test of an SSA eigenspectrum against a red noise hypothesis,
following  Shun and Duffy (1999).
}
\usage{
sdTest(dSSA, ...)
\method{sdTest}{decompSSA}(dSSA, conf = 0.95,...)
%\method{sdTest}{decompSSAM}(dSSA, conf = 0.95,...)
}
\arguments{
  \item{dSSA}{An object of class \code{\link{decompSSA}}}
  \item{conf}{Confidence limit}
  \item{\dots}{ Further arguments to be passed to different methods. }

}
\details{
  An implementation of the test introduced by
  Shun and Duffy (1999). The crucial step within the computations
  is a non linear fit by \code{\link{nls}}. In some cases \code{\link{nls}}
  throws an error which causes the test to fail. In such cases the SSA -
  spectrum does not have enough resemblance with the model to allow a
  sensible fit. Therefore, it can be assumed that the red noise null hypothesis
  is not appropriate anyway.

  The confidence level \code{conf} can be specified either as number < 1 or as
  vector of length two with \code{sum(conf)==1}.
}
\value{
  An object of class \emph{sdTest}
  \item{lambda }{The eigenvalues, ordered by decreasing frequency.}
  \item{freq }{Dominant frequency of the eigenvectors, ordered by decreasing.}
  \item{rank }{Rank of the eigenvalues, ordered by decreasing frequency.}
  \item{fittedNoise }{Fitted values of the noise model,
        ordered by decreasing frequency.}
  \item{upper }{Upper confidence limit, ordered by decreasing frequency.}
  \item{lower }{Lower confidence limit, ordered by decreasing frequency.}
  \item{conf }{ Vector of length two specifying the upper and lower confidence
      bounds of the test.}
  \item{N }{Length of the initial series}
  \item{L }{The embedding dimension}
  \item{call }{Call of the generating function}
}
\references{
Shun, T. & Duffy, C. Low-frequency oscillations in precipitation, temperature, and runoff on a west facing mountain front: A hydrogeologic interpretation. Water Resources Research, 1999, 35, 191-201
}
\author{ Lukas Gudmundsson }

\section{Warning }{In some cases the function throws an error due to a failure
of a non linear fitting procedure (see details above).}

\seealso{ \code{\link{decompSSA}}, \code{\link{decompSSAM}},
\code{\link{plot.sdTest}},\code{\link{MCSSA}} }
\examples{
n <- arima.sim(n=200,alist(ar=0.4))
n <- as.vector(scale(n))
x <- sin(seq(0,15*pi,len=200))
xn <- x+n/2
xn.dc <- decompSSA(xn,L=40)
plot(sdt<-sdTest(xn.dc))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest }
\keyword{ ts }
