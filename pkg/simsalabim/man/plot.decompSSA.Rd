\name{plot.decompSSA}
\alias{plot.decompSSA}
\alias{plot.decompSSAM}

\title{ Plot the SSA eigenspectrum }
\description{
  Plots the SSA eigenspectrum.
}
\usage{
\method{plot}{decompSSA}(x, by = "freq", normalize = FALSE, asFreq = TRUE,
      ann = TRUE, log = "xy", ...)
%\method{plot}{decompSSAM}(x, by = "freq", normalize = FALSE, asFreq = TRUE,
%      ann = TRUE, log = "xy", ...)
}
\arguments{
  \item{x}{ Object of class \code{decompSSA} }
  \item{by}{How eigenvalues shall be ordered. \code{"rank"} denotes ordering
    by decreasing order.\code{"freq"} denotes ordering by their frequency.}
  \item{normalize}{ If eigenvalues shall be normalized to their sum.}
  \item{asFreq}{Only is applied if \code{by="freq"}. Determines whether the
    frequencies shall be expressed as true \dQuote{frequency} or be inverted as
    \dQuote{period}.}
  \item{ann}{ Whether to do annotation or not.}
  \item{log}{  Use logarithmic axes, usage like \code{log} in
    \code{\link{par}}.}
  \item{\dots}{ Other arguments to be passed to graphical functions.}
}
\details{
The SSA eigenspectrum is build of the (potentially normalised) eigenvalues ordered
either by the dominant frequency of the corresponding eigenvector or by their
rank (decreasingly order).
}
\author{ Lukas Gudmundsson}
\seealso{ \code{\link{decompSSA}} , \code{\link{decompSSAM} }}
\examples{
x <- sin(seq(0,10*pi,len=200))
x <- x + rnorm(x)/2
x.dc <- decompSSA(x,L=40)
plot(x.dc,by="rank",log="")
plot(x.dc,by="freq",log="y")
}

\keyword{ hplot }