\name{plot.MCSSA}
\alias{plot.MCSSA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Monte Carlo SSA }
\description{
  Plots SSA - eigenspectrum and draws monte carlo error bars for each
  eigenvalue.
}
\usage{
\method{plot}{MCSSA}(x, by = "freq", normalize = FALSE, asFreq = TRUE,
      lam.pch = 19, lam.col = "black", lam.cex = 1, sig.col = "red",
      sig.pch = 19, sig.cex = 1, conf.col = "darkgray", log = "xy",
      ann = TRUE, legend = TRUE, axes = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ output of \code{\link{MCSSA}}.}
  \item{by}{ How to order the eigenvalues. \code{"freq"}: by frequency.
             \code{"rank"}: by their rank (decreasingly ordered). }
  \item{normalize}{ Normalise eigenvalues by their sum.}
  \item{asFreq}{ If \code{by="freq"}, shall the frequency be expressed as
                 \dQuote{frequency} or as \dQuote{period}.}
  \item{lam.pch}{ \code{pch} of the eigenvalues. }
  \item{lam.col}{ \code{col} of the eigenvalues. }
  \item{lam.cex}{ \code{cex} of the eigenvalues. }
  \item{sig.col}{ \code{col} of the significant eigenvalues. }
  \item{sig.pch}{ \code{pch} of the significant eigenvalues. }
  \item{sig.cex}{ \code{cex} of the significant eigenvalues. }
  \item{conf.col}{\code{col} of the monte carlo confidence limits.}
  \item{log}{Use logarithmic axes, usage like \code{log} in \code{\link{par}}.}
  \item{ann}{ Whether to do annotations or not.}
  \item{legend}{ Whether to produce a legend or not. }
  \item{axes}{ Whether to draw axes or not. }
  \item{\dots}{ Further arguments to be passed to graphical functions. }
}
\details{
  Eigenvalues larger than the upper bound of the monte carlo confidence limits
  are considered to be not explainable by an AR1 noise process.
}
\seealso{ \code{\link{MCSSA}}, \code{\link{par}} }
\author{ Lukas Gudmundsson }
\keyword{ hplot }