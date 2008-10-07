\name{synBoot}
\alias{synBoot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Permutation test for phase synchronisation indices }
\description{
  Computes a permutation test on phase synchronisation indices.
}
\usage{
synBoot(x1,x2,singleSide=FALSE,method="MRL",M=100,R=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1, x2}{ Numeric vectors of the same length, interpreted as time series.}
  \item{singleSide}{ If true only \code{x1} is permutated.}
  \item{method}{How to determine the phase synchronisation index.
    \code{"MRL"}: mean resultant length. \code{"SH"}: a measure based on the
    Shannon - Entropy.}
  \item{M}{ Number of bins in the histogram of the cyclic phase difference,
            needed to compute the Shannon - entropie based phase
            synchronisation index. This is only used if \code{method="SH"}
            and is ignored otherwise.}
  \item{R}{ Number of replications for random  permutation. }
}
\details{
  Permutates the values of \code{x1} and \code{x2} \code{R} times and computes
  for every replication the phase synchronisation index.
  The resulting vector of surrogate phase synchronisation
  indices may be used to compute quantiles, that can be used for significance
  testing (see example below).
}
\value{
  A numeric vector of length \code{R} containing the Phase Synchronisation
  indices between both time series for each permutation.
}
\references{
Allefeld, C. & Kurths, J. Testing for phase synchronization.
International Journal of Bifurcation and Chaos, 2004, 14, 405-416
}
\author{ Lukas Gudmundsson}
\section{Warning }{May be computational expernsive}

\seealso{ \code{\link{phaSyn}}}
\examples{

x1 <- sin(1:100)
x2 <- cos(1:100)
ps1 <- phaSyn(x1,x2)$rho
psp <- synBoot(x1,x2)
quantile(psp,0.95) < ps1 # Significance at the 95\% confidence level?

}

\keyword{ ts }
\keyword{ htest }
