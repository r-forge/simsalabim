\name{w.cor}
\alias{w.cor}

\title{ w - correlation}
\description{
  Computes w - correlations (Golyandina, et.al. 2001) from a matrix
  containing reconstructed components in its columns.
}
\usage{
w.cor(R, L)
}
\arguments{
  \item{R}{ A matrix containing reconstructed components in its columns.}
  \item{L}{ The window length, may be omitted in case of complete set of
            reconstructed components.}
}
\details{
  W - correlation can be used to examine whether reconstructed components (RC) are
  independent (uncorrelated) or not. This is useful in signal to noise enhancement.
  RC being correlated to a high degree are assumed to correspond to the same signal.
}
\value{
  A suqare matrix containing the w - correlations.
}
\references{
Golyandina, N.; Nekrutkin, V. & Zhiglkilavskifi, A. Analysis of Time Series Structure: SSA and Related Techniques. CRC Press, 2001
}
\author{ Lukas Gudmundsson}

\seealso{ \code{\link{decompSSA}}, \code{\link{decompSSAM}} }
\examples{
x <- sin(seq(0,10*pi,len=200))
x <- x + rnorm(x)/2
x.dc <- decompSSA(x,L=40)
x.rc <- reconSSA(x.dc,x,as.list(1:40))
x.wc <- w.cor(x.rc,40)
image(x.wc,col=gray(100:0/100))
}

\keyword{ ts }
