\name{decompSSAM}
\alias{decompSSAM}
\alias{reconSSAM}

\title{ Singular System Analysis for time series with missing values}
\description{
  Implementation of the SSA algorithm for time series with missing values
  following Golyandina and Osipov (2007).
}
\usage{
decompSSAM(x, L, tau = 0, toeplitz = FALSE, getFreq = TRUE)
reconSSAM(dSSAM, x, groups,
      method = list(alpha = "PI", beta = "simultaneous"))

}
\arguments{
  \item{x}{A vector representing the time series.}
  \item{L}{ Embedding dimension.}
  \item{tau}{ A number ranging from 0 to 1, indicating the maximal proportion
    of missing values within columns of the trajectory matrix for the usage
    of a modification of the inner product. If tau is negative it is treated
    as if it is zero. Only applies if \code{toeplitz = FALSE}.}
  \item{toeplitz}{ Whether to use the Toeplitz modification of SSA for
        stationary time series or not.}
  \item{getFreq}{ Whether dominant frequencies of the eigenvectors shall
      be determined.}
  \item{dSSAM}{Output of decompSSAM }
  \item{groups}{ A list of vectors. Each vector is representing a
    selection of eigenvalues and eigenvectors which shall be used to
    compute reconstructed components.}
  \item{method}{Method or combination of methods that shall be used to estimate
    the reconstructed components at the position of missing values. Currently
    one of \code{"PC"} or \cr
    \code{list(alpha="PI",beta="simultaneous")}.}
}
\details{
  \code{decompSSAM} performs the SSA decomposition whereas \code{reconSSAM}
  performs the SSA reconstruction, which is a type of band pass filtering.

  In general the application of SSA for time series with missing values follows
  the same principals as for standard SSA. For general comments on the
  application of these functions see the documentation of
  \code{\link{decompSSA}}.

  SSA embeds lagged copies of a time series \eqn{x} into a augmented matrix
  \eqn{X} (trajectory matrix ). In a second step the orthonormal basis
  of \eqn{X} is found via singular value decomposition (SVD). One of the
  internal steps of SVD is the computation of \eqn{XX^T}. There are two
  strategies available to obtain this product if values of \eqn{x} are
  missing. The first possibility is to omit any column of \eqn{X} and to
  compute \eqn{XX^T} for that reduced matrix. The other possibility is
  to use a modification of the inner product for vectors containing
  missing values. The threshold parameter \code{tau} controls the
  computation of \eqn{XX^T}. It gives the acceptable proportion of
  missing values whithin a vector for the applicatoin of the modified
  inner product. If \code{tau} is \eqn{<=0}, all columns of \eqn{X}
  containing any missing value are ommited. If \code{tau} \eqn{>=1} no
  column of \code{X} is ommited. For any value of \code{tau} in between
  all columns of \code{X} having a propotion of more or equal than
  \code{tau} are omitted. There are several possibilities to capture the
  reconstructed components. One is the recovery by \emph{means of
  principal component}, which has first been introduced by Schoellhamer
  (2001). Golyandina and Osipov split the recovery of the reconstructed
  signals into two steps \code{alpha} and \code{beta} (see reference for
  more details). Currently only one method for each step is implemented,
  namely  \code{alpha="PI"}, denoting the usage of the \emph{\eqn{Pi} -
  projector} and \code{beta="simultaneous"}, denoting \emph{simultaneous
  filling in}.
}
\value{
  The output of \code{decompSSAM} is an object of class
  \emph{decompSSAM} inheriting from \emph{decompSSA} with following
  items:
  \item{lambda }{The eigenvalues, ordered decreasing.}
  \item{U }{The eigenvectors (columns), ordered by decreasing eigenvalues.}
  \item{freq }{Dominant frequency of the eigenvectors,
    ordered by decreasing eigenvalues.}
  \item{rank }{Rank of the eigenvalues, ordered by decreasing eigenvalues.}
  \item{N }{Length of the input series.}
  \item{L }{Embedding dimension or window length.}
  \item{toeplitz}{Logical, indicates if Toeplitz modification has been used.}
  \item{numMisssing}{The number of missing entries in the input.}
  \item{tau}{The user specification of \code{tau}.}
  \item{seriesName}{Name of input series.}
  \item{call }{Call of generating function.}

  The output of \code{reconSSAM} is a matrix with \code{length(groups)} column
  and \code{length(x)} rows. Each columns represents the sum of the
  reconstructed components defined by the list entries of groups.

}
\references{
Golyandina, N. & Osipov, E.
The "Caterpillar"-SSA method for analysis of time series with missing values.
Journal of Statistical Planning and Inference, 5th St.
Petersburg Workshop on Simulation, 2007, 137, 2642-2653
\url{http://www.gistatgroup.com/cat/mvssa1en.pdf}

Schoellhamer, D.
Singular spectrum analysis for time series with missing data.
Geophysical Research Letters, 2001, 28, 3187-3190
}
\author{Lukas Gudmundsson}
\section{Warning }{
  May cause extreme memory demands. \code{reconSSAM} is computionally expensive.
}
\seealso{
  \code{\link{plot.decompSSAM}}, \code{\link{decompSSA}}, \code{\link{sdTest}}
}
\examples{
x <- sin(seq(0,10*pi,len=200))
x <- x + rnorm(x)/2
x[100:120] <- NA
x.dc <- decompSSAM(x,L=40)
plot(x.dc,by="rank",log="")
x.rc1 <- reconSSAM(x.dc, x, groups=list(1:2),
  method=list(alpha = "PI", beta = "simultaneous"))
x.rc2 <- reconSSAM(x.dc, x, groups=list(1:2),
  method="PC")
# compare result with input signal
plot(x,type="l")
lines(x.rc1,col="red",lwd=2,lty=1)
lines(x.rc2,col="black",lwd=2, lty=2)
lines(sin(seq(0,10*pi,len=200)),col="blue",lwd=2)
}
\keyword{ ts }
