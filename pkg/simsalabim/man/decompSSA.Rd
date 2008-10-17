\name{decompSSA}
\alias{decompSSA}
\alias{reconSSA}

\title{ Singular Spectrum Analysis }
\description{
  Performs the two stages of Singular System Analysis (SSA). Decomposition
  and reconstruction.
}
\usage{
  decompSSA(x, L, toeplitz = FALSE, getFreq = TRUE)
  reconSSA(dSSA,x, groups)
}
\arguments{
  \item{x}{A vector representing the time series.}
  \item{L}{ Embedding dimension also referred to as window length.}
  \item{toeplitz}{ Whether to use the Toeplitz modification of SSA for
        stationary time series or not.}
  \item{getFreq}{ Whether dominat frequencies of the eigenvectors shall
      be determined. }
  \item{dSSA}{Output of decompSSA. }
  \item{groups}{ A list of vectors. Each vector is representing a
    selection of eigenvalues and eigenvectors which shall be used to
    compute reconstructed components.}
}
\details{
	\code{decompSSA} performs the SSA decomposition whereas \code{reconSSA}
  performs SSA reconstruction.

  Singular Spectrum Analysis (SSA) embeds
  lagged copies of a time series \eqn{x} into a augmented (trajectory)
  matrix \eqn{X}. The orthonormal basis of \eqn{X} is  via  singular
  value decomposition (SVD). Having obtained this decomposition, a
  selection  of eigenvectors (or left singular vectors) can be used to
  filter the  time series. SSA has one free parameter, the window length
  \code{L} which determines how the time series is embedded. A general
  advice is to choose \code{L} close to half of the series length.

  SVD of an matrix \eqn{X} can be split up into several steps. One of them
  beeing the eigenvalue decomposition of \eqn{XX^T}.
  The product \eqn{XX^T } has simmilarities whith
  the covariance matrix of \eqn{X}. Moreover if \eqn{X} is set to have zero
  mean and \eqn{XX^T} is normalised, it is in fact a covariance matrix.
  In case of SSA \eqn{XX^T} containes the inner product of time lagged
  copies of the time series. The Toeplitz modification of SSA puts
  stationarity assumptions into the computation of \eqn{XX^T}.
  In such a case the secondary diagonals of \eqn{XX^T} only depend on the
  time lag. This can be enforced by giving these entries the values of the
  autocorrelation function at the corresponding lag.

  It is possible to assign each eigenvector a dominat frequency. This
  is found by applying \code{\link{fft}} taking the frequency of the largest
  Fourier component.

  Beside the results of the output of \code{decompSSA} and the original
  time series \code{x}, \code{reconSSA} has a third input called
  groups. This is implemented a list of vectors where each vector should
  contain integers representing the index of the eigenvectors  to be
  used for reconstruction of a sub signal of the initial time
  series. The index corresponds to the index of the decreasing ordered
  eigenvalues (\code{rank}).

  For reasons computational efficiency the sum of several components (represented
  by the single list entries) can be obtained in one step. This procedure is
  recommended if the interest lies on the shape of the entire reconstructed
  signal. However, if the focus of the analysis is to examine the shape of
  single reconstructed components one has to obtain each each reconstructed
  compontent seperately.
}
\value{
  The output of \code{decompSSA} is an object of class \emph{decompSSA}
  with following items:
  \item{lambda }{The eigenvalues, ordered by decreasing value.}
  \item{U }{The eigenvectors (columns), ordered by decreasing eigenvalues.}
  \item{freq }{Dominant frequency of the eigenvectors,
    ordered by decreasing eigenvalues.}
  \item{rank }{Rank of the eigenvalues, ordered by decreasing eigenvalues.}
  \item{N }{Length of the input series}
  \item{L }{The embedding dimension}
  \item{toeplitz}{Logical, indicates if Toeplitz modification has been used.}
  \item{seriesName}{Name of input series.}
  \item{call }{Call of the generating function.}

  The output of \code{reconSSA} is a matrix with \code{length(groups)} columns
  and \code{length(x)} rows. Each column represents the sum of
  the reconstructed components defined by the list entries of \code{groups}.

}
\references{
  
  Elsner, J. B. & Tsonis, A. A. Singular Spectrum Analysis: A New Tool in
  Time Series Analysis Springer, 1996 
  
  Ghil, M.; Allen, M.; Dettinger, M.; Ide, K.; Kondrashov, D.; Mann, M.;
  Robertson, A.; Saunders, A.; Tian, Y.; Varadi, F. & others Advanced
  spectral methods for climatic time series. Rev. Geophys, 2002, 40, 1003 
  
  Golyandina, N.; Nekrutkin, V. & Zhiglkilavskifi, A. Analysis of Time
  Series Structure: SSA and Related Techniques. CRC Press, 2001 
}
\author{ Lukas Gudmundsson }
\section{Warning }{
  May cause extreme memory demands. \code{reconSSA} is computational expensive.
}

\seealso{ \code{\link{plot.decompSSA}}, \code{\link{decompSSAM}},
  \code{\link{kgSSAM}}, \code{\link{sdTest}}, \code{\link{MCSSA}},
  \code{\link{svd}}, \code{\link{eigen}}, \code{\link{toeplitz}},
  \code{\link{fft}}
}

\examples{
x <- sin(seq(0,10*pi,len=200))
x <- x + rnorm(x)/2
x.dc <- decompSSA(x,L=40)
plot(x.dc,by="rank",log="") # the first two elements contain the signal
pairs(x.dc$U[,1:4])  # the first two elements contain the signal

# Obtaining reconstructed components(RC)
x.rc1 <- reconSSA(x.dc,x,list(1:2)) # the signal
x.rc2 <- reconSSA(x.dc,x,list(1:40)) # full reconstruction of the series
plot(x,type="l")
lines(x.rc1,col="red",lwd=2)
points(x.rc2,col="blue")


# Obtaining several RC at once
# first colum: signal
# second column: noise
x.rc3 <- reconSSA(x.dc,x,list(c(1,2),3:40))
dim(x.rc3)
matplot(x.rc3,type="l")

# The sum of RC can be obtained in one step.
# This avoids costly computations that would occure if all
# RC would be obtained seperately.
x.rc4 <- reconSSA(x.dc,x,as.list(1:40)) # separate RCs for all eigenvectors
sum(rowSums(x.rc4)-x.rc2) # identical up to numerical accuracy
}

\keyword{ ts }

