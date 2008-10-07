\name{kgSSAM}
\alias{kgSSAM}
\alias{CV.kgSSAM}
\title{Itterative filling in of missing values in time series, based on SSA}
\description{
Implementation of the iterative SSA algorithm for time series with missing values introduced by Kondrashov \& Ghil (2006).
}
\usage{
kgSSAM(x, L, toeplitz = FALSE, K, max.inner = 10,
       eps = 0.98, getFreq = TRUE)
CV.kgSSAM(x, L, K_min = 1, K_max = 10, cv.gap = 0.05, runs = 30, ...)

}
\arguments{
  \item{x}{A vector representing the time series.}
  \item{L}{ Embedding dimension or window length. }
  \item{toeplitz}{ Use toeplitz variant of SSA for stationary time series.}
  \item{getFreq}{ Whether dominant frequencies of the eigenelements shall
      be determined. }
  \item{K}{ Vector with index number (rank) of the eigenvectors to be used for 
  		the reconstruction. If it is a scalar the first \code{K} eigenvectors will be used.}
  \item{max.inner}{ Maximum number of iterations for the estimation of each 
  		component represented by \code{K}. This is only used the error of reconstruction
  		is not converging. }
  \item{eps}{Convergence criterion for the iteration which is estimating the 
  		reconstructed components specified by \code{K}. Currently this is implemented 
  		as the minimal coefficient of determination (\eqn{R^2}) value
      between the corresponding reconstructed components between two iteration steps.}
  \item{K_min}{The minimal number of reconstructed components to be used within 
  		cross validation.}
  \item{K_max}{ The maximal number of reconstructed components to be used within 
  		cross validation.}
  \item{cv.gap}{ Proportion of artificially introduced gaps for the cross validation 
  		procedure. }
  \item{runs}{ Number of runs for the cross validation for each parameter combination.}
  \item{\dots}{ arguments to be passed to \code{kgSSAM}.}
}
\value{
	The output of \code{kgSSAM} is an object of class \emph{kgSSAM}
  containing following items:
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
	\item{R2}{The value of the coefficient of determination (\eqn{R^2}{R^2})
            convergence criterion of the last iteration.}
	\item{seriesName}{Name of input series.}
  \item{call }{Call of the generating function}

  The output of \code{CV.kgSSAM} is a \code{runs} times \code{(K_max - K_min)}
  matrix containing the root mean square errors of the reconstructed
  components at the position of artificially introduced gaps for each cross
  validation run. Columns represent the runs.
  Rows represent the additive components used for the reconstruction.
}
\references{ 
Kondrashov, D. & Ghil, M.
Spatio-temporal filling of missing points in geophysical data sets.
Nonlinear Processes in Geophysics, 2006, 13, 151-159
}
\author{ Lukas Gudmundsson }
\section{Warning }{Extreme computational demands.}
  
\seealso{ \code{\link{decompSSAM}} and \code{\link{decompSSA}}}
\keyword{ ts }

