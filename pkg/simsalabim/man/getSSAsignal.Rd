\name{getSignal}
\alias{getSignal}
\alias{getSignal.decompSSA}
\alias{getTrend}
\alias{getPeriod}
\title{Parametric method for identifying trends and periodic components
  in SSA decompositions. }
\description{
Signal detection algorithm based on ideas of Vautard et.al. (1992) and
Alexandrov (2008a,b)
}
\usage{
getSignal(dSSA, ...)
\method{getSignal}{decompSSA}(dSSA, omega0=1/dSSA$L,C0,r0,...)
getTrend(U,omega0=1/nrow(U),C0)
getPeriod(U,r0)

}
\arguments{
  \item{dSSA}{An object of class \code{\link{decompSSA}}.}
  \item{omega0}{Low frequency boundary.}
  \item{C0}{Low frequency contribution.}
  \item{r0}{Dominant frequency contribution.}
  \item{U}{Matrix of SSA -  Eigenvectors.}
  \item{\dots}{ Further arguments to be passed to different methods. }

}
\details{
  This signal selection algorithm is based on the fact that the SSA -
  Eigenvectors \eqn{U^{(i)}}{U^i} for trends and periodic components have a simple
  shape and their periodogram
  %%
  \deqn{
    \Pi^{(i)}(\omega)=
    \frac{1}{L}\left|\sum_{n=0}^{L-1}\exp\{-i2\pi \omega n\}
    U_n^{(i)}\right|^2
  }{
    Pi^i(omega) = 1/L |sum_n{exp(-i 2 pi omega n)}
  }
  %%
  is governed by a dominant frequency. (The periodogram is determined by
  means of \code{\link{fft}} and therefore the spacing of \eqn{\omega}{omega} is
  the Nyquist frequency \eqn{1/L}.)

  The \strong{\emph{trend}} extraction is based on the idea that the $U^{(i)}$
  contributing to a trend have large periodogram values only for low
  frequencies \eqn{\leq \omega_0}{<= omega_0}. The contribution of low frequencies is
  determined by
  %%
  \deqn{
    C(U^{(i)},\omega_0)=\frac{
      \sum_{\omega \leq \omega_0}\Pi(\omega)}{
      \sum_{\omega \leq \frac{1}{2}}\Pi(\omega)}%% \text{.}
  }{
    C(U^i,omega) = sum_{omega <= omega_0}{Pi(omega)} / sum_{omega <=
      0.5}{Pi(omega)}
  }
%   %%
  \eqn{U^{(i)}}{U^i} are considered to represent a trend component if
  \eqn{C(U^{(i)},\omega_0) \leq C_0 ,\; \epsilon \:
  [0,1]}{C(U^i,omega_0) <= C0} (see argument \code{C0}).

  The selection of \strong{\emph{periodic}} components is based on the fact that
  such components produce pairs of Eigenvalues and Eigenvectors
  (\eqn{U^{(i)}}{U^i}, \eqn{U^{(i+1)}}{U^i+1}). Such pairs have almost
  equal eigenvectors and the dominant frequency
  %%
  \deqn{
    \omega_*^{(i)} = \arg\max_\omega\{\Pi(\omega)\}
  }{
    omega*^i = argmax {Pi(omega)}
  }
  %%
  of the Eigenvectors is similar. This is the case if the following
  condition is met:
  %%
  \deqn{
    \delta_{\omega_*^{(i)}} = \|\omega_*^{(i)} - \omega_*^{(i+1)} \|
    \leq \frac{1}{L^2}
    %%\text{.}
  }{
    delta^i = | omega*^i - omega*^(i+1) |
  }
  %%
  (As the \eqn{\Pi(\omega)}{Pi(omega)} is estimated via \code{\link{fft}},
  \eqn{\delta_{\omega_*^{(i)}}}{delta^i} has to be \eqn{0} to fulfill
  that criterion, meaning that \eqn{\omega_*^{(i)}}{omega*^i} and
  \eqn{\omega_*^{(i+1)}}{omega*^(i+1)} are equal.)

  The above criterion may produce so called "false alarms". Therefore it
  is in addition necessary to check whether the Eigenvectors \eqn{U^{(i)}}{U^i}
  are really governed by one single frequency. This is done by checking
  whether the dominant frequency explains enough of the variance of
  \eqn{\Pi(\omega)}{Pi(omega)}. This is checked by the criterion
  %%
  \deqn{
    R^{(i)}=\max_\omega \left\{\frac{1}{2}(\gamma^{(i)}(\omega) +
    \gamma^{(i)}(\omega + \frac{1}{L}))\right\}
  }{
    R^i = max{0.5( gamma^i(omega) + gamma^i(omega + 1/L))}
  }
  %%
  where
  %%
  \deqn{
    \gamma^{(i)}(\omega_k) =
    \frac{\Pi^{(i)}(\omega_k)}{\sum_{\omega} \Pi^{(i)}(\omega)} +
    \frac{\Pi^{(i+1)}(\omega_{k})}{\sum_{\omega} \Pi^{(i+1)}(\omega)}
    %%\text{.}
  }{
    gamma^i(omega_k) = Pi^i(omega_k)/sum{Pi^i(omega)}   +
    Pi^(i+1)(omega_k)/sum{Pi^(i+1)(omega)}
  }
%   %%
  The term \eqn{\gamma^{(i)}(\omega + \frac{1}{L}}{gamma^i(omega + 1/L)}
  is to account for smearing of the periodogram in the case of signals
  with frequencies that are not multiples of
  \eqn{1/L}. \eqn{U^{(i)}}{U^i} and \eqn{U^{(i+1)}}{U^(i+1)} are
  considered to represent a quasi periodic  component if \eqn{R^{(i)}
  \leq  r_0,\; \epsilon \: [0,1]}{ R^i <= r_0} (see argument \code{r0}).

  If SSA Eigenelements are detected as contributing to the trend and the
  (quasi) periodic components they are only assigned to the periodicities.
}
\value{
  An object of class \emph{SSAsignal}
  \item{periodic}{A (possibly empty) data.frame. Each row is associated
  with a quasi periodic component, containing the following informations in its columns:
    \begin{itemize}
    \item{\strong{\code{k1}},\strong{\code{k2}}}{The rank of the Eigenelements
      contributing to the quasi periodic component}
    \item{\strong{\code{freq}}}{The frequency of the maximal value of the power
      spectrum of the corresponding Eigenvectors (dominant frequency).}
    \item{\strong{\code{smoothness}}}{Approximately the proportion of variance explained
      by the dominant frequency.}
    \item{\strong{\code{explainedVariance}}}{Percentage of the variance explained by the
      Eigenvalues with the rank \code{k1},\code{k2}}
    \end{itemize}
  }
  %%
  \item{trend }{A (possibly empty) data.frame. Each row is associated
    with a trend component, containing the following
    informations in its columns:
    \begin{itemize}
    \item{\strong{\code{k}}}{The rank of the Eigenelements contributing to the trend
      component}
    \item{\strong{\code{smoothness}}}{The proportion of variance explained by the
      frequencies \eqn{\leq} \code{omega0}.}
    \item{\strong{\code{explainedVariance}}}{Percentage of the variance explained by the
      Eigenvalues with the rank \code{k}}
    \end{itemize}
  }
%%
  \item{trendAndPeriod }{A (possibly empty) vector kontaining the rank
    numbers of the Eigenelements that have been identified as
    contributing to trends and quasi periodic components. These elements
    are all considered to correspond to periodicities and are therefore
    not listed in \code{trend}.}
  \item{call }{Call of the generating function}
}
\references{
  Alexandrov, T. A Method of Trend Extraction Using Singular Spectrum
  Analysis ArXiv e-prints, 2008, 804.

  Alexandrov, T. A method of extraction of quasi-periodic time series
  components using Singular Spectrum Analysis 2008.

  Vautard, R.; Yiou, P. & Ghil, M. Singular-spectrum analysis: a toolkit
  for short, noisy chaotic signals Physica D-Nonlinear Phenomena,
  Elsevier North-Holland, Inc. New York, NY, USA, 1992, 58, 95-126.
}
\author{ Lukas Gudmundsson }


\seealso{ \code{\link{decompSSA}}, \code{\link{decompSSAM}},
\code{\link{plot.sdTest}},\code{\link{MCSSA}} }
\examples{

## Artifitial Time series
set.seed(120)
n <- arima.sim(n=200,alist(ar=0.4))
n <- as.vector(scale(n))
x <- sin(seq(0,20*pi,len=200)) + 0.01*(1:200)
## To enable the interpretation of the eigenvectors
## as a representation of explained variance: remove the mean
xn <- scale(x+n/2,center=TRUE,scale=FALSE)
xn.dc <- decompSSA(xn,L=40)

## Identify signal
## It may happen that no or spurious signals are detected due to noise
(xn.signal <- getSignal(xn.dc,C0=0.85,r0=0.85))

## reconstruct signal (trend and period)
trendInd <- as.numeric(xn.signal$trend[,1])
perInd <- as.numeric(as.matrix(xn.signal$period[,1:2]))
xn.trend <- reconSSA(xn.dc,xn,list(trendInd))
xn.period <- reconSSA(xn.dc,xn,list(perInd))

plot(xn,t="l")
lines(xn.trend,lwd=2,col="red")
lines(xn.period,lwd=2,col="blue")
legend("bottomright",legend=c("Trend","Period","xn"),
lwd=c(2,2,1),col=c("red","blue","black")) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
