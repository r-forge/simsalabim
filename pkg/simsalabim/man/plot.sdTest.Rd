\name{plot.sdTest}
\alias{plot.sdTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots the results of the Shun and Duffy (1999) significance test.}
\description{
  Draws the SSA - Spectrum (eigenvalues versus the dominant frequency of the
  eigenvectors) and draws the the fitted noise model as well as the confidence
  bounds of the model.
}
\usage{
\method{plot}{sdTest}(x, asFreq = TRUE, normalize = TRUE, lam.pch = 19,
      lam.col = "black", lam.cex = 1, sig.pch = 19,
      sig.col = "red", sig.cex = 1, conf.col = "lightgray",
      conf.border = "lightgray", noise.col = "black",
      noise.lwd = 1, ann = TRUE, legend = TRUE,
      axes = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Output of \code{\link{sdTest}} }
  \item{asFreq}{  If \code{by="freq"}, whether frequency shall be plotted as
    frequency or as period.}
  \item{normalize}{Whether eigenvalues shall be normalised to their sum
                  or not.}
  \item{lam.pch}{\code{pch} of the eigenvalues. }
  \item{lam.col}{\code{col} of the eigenvalues. }
  \item{lam.cex}{\code{cex} of the eigenvalues.}
  \item{sig.pch}{\code{pch} of the eigenvalues, that are significantly different
                  from the noise spectrum.}
  \item{sig.col}{ \code{col} of the eigenvalues, that are significantly different
                  from the noise spectrum. }
  \item{sig.cex}{ \code{cex} of the eigenvalues, that are significantly different
                  from the noise spectrum.}
  \item{conf.col}{\code{col} of the confidence bounds of the noise spectrum,
                  as defined in \code{\link{polygon}}.}
  \item{conf.border}{ \code{border} of the confidence bounds of the noise spectrum,
                  as defined in \code{\link{polygon}}.}
  \item{noise.col}{ Colour of the line representing the fitted noise model.}
  \item{noise.lwd}{ \code{lwd} of the line representing the fitted noise model.}
  \item{ann}{ Whether to do annotations or not. }
  \item{legend}{ Whether to add a legend or not }
  \item{axes}{ Whether to add axes or not.}
  \item{\dots}{ Further arguments to be passed to other graphic functions. }
}
\author{Lukas Gudmundsson }
\seealso{ \code{\link{sdTest}}, \code{\link{par}} }
\keyword{ hplot }