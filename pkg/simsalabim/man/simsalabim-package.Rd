\name{simsalabim-package}
\alias{simsalabim-package}
\alias{ssaps}
\docType{package}
\title{
A collection of methods for time series analysis and signal detection.
}
\description{
This package implements Singular System Analysis and Phase  Synchronization Indices.
}
\details{
The functions bundled in this package focus on two objectives: The application
of Singular Spectrum Analysis (SSA) and the computation of phase
synchronization measures. The working horses of this package are the functions
\code{\link{decompSSA}} and \code{\link{phaSyn}}. For both applications some
statistical testing procedures are provided. SSA is implemented in a standard
form and in a variant able to handle missing data points
\code{\link{decompSSAM}}.
Currently the entire package is implemented in R which leads to some performance issues.

}

\author{
Lukas Gudmundsson

Maintainer: Lukas Gudmundsson \email{lukas.gudmundsson@geo.uio.no}
}
\references{
Golyandina, N.; Nekrutkin, V. & Zhiglkilavskifi, A.
Analysis of Time Series Structure: SSA and Related Techniques. CRC Press, 2001

Ghil, M.; Allen, M.; Dettinger, M.; Ide, K.; Kondrashov, D.; Mann, M.;
Robertson, A.; Saunders, A.; Tian, Y.; Varadi, F. & others.
Advanced spectral methods for climatic time series. Rev. Geophys, 2002, 40, 1003

Pikovsky, A.; Rosenblum, M. & Kurths, J.
Synchronization: A Universal Concept in Nonlinear Sciences.
(Cambridge Nonlinear Science Series) Cambridge University Press, 2003
}
\keyword{ package }
\keyword{ ts }
