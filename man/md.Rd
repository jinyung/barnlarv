\name{md}
\alias{md}
\alias{pd}
\title{Calculate morphological disparity (MD) and partial disparity (PD)}
\description{MD and PD as defined in Foote (1993)}
\usage{
  md(coeff)
  pd(coeff, idx)
}
\arguments{
  \item{coeff}{(matrix) shape data, can be coefficients from normalized elliptic     Fourier analysis or Procrustes aligned coordinates}
  \item{idx}{(integer) index to subset rows of shape data to calculate the
    partial disparity}
}
\details{
  \describe{
    \item{\code{md}}{calculate morphological disparity (MD)}
    \item{\code{pd}}{calculate partial disparity (PD)}
  }
}
\value{
  \code{md} returns a variance measure and \code{pd} returns a
  list of:
    \item{pd}{absolute PD}
    \item{prop.md}{relative PD in \% (proportion of MD)}
}
\references{
Foote, M. (1993). Contributions of individual taxa to overall morphological disparity. \emph{Paleobiology}, 19(4), 403-419.}
