\name{shapesc}
\alias{shapesc}
\title{Calculate shape score}
\description{one-dimensional shape score from multivariate shape regression, as from Drake & Klingenberg (2008)}
\usage{
shapesc(model, pred = FALSE)
}
\arguments{
  \item{model}{(class "lm" object) regression model from
    \code{\link[stats]{lm}}}
  \item{pred}{(logical) if \code{TRUE} returns projection of model predicted
    values}
}
\details{
  Originally for allometric regression of multivariate shape data,
  e.g. \code{lm} model of \code{lm(shape~size)}, but can be used for any other
  covariate of interest, see reference for more explanation.}
\value{Returns the projected shape scores}
\references{
  Drake, A. G., & Klingenberg, C. P. (2008). The pace of morphological change: historical transformation of skull shape in St Bernard dogs. \emph{Proceedings of the Royal Society of London B: Biological Sciences}, 275(1630), 71-76.
  }
