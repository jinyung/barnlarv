\name{rv}
\alias{vector_angle}
\alias{rv}
\alias{rv.test}
\title{Vector correlation}
\description{
  Calculate correlation or angle between vectors from shape regression on
  covariate and principal component analysis/ other type of analysis
}
\usage{
vector_angle(V1, V2)
rv(V1, V2)
rv.test(shape, V1, pc, V2, perm = 999, seed = 8888)
}
\arguments{
  \item{shape}{(matrix) shape data, individuals in rows, see details}
  \item{V1}{(numeric) numeric vector, see details}
  \item{pc}{(integer) the n-th principal component used to find correlation}
  \item{V2}{(numeric) numeric vector, see details}
  \item{perm}{(integer) number of permutation to perform}
  \item{seed}{(integer) seed used for random number generation for permuation}
}
\details{
  \code{vector_angle} and \code{rv} are generalized functions to calculate angle
  and correlation between two vectors, respectively. The vectors can be anything   but should be of same length.
  \code{rv.test} was written specifically to test the correlation between a
  regression vector of \code{shape~V1} and a specific principal component (PC)
  from principal component analysis of the shape data. \code{shape} can be
  coefficients from normalized elliptic Fourier analysis or Procrustes aligned
  coordinates). In \code{shape~V1}, \code{V1} is a covariate which one wish to
  check its correlation with a specific \code{pc}. Optionally, if \code{pc} is
  not provided, \code{V2} should be provided, such as regression vector from
  another \code{shape~covariate} regression or residual shape component
  from common allometric component analysis (Mitteroecker et al., 2004).
  Note that absolute
  value of correlation is returned by \code{rv.test}, as direction of PC is not
  considered when testing correlation with PC. Permutational test
  is used for testing the statistical significance of correlation.}
\value{
  \code{vector_angle} returns angle in radian, \code{rv} returns correlation in   range of -1 to 1. \code{rv.test} returns a list of:
    \item{rv.obs}{correlation in range of 0-1}
    \item{theta.obs}{angle in degree}
    \item{p.val}{p-value from permutation test}
}
\references{
  Mitteroecker, P., Gunz, P., Bernhard, M., Schaefer, K., & Bookstein, F. L. (2004). Comparison of cranial ontogenetic trajectories among great apes and humans. \emph{Journal of Human Evolution}, 46(6), 679-698.
}
