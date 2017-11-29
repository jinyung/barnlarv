\name{semilandmark}
\docType{data}
\alias{semilandmark}
\alias{stage2landmark}
\alias{ontogeny_landmark}
\title{Sampled semi-landmarks on outlines of barnacle nauplii}
\description{
  There are two datasets of semi-landmarks:
  \describe{
    \item{\code{stage2landmark}}{semi-landmarks sampled from
      \code{\link{stage2outline}}, outline from 102 species of
        stage II barnacle nauplii, gathered from publised description and lab
        collection.}
    \item{\code{ontogeny_landmark}}{semi-landmarks sampled from
      \code{\link{ontogeny_outline}}, outline from 12 species of
        barnacle nauplii across developmental stages.}
  }
}
\usage{
  data(stage2landmark)
  data(ontogeny_landmark)}
\format{
  \describe{
    \item{\code{stage2landmark}}{an array with dimensions of (200, 2, 102),
      representing (200 semi-landmarks, 2 columns of xy coordinates,
      102 species)}
    \item{\code{ontogeny_landmark}}{an array with dimensions of (200, 2, 67),
      67 represents 9 species with 6 naupliar stages, 2 with 4 stages
      and 1 with 5 stages.}
  }
}
\seealso{\code{\link{outline}}}
\keyword{datasets}
