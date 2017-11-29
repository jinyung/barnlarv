\name{barnlarv-sources}
\alias{barnlarv-sources}
\alias{.aligne}
\alias{.centsiz}
\alias{.dotchart2}
\alias{.ed}
\alias{.efourier}
\alias{.extract_landmark}
\alias{.iefourier}
\alias{.iefourier_co}
\alias{.loadimg}
\alias{.plotf}
\alias{.removeNA}
\alias{.tps2}
\alias{.tps2d}
\title{Sourced/internal functions}
\description{
  Functions imported from other sources and small miscellaneous functions,
  including:
  \describe{
    \item{\code{.aligne}}{align outline/landmark configuration to first
      principal component}
    \item{\code{.centsiz}}{calculate centroid size}
    \item{\code{.dotchart2}}{modified \code{\link[graphics]{dotchart}}function}
    \item{\code{.ed}}{calculate euclidean distance}
    \item{\code{.efourier}}{elliptic Fourier analysis}
    \item{\code{.extract_landmark}}{sample semi-landmark along outline}
    \item{\code{.iefourier}}{inverse elliptic Fourier analysis}
    \item{\code{.iefourier_co}}{a wrapper for \code{.iefourier}}
    \item{\code{.loadimg}}{read image}
    \item{\code{.plotf}}{plot arrays of outlines from coefficients of EFA
      harmonics }
    \item{\code{.removeNA}}{for markdown table formatting}
    \item{\code{.tps2d}}{calculate thin plate spline deformation}
    \item{\code{.tps2}}{plot thin plate spline deformation}
  }
}
% \usage{
% .aligne(A)
% .centsiz(M)
% .dotchart2(x, labels = NULL, groups = NULL, gdata = NULL, cex = par("cex"),
%                       pch = 21, gpch = 21, bg = par("bg"), color = par("fg"),
%                       gcolor = par("fg"), lcolor = "gray",
%                       xlim = range(x[is.finite(x)]), main = NULL,
%                       xlab = NULL, ylab = NULL, ...)
% .ed(E, F)
% .efourier(M, n = dim(M)[1]/2)
% .extract_landmark(outline, p = 200, k = 2, star)
% .iefourier(an, bn, cn, dn, k, n, ao = 0, co = 0)
% .iefourier_co(coef, har, n, ...)
% .loadimg(path, resize = FALSE, width)
% .plotf(coef, har, n, meancol, indcol, plot = TRUE)
% .removeNA(x)
% .tps2(matr, matt, n, master, ...)
% .tps2d(M, matr, matt)
% }
\details{
Functions that starts with a period belong here. Most of the undocumented
functions here are sourced from 1) Claude (2008) book and 2) Wong et. al. (2016)
package '\code{otolith}' and the derivatives from these sources. Explanations of these function are provided in the comments of the source codes in

\code{file.path(system.file("R", package = "barnlarv"), "barnlarv-sources.R")}

and further explanations should refer to the references stated in the comments.
}

\references{
  Claude, J. (2008). \emph{Morphometrics with R}. New York:Springer.

  Wong, J. Y., Chu, C., Chong, V. C., Dhillon, S. K., & Loh, K. H. (2016). Automated otolith image classification with multiple views: an evaluation on Sciaenidae. \emph{Journal of Fish Biology}, 89(2), 1324-1344.

}
%\keyword{internal}
