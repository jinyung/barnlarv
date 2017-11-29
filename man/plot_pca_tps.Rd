\name{plot_pca_tps}
\alias{plot_pca_tps}
\title{outline reconstruction from EFA cofficients predicted by PCA}
\description{Reconstruct and plot the shape changes in outline along
  principle components (PC)}
\usage{
plot_pca_tps(nef.coef, pca, har, n = 40, fun = c("sd", "2sd", "range"), pc = 2,
             ev, score, plotlayout = TRUE, ...)
}
\arguments{
  \item{nef.coef}{(matrix) coefficients from normalized elliptic Fourier
    analysis}
  \item{pca}{(class "prcomp" object) result from \code{\link[stats]{prcomp}}
    principal component analysis}
  \item{har}{(integer) number of harmanics}
  \item{n}{(integer) number of coordinates of outline to interpolate}
  \item{fun}{(character) type of range along PC}
  \item{pc}{(integer) number of PC to reconstruct, 1 or 2}
  \item{ev}{(numeric, optional) the eigennvector if \code{pca} not provided}
  \item{score}{(numeric, optional) the rotated data/PC scores if
    \code{pca} not provided}
  \item{plotlayout}{(logical) if TRUE plot in single plot window with default
  layout, else multiple single plots}
  \item{...}{(optional) optional argument to be passed to
    \code{\link[graphics]{polygon}} to control look of the reconstructed
    outline}
}
\details{
  A wrapper for functions developed in Claude (2008) for visualization of
  principal component analysis of normalized elliptic Fourier analysis of
  outline shapes. Note that for inverse elliptic Fourier
  transform, first three coefficients of the first harmonics are not used
  (contain no meaningful information), but the \eqn{d_1} coefficient is
  retained.

  This function can also be used to
  reconstruct outline shape changes along common allometric component/ residual
  shape components described in Mitteroecker et al. (2004).}
\value{Plot(s) of reconstructed outlines along extremes of PC(s) shown as
  thin-plate-spline deformations from the overall mean shape.}
\references{
  Claude, J. (2008). \emph{Morphometrics with R}. New York:Springer.

  Mitteroecker, P., Gunz, P., Bernhard, M., Schaefer, K., & Bookstein, F. L. (2004). Comparison of cranial ontogenetic trajectories among great apes and humans. \emph{Journal of Human Evolution}, 46(6), 679-698.}
