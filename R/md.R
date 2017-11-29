# calculation of morphological disparity
md <- function(coeff) {
  N <- dim(coeff)[1]
  mean.coeff <- apply(coeff, 2, mean)
  vd <- apply(coeff, 1, `-`, mean.coeff)  # vector diff
  d.square <- apply(vd^2, 2, sum)  # squared distance
  MD <- sum(d.square) / (N - 1)
  MD
}

pd <- function(coeff, idx) {  # inherit all codes from .md
  N <- dim(coeff)[1]
  mean.coeff <- apply(coeff, 2, mean)
  vd <- apply(coeff, 1, `-`, mean.coeff)  # vector diff
  d.square <- apply(vd^2, 2, sum)  # squared distance
  MD <- sum(d.square) / (N - 1)
  PD <- sum(d.square[idx]) / (N - 1)
  return(list(pd = PD, prop.md = PD / MD))
}
