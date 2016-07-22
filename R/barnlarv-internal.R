# =================== NOTES =============================
# data-raw/ cannot call internal functions, thus exported
# so this is not really internal functions now
# =======================================================

# modified from 'otolith' package
.loadimg <- function (path, resize = FALSE, width) {
  require(tiff)
  img <- readTIFF(path)
  if (numberOfFrames(img) > 1)
    img <- img[, , 1]
  img <- t(img)
  img <- EBImage::flip(img)
  invisible(img)
}

# modified from Julien Claude's "Morphometrics with R" (2008)
.aligne <- function (A) {
    B <- A
    Ms <- scale(A, scale = FALSE)
    sv <- eigen(var(Ms))
    M <- Ms %*% sv$vectors
    B <- M
    B
}

# extract semi-landmarks, modified from 'otolith' package
.extract_landmark <- function (outline, p = 200, k = 2, star) {
  # check
  if (length(star) != length(outline))
    stop("length of outline not match with length of star")

  # initialize empty array
  landmark <- array(data = NA, dim = c(p, k, length(outline)),
                    dimnames = list(NULL, c("x", "y"), names(outline)))

  # loop thru to get the semi-landmarks
  for (i in 1:length(outline)) {
    outlength <- dim(outline[[i]])[1]
    landidx <- star[i] + round((1:p) * (outlength / p))
    landidx[landidx > outlength] <- landidx[landidx > outlength] - outlength
    landmark[, , i] <- outline[[i]][landidx, ]
  }
  landmark
}

# euclidean distance, from Julien Claude's "Morphometrics with R" (2008)
.ed <- function (E, F) {
  sqrt(sum((E-F)^2))
}

# centroid size, from Julien Claude's "Morphometrics with R" (2008)
.centsiz <- function(M) {
  p <- dim(M)[1]
  size <- sqrt(sum(apply(M, 2, var)) * (p - 1))
  size
}
