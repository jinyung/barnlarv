vector_angle <- function(V1, V2) {
  V1.norm <- sqrt(sum(V1^2))
  V2.norm <- sqrt(sum(V2^2))
  acos((V1 %*% V2) / (V1.norm * V2.norm))
}

rv <- function(V1, V2) {
  V1.norm <- sqrt(sum(V1^2))
  V2.norm <- sqrt(sum(V2^2))
  V1 <- V1 / V1.norm
  V2 <- V2 / V2.norm
  V1 %*% V2
}

rv.test <- function(shape, V1, pc, V2, perm = 999, seed = 8888) {
  lm.obs <- lm(shape ~ V1)
  reg.V1 <- coefficients(lm.obs)[2, ]
  if (!missing(pc)){
    pca.shape <- prcomp(shape)
    rot <- pca.shape$rotation[, pc]
  } else if (!missing(V2)) {
    rot <- V2
  } else {
    stop("provide pc or V2")
  }
  rv.obs <- abs(rv(reg.V1, rot))
  rv.perm <- NULL
  set.seed(seed)
  for (i in 1:perm) {
    V1.rand <- sample(V1)
    lm.perm <- lm(shape ~ V1.rand)
    reg.V1.perm <- coefficients(lm.perm)[2, ]
    rv.perm[i] <- abs(rv(reg.V1.perm, rot))
  }

  return(list(rv.obs = format(rv.obs, digits = 3, nsmall = 3),
              theta.obs = format(acos(rv.obs)* 180 / pi,
                                 digits = 1, nsmall = 1),
              p.val = sum(sapply(c(rv.perm, rv.obs), `>=`, rv.obs)) /
                (perm + 1)))
}
