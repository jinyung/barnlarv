plot_pca_tps <- function(nef.coef, pca, har, n = 40,
                          fun = c("sd", "2sd", "range"), pc = 2,
                          ev, score, plotlayout = TRUE, ...) {
  fun <- match.arg(fun)
  coef.len <- dim(nef.coef)[2] / 4

  # coefficient of mean shape
  mshape.coef <- apply(nef.coef[, c(1:har, 1:har + coef.len,
                                    1:har + coef.len * 2,
                                    1:har + coef.len * 3)], 2,  mean)

  # get the loadings of the PC
  if(!missing(pca)) {
    rot <- pca$rotation
    x <- pca$x
    har <- (dim(pca$rotation)[1] + 3) /4
  } else {
    rot <- ev
    x <- score
  }

  if(!missing(pca)) {
    # compute the extreme or sd coefficient values on PC
    nhar <- har - 1
    pc1.rot <- c(0, rot[1:nhar, 1], 0, rot[1:nhar+nhar, 1],
                 0, rot[1:nhar+nhar*2, 1], rot[1:har+nhar*3, 1])
    if (pc >= 2)
      pc2.rot <- c(0, rot[1:nhar, 2], 0, rot[1:nhar+nhar, 2],
                   0, rot[1:nhar+nhar*2, 2], rot[1:har+nhar*3, 2])
  } else {
    pc1.rot <- c(rot[1:har, 1], rot[1:har+coef.len, 1],
                 rot[1:har+coef.len*2, 1], rot[1:har+coef.len*3, 1])
    if (pc >= 2)
      pc2.rot <- c(rot[1:har, 2], rot[1:har+coef.len, 2],
                   rot[1:har+coef.len*2, 2], rot[1:har+coef.len*3, 2])
  }

  if (fun == "range") {
    pc1.plus.score <- max(x[, 1])
    pc1.minus.score <- min(x[, 1])
    if (pc >= 2) {
      pc2.plus.score <- max(x[, 2])
      pc2.minus.score <- min(x[, 2])
    }
  } else {
    times <- ifelse (fun == "2sd", 2, 1)
    pc1.plus.score <- mean(x[, 1]) + times * sd(x[, 1])
    pc1.minus.score <- mean(x[, 1]) - times * sd(x[, 1])
    if (pc > 1) {
      pc2.plus.score <- mean(x[, 2]) + times * sd(x[, 2])
      pc2.minus.score <- mean(x[, 2]) - times * sd(x[, 2])
    }
  }

  pc1.plus.coef <- mshape.coef + pc1.plus.score * pc1.rot
  pc1.minus.coef <- mshape.coef + pc1.minus.score * pc1.rot

  if (pc >= 2) {
    pc2.plus.coef <- mshape.coef + pc2.plus.score * pc2.rot
    pc2.minus.coef <- mshape.coef + pc2.minus.score * pc2.rot
  }

  if (pc >= 3) {
    topleft.coef <- mshape.coef + pc1.minus.score * pc1.rot + pc2.plus.score * pc2.rot
    topright.coef <- mshape.coef + pc1.plus.score * pc1.rot + pc2.minus.score * pc2.rot
    bottomleft.coef <- mshape.coef + pc1.minus.score * pc1.rot + pc2.minus.score * pc2.rot
    bottomright.coef <- mshape.coef + pc1.plus.score * pc1.rot + pc2.minus.score * pc2.rot
  }

  # inverse fourier to get the outline
  n.point <- dim(nef.coef)[2] / 2
  if (missing(har))
    har <- dim(nef.coef)[2] / 4
  mshape_co <- .iefourier_co(mshape.coef, har, n = n.point)
  pc1_max_co <- .iefourier_co(pc1.plus.coef, har, n = n.point)
  pc1_min_co <- .iefourier_co(pc1.minus.coef, har, n = n.point)
  if (pc >= 2) {
    pc2_max_co <- .iefourier_co(pc2.plus.coef, har, n = n.point)
    pc2_min_co <- .iefourier_co(pc2.minus.coef, har, n = n.point)
  }
  if(pc == 3) {
    topleft.co <- .iefourier_co(topleft.coef, har, n = n.point)
    topright.co <- .iefourier_co(topright.coef, har, n = n.point)
    bottomleft.co <- .iefourier_co(bottomleft.coef, har, n = n.point)
    bottomright.co <- .iefourier_co(bottomright.coef, har, n = n.point)
  }

  # plot
  master <- .plotf(nef.coef, plot = FALSE)
  empty.plot <- expression(plot(1, 1, type = "n", axes = F, xlab = "",
                                ylab = ""))
  pc1 <- expression({
    .tps2(mshape_co, pc1_min_co, n, master, main = "-PC1")
    polygon(pc1_min_co, ...)
    .tps2(mshape_co, mshape_co, n, master, main = "Mean")
    polygon(mshape_co, ...)
    .tps2(mshape_co, pc1_max_co, n, master, main = "+PC1")
    polygon(pc1_max_co, ...)
  })

  pc2 <- expression({
    .tps2(mshape_co, pc2_max_co, n, master, main = "+PC2")
    polygon(pc2_max_co, ...)
    .tps2(mshape_co, pc2_min_co, n, master, main = "-PC2")
    polygon(pc2_min_co, ...)
  })

  par(mar = c(0, 0, 3, 0))
  if (pc > 1) {
    if(plotlayout)
      layout(matrix(c(6, 4, 7, 1:3, 8, 5, 9), 3, 3, byrow = TRUE))
    eval(pc1)
    eval(pc2)
    if (pc == 2) {
      eval(rep(empty.plot, 4))
    } else {
      .tps2(mshape_co, topleft.co, n, master, main = "-PC1/+PC2")
      polygon(topleft.co, ...)
      .tps2(mshape_co, topright.co, n, master, main = "+PC1/+PC2")
      polygon(topright.co, ...)
      .tps2(mshape_co, bottomleft.co, n, master, main = "-PC1/-PC2")
      polygon(bottomleft.co, ...)
      .tps2(mshape_co, bottomright.co, n, master, main = "+PC1/-PC2")
      polygon(bottomright.co, ...)
    }
  } else {
    layout(pc.layout <- matrix(c(1:3), 1, 3))
    eval(pc1)
  }
}
