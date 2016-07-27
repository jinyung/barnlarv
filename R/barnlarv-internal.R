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

# inverse fourier, taken from Julien Claude's "Morphometrics with R" (2008)
.iefourier <- function(an,bn,cn,dn,k,n,ao=0,co=0) {
  theta<-seq(0,2*pi, length=n+1)[-(n+1)]
  harmx <- matrix (NA, k, n)
  harmy <- matrix (NA, k, n)
  for (i in 1:k){
    harmx[i,]<-an[i]*cos(i*theta)+bn[i]*sin(i*theta)
    harmy[i,]<-cn[i]*cos(i*theta)+dn[i]*sin(i*theta)}
  x<-(ao/2) + apply(harmx, 2, sum)
  y<-(co/2) + apply(harmy, 2, sum)
  list(x=x, y=y)
}

# wrapper for iefourier
.iefourier_co <- function (coef, har, n) {
  coef.len <- length(coef) / 4  # length of each harmonic series
  out <- .iefourier(coef[1:har], coef[1:har + coef.len],
                    coef[1:har + coef.len * 2],
                    coef[1:har + coef.len * 3], har, n)
  # iefourier output is xy coords, changed to matrix
  out <- matrix(unlist(out), length(out$x), 2)
}

# plot thin plate spline deformation, modified from Julien Claude's
# "Morphometrics with R" (2008)
.tps2 <-function(matr, matt, n, master, ...) {
  xm<-min(matt[,1])
  ym<-min(matt[,2])
  xM<-max(matt[,1])
  yM<-max(matt[,2])
  rX<-xM-xm; rY<-yM-ym
  a<-seq(xm-1/5*rX, xM+1/5*rX, length=n)
  b<-seq(ym-1/5*rX, yM+1/5*rX,by=(xM-xm)*7/(5*(n-1)))
  m<-round(0.5+(n-1)*(2/5*rX+ yM-ym)/(2/5*rX+ xM-xm))
  M<-as.matrix(expand.grid(a,b))
  ngrid<-.tps2d(M,matr,matt)
  r1<-range(master[,1,]); r2<-range(master[,2,])
  r1.1<- abs(r1[2]-r1[1]); r2.1<- abs(r2[2]-r2[1])
  xlim<- c(r1[1]-r1.1/4, r1[2]+r1.1/4)
  ylim<- c(r2[1]-r2.1/4, r2[2]+r2.1/4)
  plot(ngrid, cex=0.2,asp=1,axes=F,xlab="",ylab="", col="gray",
       xlim=xlim, ylim=ylim, ...)
  for (i in 1:m){lines(ngrid[(1:n)+(i-1)*n,], col="gray")}
  for (i in 1:n){lines(ngrid[(1:m)*n-i+1,], col="gray")}
}

# calculate thin plate spline deformation, taken from Julien Claude's
# "Morphometrics with R" (2008)
.tps2d <- function(M, matr, matt) {
  p<-dim(matr)[1]; q<-dim(M)[1]; n1<-p+3
  P<-matrix(NA, p, p)
  for (i in 1:p)
  {for (j in 1:p){
    r2<-sum((matr[i,]-matr[j,])^2)
    P[i,j]<- r2*log(r2)}}
  P[which(is.na(P))]<-0
  Q<-cbind(1, matr)
  L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,3,3)))
  m2<-rbind(matt, matrix(0, 3, 2))
  coefx<-solve(L)%*%m2[,1]
  coefy<-solve(L)%*%m2[,2]
  fx<-function(matr, M, coef)
  {Xn<-numeric(q)
  for (i in 1:q)
  {Z<-apply((matr-matrix(M[i,],p,2,byrow=T))^2,1,sum)
  Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*(Z*log(Z)))}
  Xn}
  matg<-matrix(NA, q, 2)
  matg[,1]<-fx(matr, M, coefx)
  matg[,2]<-fx(matr, M, coefy)
  matg
}

# plot array of outlines from coefficients of EFA harmonics
.plotf <- function(coef, har, n, meancol, plot = TRUE) {
  coef.len <- dim(coef)[2] / 4  # length of each harmonic series
  out.len <- dim(coef)[1] # number of outlines
  if (missing(n))
    n <- dim(coef)[2] / 2 # number of points to be interpolated
  if (missing(har))
    har <- coef.len
  if (missing(meancol))
    meancol = par()$col

  # compute mean shape from mean coefficients
  meancoef <- apply(coef, 2, mean)
  meansh <- .iefourier_co(coef = meancoef, har = har, n = n)

  # compute outlines
  all.out <- array(NA, c(dim(coef)[2] / 2, 2, out.len))
  for(i in 1:out.len) {
    all.out[, , i] <- .iefourier_co(coef = coef[i, ], har = har, n = n)
    # all.out[, 1, i] <- temp$x; all.out[, 2, i] <- temp$y
  }

  # plot
  if (plot == TRUE) {
    xlim <- c(min(all.out[, 1, ]), max(all.out[, 1, ]))
    ylim <- c(min(all.out[, 2, ]), max(all.out[, 2, ]))
    par(mar = rep(0, 4))
    plot(meansh[, 1], meansh[, 2], asp = 1, type = "n",
         axes = F, xlim = xlim, ylim = ylim)
    for (i in 1:out.len)
      polygon(all.out[, 1, i], all.out[, 2, i], border = rgb(.5,.5,.5,.5))

    # add on the meanshape outline
    polygon(meansh[, 1], meansh[, 2], border = meancol, lwd = 3)
  }

  # output
  invisible(all.out)
}

# Note: this works specifically for outline reconstruction without first
# harmonics BUT retain the d_1 coefficient
.plot_pca_tps <- function(nef.coef, pca, har, n = 40) {
  coef.len <- dim(nef.coef)[2] / 4

  # coefficient of mean shape
  mshape_coef <- apply(nef.coef[, c(1:har, 1:har + coef.len,
                                    1:har + coef.len * 2,
                                    1:har + coef.len * 3)], 2,  mean)

  # get the loadings of the PC
  rot <- pca$rotation

  # compute the extreme coefficient values on PC
  nhar <- har - 1
  pc1_max_coef <- mshape_coef + max(pca$x[, 1]) * c(0, rot[1:nhar, 1],
                                                    0, rot[1:nhar+nhar, 1],
                                                    0, rot[1:nhar+nhar*2, 1],
                                                    rot[1:har+nhar*3, 1])
  pc1_min_coef <- mshape_coef + min(pca$x[, 1]) * c(0, rot[1:nhar, 1],
                                                    0, rot[1:nhar+nhar, 1],
                                                    0, rot[1:nhar+nhar*2, 1],
                                                    rot[1:har+nhar*3, 1])
  pc2_max_coef <- mshape_coef + max(pca$x[, 2]) * c(0, rot[1:nhar, 2],
                                                    0, rot[1:nhar+nhar, 2],
                                                    0, rot[1:nhar+nhar*2, 2],
                                                    rot[1:har+nhar*3, 2])
  pc2_min_coef <- mshape_coef + min(pca$x[, 2]) * c(0, rot[1:nhar, 2],
                                                    0, rot[1:nhar+nhar, 2],
                                                    0, rot[1:nhar+nhar*2, 2],
                                                    rot[1:har+nhar*3, 2])

  # inverse fourier to get the outline
  n_point <- dim(nef.coef)[2] / 2
  if (missing(har))
    har <- dim(nef.coef)[2] / 4
  mshape_co <- .iefourier_co(mshape_coef, har, n = n_point)
  pc1_max_co <- .iefourier_co(pc1_max_coef, har, n = n_point)
  pc1_min_co <- .iefourier_co(pc1_min_coef, har, n = n_point)
  pc2_max_co <- .iefourier_co(pc2_max_coef, har, n = n_point)
  pc2_min_co <- .iefourier_co(pc2_min_coef, har, n = n_point)

  # plot
  master <- .plotf(nef.coef, plot = FALSE)
  layout(matrix(c(rep(rep(c(1, 2, 3), each = 2), 2),
                  rep(rep(c(4, 5, 6), each = 2), 2),
                  rep(rep(c(7, 8, 9), each = 2), 2)), 6, 6))
  par(mar = c(0, 0, 3, 0))
  plot(1, 1, type = "n", axes = F, xlab = "", ylab = "")
  .tps2(mshape_co, pc1_min_co, n, master, main = "-PC1")
  polygon(pc1_min_co, lwd = 2)
  plot(1, 1, type = "n", axes = F, xlab = "", ylab = "")
  .tps2(mshape_co, pc2_max_co, n, master, main = "+PC2")
  polygon(pc2_max_co, lwd = 2)
  .tps2(mshape_co, mshape_co, n, master, main = "Mean")
  polygon(mshape_co, lwd = 2)
  .tps2(mshape_co, pc2_min_co, n, master, main = "-PC2")
  polygon(pc2_min_co, lwd = 2)
  plot(1, 1, type = "n", axes = F, xlab = "", ylab = "")
  .tps2(mshape_co, pc1_max_co, n, master, main = "+PC1")
  polygon(pc1_max_co, lwd = 2)
  plot(1, 1, type = "n", axes = F, xlab = "", ylab = "")
}
