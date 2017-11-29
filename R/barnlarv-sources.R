# modified from 'otolith' package
.loadimg <- function (path, resize = FALSE, width) {
  img <- tiff::readTIFF(path)
  if (EBImage::numberOfFrames(img) > 1)
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
    landidx <- star[i] + round((0: (p-1)) * (outlength / p))
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
.iefourier_co <- function (coef, har, n, ...) {
  coef.len <- length(coef) / 4  # length of each harmonic series
  out <- .iefourier(coef[1:har], coef[1:har + coef.len],
                    coef[1:har + coef.len * 2],
                    coef[1:har + coef.len * 3], k = har, n = n, ...)
  # iefourier output is xy coords, changed to matrix
  out <- matrix(unlist(out), length(out$x), 2)
}

# plot array of outlines from coefficients of EFA harmonics
.plotf <- function(coef, har, n, meancol, indcol, plot = TRUE) {
  coef.len <- dim(coef)[2] / 4  # length of each harmonic series
  out.len <- dim(coef)[1]  # number of outlines
  if (missing(n))
    n <- dim(coef)[2] / 2  # number of points to be interpolated
  if (missing(har))
    har <- coef.len
  if (missing(meancol))
    meancol = par()$col
  if (missing(indcol))
    indcol = adjustcolor("gray40", 0.01)

  # compute mean shape from mean coefficients
  meancoef <- apply(coef, 2, mean)
  meansh <- .iefourier_co(coef = meancoef, har = har, n = n)

  # compute outlines
  all.out <- array(NA, c(n, 2, out.len))
  for(i in 1:out.len) {
    all.out[, , i] <- .iefourier_co(coef = coef[i, ], har = har, n = n)
  }

  # plot
  if (plot == TRUE) {
    xlim <- c(min(all.out[, 1, ]), max(all.out[, 1, ]))
    ylim <- c(min(all.out[, 2, ]), max(all.out[, 2, ]))
    orimar <- par()$mar
    par(mar = rep(0, 4))
    plot(meansh[, 1], meansh[, 2], asp = 1, type = "n",
         axes = F, xlim = xlim, ylim = ylim)
    for (i in 1:out.len)
      polygon(all.out[, 1, i], all.out[, 2, i],
              border = indcol)
    par(mar = orimar)
    # add on the meanshape outline
    polygon(meansh[, 1], meansh[, 2], border = meancol, lwd = 1.5)
  }

  # output
  invisible(all.out)
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
  plot(ngrid, cex=0.2,asp=1,axes=FALSE,xlab="",ylab="", col="gray",
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

# for Rmarkdown table format
# https://stackoverflow.com/a/29953844
.removeNA <- function(x) {
  x <- gsub('\\bNA\\b', '  ', x)
  cat(x, sep='\n')
}

# calculate elliptic Fourier, from Julien Claude's "Morphometrics with R" (2008)
.efourier<-function(M, n=dim(M)[1]/2)
{p<-dim(M)[1]
Dx<-M[,1]-M[c(p,(1:p-1)),1]
Dy<-M[,2]-M[c(p,(1:p-1)),2]
Dt<-sqrt(Dx^2+Dy^2)
t1<-cumsum(Dt)
t1m1<-c(0, t1[-p])
T<-sum(Dt)
an<-bn<-cn<-dn<-numeric(n)
for (i in 1:n){
  an[i]<- (T/(2*pi^2*i^2))*sum((Dx/Dt)*
                                 (cos(2*i*pi*t1/T)-cos(2*pi*i*t1m1/T)))
  bn[i]<- (T/(2*pi^2*i^2))*sum((Dx/Dt)*
                                 (sin(2*i*pi*t1/T)-sin(2*pi*i*t1m1/T)))
  cn[i]<- (T/(2*pi^2*i^2))*sum((Dy/Dt)*
                                 (cos(2*i*pi*t1/T)-cos(2*pi*i*t1m1/T)))
  dn[i]<- (T/(2*pi^2*i^2))*sum((Dy/Dt)*
                                 (sin(2*i*pi*t1/T)-sin(2*pi*i*t1m1/T)))}
ao<-2*sum(M[,1]*Dt/T)
co<-2*sum(M[,2]*Dt/T)
# c(an, bn,cn,dn)
list(ao=ao,co=co,an=an,bn=bn,cn=cn,dn=dn)
}



# dotchart2() modified from dotchat()
.dotchart2<- function (x, labels = NULL, groups = NULL, gdata = NULL, cex = par("cex"),
                      pch = 21, gpch = 21, bg = par("bg"), color = par("fg"), gcolor = par("fg"),
                      lcolor = "gray", xlim = range(x[is.finite(x)]), main = NULL,
                      xlab = NULL, ylab = NULL, ...)
{
  opar <- par("mai", "mar", "cex", "yaxs")
  on.exit(par(opar))
  par(cex = cex, yaxs = "i")
  if (!is.numeric(x))
    stop("'x' must be a numeric vector or matrix")
  n <- length(x)
  if (is.matrix(x)) {
    if (is.null(labels))
      labels <- rownames(x)
    if (is.null(labels))
      labels <- as.character(1L:nrow(x))
    labels <- rep_len(labels, n)
    if (is.null(groups))
      groups <- col(x, as.factor = TRUE)
    glabels <- levels(groups)
  }
  else {
    if (is.null(labels))
      labels <- names(x)
    glabels <- if (!is.null(groups))
      levels(groups)
    if (!is.vector(x)) {
      warning("'x' is neither a vector nor a matrix: using as.numeric(x)")
      x <- as.numeric(x)
    }
  }
  plot.new()
  linch <- if (!is.null(labels))
    max(strwidth(labels, "inch"), na.rm = TRUE)
  else 0
  if (is.null(glabels)) {
    ginch <- 0
    goffset <- 0
  }
  else {
    ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
    goffset <- 0.4
  }
  if (!(is.null(labels) && is.null(glabels))) {
    nmai <- par("mai")
    nmai[2L] <- nmai[4L] + max(linch + goffset, ginch) +
      0.1
    par(mai = nmai)
  }
  if (is.null(groups)) {
    o <- 1L:n
    y <- o
    ylim <- c(0, n + 1)
  }
  else {
    o <- sort.list(as.numeric(groups), decreasing = TRUE)
    x <- x[o]
    groups <- groups[o]
    color <- rep_len(color, length(groups))[o]
    lcolor <- rep_len(lcolor, length(groups))[o]
    offset <- cumsum(c(0, diff(as.numeric(groups)) != 0))
    y <- 1L:n + 2 * offset
    ylim <- range(0, y + 2)
  }
  plot.window(xlim = xlim, ylim = ylim, log = "")
  lheight <- par("csi")
  if (!is.null(labels)) {
    linch <- max(strwidth(labels, "inch"), na.rm = TRUE)
    loffset <- (linch + 0.1)/lheight
    labs <- labels[o]
    mtext(labs, side = 2, line = loffset, at = y, adj = 0,
          col = color, las = 2, cex = cex, ...)
  }
  abline(h = y, lty = "dotted", col = lcolor)
  points(x, y, pch = pch, col = color, bg = bg)
  if (!is.null(groups)) {
    gpos <- rev(cumsum(rev(tapply(groups, groups, length)) +
                         2) - 1)
    ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
    goffset <- (max(linch + 0.2, ginch, na.rm = TRUE) + 0.1)/lheight
    mtext(glabels, side = 2, line = goffset, at = gpos, adj = 0,
          col = gcolor, las = 2, cex = cex, ...)
    if (!is.null(gdata)) {
      abline(h = gpos, lty = "dotted")
      points(gdata, gpos, pch = gpch, col = gcolor, bg = bg,
             ...)
    }
  }

  box()
  title(main = main, xlab = xlab, ylab = ylab, ...)
  return(y)
}
