# =================== NOTES =========================================
# All the raw data are available and can be accessed at:
# path/to/package/inst/extdata
# ===================================================================

# ========== required packages ===============
require(EBImage)
require(otolith)
# refer to https://github.com/jinyung/otolith/wiki/Installation-Guide
# for installation of these two packages
require(devtools)

# ========== reading the files ===============
# morphological data
stage2mdat <- read.csv(system.file("extdata/stage2", "frontal-horn-summary.csv",
                                   package = "barnlarv"))

# img list
drawing_dir <- system.file("extdata/stage2", "larvae_drawing", package = "barnlarv")
stage2list <- list.files(drawing_dir, full.names = TRUE, pattern = ".tif")

# ========== initialize ===============
# number of images
imgnum <- length(stage2list)
# empty list
stage2outline <- vector("list", imgnum)
# Change names in list
names(stage2outline) <- paste(stage2mdat$genus, stage2mdat$species)

# ========== extracting stage2outline from images ===============
# loop thru to extract stage2outline
for (i in 1:imgnum) {
  cat("\rphoto no:", i, "/", imgnum)
  im <- .loadimg(stage2list[i])
  im <- abs(im - 1) # inverse the image
  x <- dim(im)[1]
  y <- dim(im)[2]
  #rscale <- stage2mdat$pic.scale[i] / 0.273822563 # all enlarged
  # using the scale provided in frontal-horn-summary.csv to standardize size
  rscale <- stage2mdat$pic.scale[i] / 0.25 # all enlarged 4 times
  im <- resize(im, w = x * rscale)
  stage2outline[[i]] <- otolith::extractout(im, threshold = 0.3, plot = "no")
}

# ========== post-processing of outlines ===============
# PC align the stage2outline
stage2outline <- lapply(stage2outline, .aligne)

# change the direction because pc align can rotate it
# algorithm: posterior part of stage2outline is smaller than anterior part
for (i in 1:imgnum) {
  mat <- stage2outline[[i]]
  centroid <- apply(stage2outline[[i]], 2, mean)
  if (sum(mat[, 1] < centroid[1]) > sum(mat[, 1] > centroid[1])) {
    stage2outline[[i]] <- -stage2outline[[i]]
  }
}

# ========== output outlines ===============
# save to data/
devtools::use_data(stage2outline, overwrite = TRUE)
