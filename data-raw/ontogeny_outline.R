# =================== NOTES =========================================
# This is basically same as inst/extdata/stage2outlines.R,
# just a different dataset with minor modifications
# Remember to checkout stage2outlines.R first!
# ===================================================================

# ========== required packages ===============
require(EBImage)
require(otolith)
require(devtools)
require(splancs)

# ========== reading the files ===============
# morphological data
ontogeny_mdat <- read.csv(system.file("extdata/ontogeny",
                                      "frontal-horn-summary.csv",
                                      package = "barnlarv"))

# img list of the raw digitized images
drawing_dir <- system.file("extdata/ontogeny", "larvae_drawing",
                           package = "barnlarv")
ontogeny_list <- list.files(drawing_dir, full.names = TRUE, pattern = ".tif")

# number of images
imgnum <- length(ontogeny_list)

# initialize
ontogeny_outline <- vector("list", imgnum)

# Change names in list
names(ontogeny_outline) <- paste(ontogeny_mdat$genus, ontogeny_mdat$species,
                                 ontogeny_mdat$stage, sep = "-")

# ========== extracting outline from images ===============
# loop thru to extract outline
for (i in 1:length(ontogeny_list)) {
  cat("\rpicture no:", i, "/", imgnum)
  im <- .loadimg(ontogeny_list[i])
  im <- abs(im - 1) # inverse the image
  x <- dim(im)[1]
  y <- dim(im)[2]
  # using the scale provided in frontal-horn-summary.csv to standardize size
  rscale <- ontogeny_mdat$pic.scale[i] / 0.25 # all enlarged 4 times
  im <- resize(im, w = x * rscale)
  ontogeny_outline[[i]] <- otolith::extractout(im, threshold = 0.3, plot = "no")
}

# ========== post-processing of outlines ===============
# PC align the outline
ontogeny_outline <- lapply(ontogeny_outline, .aligne)

# change the direction because pc align can rotate it
# algorithm: posterior part of outline is smaller than anterior part
for (i in 1:imgnum) {
  mat <- ontogeny_outline[[i]]
  centroid <- apply(ontogeny_outline[[i]], 2, mean)
  # a different algorithm used here
  left <- mat[mat[, 1] < centroid[1], ]
  right <- mat[mat[, 1] > centroid[1], ]
  if (diff(range(left[, 2])) > diff(range(right[, 2]))) {
    ontogeny_outline[[i]] <- -ontogeny_outline[[i]]
  }
}

# Above algorithm did not work for outline of stage I O.cor (horns not extended)
# Need to change manually
idx <- which(names(ontogeny_outline) %in% "Octolasmis-cor-1")
ontogeny_outline[[idx]] <- -ontogeny_outline[[idx]]

# ========== output outlines ===============
# save the outline to data/
devtools::use_data(ontogeny_outline, overwrite = TRUE)
