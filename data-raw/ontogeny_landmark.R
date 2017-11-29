# =================== NOTES =========================================
# data-raw/ontogeny_outline.R was run prior to this
# This is basically same as inst/extdata/stage2landmark.R,
# just a different dataset, I was lazy to write a all-in-one wrapper,
# so you see a copy and paste here
# Remember to checkout stage2_landmark.R first!
# ===================================================================

# ========== required packages ===============
require(geomorph)
require(devtools)

# ========== call outlines for processing ===============
data(ontogeny_outline)

# ========== Get the starting point for semi-landmark extraction ===============
# get the starting point
ontogeny_star <- integer(length(ontogeny_outline))
for (i in 1:length(ontogeny_outline)) {
  ontogeny_star[i] <- which.min(ontogeny_outline[[i]][, 1])
}

# ***** NOTES: BELOW PARTS REQUIRE MANUAL ADJUSTMENT *****
# check for starting points
for (i in 1:length(ontogeny_outline)) {
  plot(ontogeny_outline[[i]], asp = 1, type = "l", xlab = "x", ylab = "y")
  points(ontogeny_outline[[i]][ontogeny_star[i], 1],
         ontogeny_outline[[i]][ontogeny_star[i], 2], col =3)
  centroid <- apply(ontogeny_outline[[i]], 2, mean)
  text(centroid[1], centroid[2], paste(i, names(ontogeny_outline[i])))
}

# ontogeny_outline with wrong starting points
# ***** manual adjustments commented out
# following outline starting point needed to be adjusted
# spidx <- c(grep("pilosella", names(ontogeny_outline)),
#            grep("polygenea", names(ontogeny_outline)),
#            grep("carcini", names(ontogeny_outline)),
#            grep("papillosus", names(ontogeny_outline)))

# correction: define start point as left half closest to the x-axis
# plot to check again
# for (i in spidx){
#   lefthalf <- which(ontogeny_outline[[i]][, 1] < 0)
#   ontogeny_star[i] <- lefthalf[which.min(abs(
#                                      ontogeny_outline[[i]][lefthalf, 2] - 0))]
#   plot(ontogeny_outline[[i]], asp = 1 , type = "l")
#   points(ontogeny_outline[[i]][ontogeny_star[i], 1],
#            ontogeny_outline[[i]][ontogeny_star[i], 2], col = 4)
# }
# ***** MANUAL ADJUSTMENT END *****

# ========== semi-landmark extraction ===============
ontogeny_landmark <- .extract_landmark(ontogeny_outline, p = 200,
                                       star = ontogeny_star)

# check landmarks
# for (i in 1:length(ontogeny_outline)){
#  plot(ontogeny_landmark[, ,i], asp = 1, type = "l")
#  points(ontogeny_landmark[, 1,i], ontogeny_landmark[, 2, i],
#         pch = 21, bg = 2, cex = 0.4)
# }

# ========== output semi-landmarks ===============
geomorph::writeland.tps(ontogeny_landmark,
                        "data-raw/nauplii-ontogeny-semilandmarks.tps")
devtools::use_data(ontogeny_landmark, overwrite = TRUE)
