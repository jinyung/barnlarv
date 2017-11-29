# =================== NOTES =============================
# data-raw/stage2outline.R was run prior to this
# =======================================================

# ========== required packages ===============
require(geomorph)
require(devtools)

# ========== call outlines for processing ===============
data(stage2outline)

# ========== Get the starting point for semi-landmark extraction ===============
# get the starting point
stage2star <- integer(length(stage2outline))
for (i in 1:length(stage2outline)) {
  stage2star[i] <- which.min(stage2outline[[i]][, 1])
}

# ***** NOTES: BELOW PARTS REQUIRE MANUAL ADJUSTMENT *****
# visually check for starting points
for (i in 1:length(stage2outline)) {
  plot(stage2outline[[i]], asp = 1, type = "l")
  points(stage2outline[[i]][stage2star[i], 1], stage2outline[[i]][stage2star[i],
                                                                  2], col =3)
  centroid <- apply(stage2outline[[i]], 2, mean)
  text(centroid[1], centroid[2], paste(i, names(stage2outline[i])))
}

# # stage2outline with wrong starting points
# # ***** manual adjustments commented out
# # following species starting point needed to be adjusted
# idx <- c("Sacculina pilosella", "Sacculina carcini",
#          "Polyascus polygenea", "Polyascus plana", "Heterosaccus papillosus",
#          "Heterosaccus lunatus", "Briarosaccus callosus")
# spidx <- which(names(stage2outline) %in% idx)

# # correction: define start point as left half closest to the x-axis
# # plot to check again
# for (i in spidx){
#   lefthalf <- which(stage2outline[[i]][, 1] < 0)
#   stage2star[i] <- lefthalf[which.min(abs(stage2outline[[i]][lefthalf, 2] -
#                                            0))]
#   plot(stage2outline[[i]], asp = 1 , type = "l")
#   points(stage2outline[[i]][stage2star[i], 1],
#            stage2outline[[i]][stage2star[i], 2], col = 4)
#   text(centroid[1], centroid[2], paste(i, names(stage2outline[i])))
# }

# # Polyascus plana still has problem
# spidx <- "Polyascus plana"
# for (i in which(names(stage2outline) %in% spidx)) {
#   plot(stage2outline[[i]], asp = 1, type = "l")
#   point <- unlist(locator(1))
#   stage2star[i] <- which.min(apply(stage2outline[[i]], 1,
#                                    .ed, point))
#                 # stage2outline point nearest to my click
# }
# ***** MANUAL ADJUSTMENT END *****

# ========== semi-landmark extraction ===============
stage2landmark <- .extract_landmark(stage2outline, p = 200, star = stage2star)

# # visually check landmarks
# for (i in 1:length(stage2outline)){
#  plot(stage2landmark[, ,i], asp = 1, type = "l")
#  points(stage2landmark[, 1, i], stage2landmark[, 2, i], pch = 21,
#         bg = 2, cex = 0.4)
# }

# ========== output semi-landmarks ===============
# save to data-raw/ or data/
geomorph::writeland.tps(stage2landmark,
                        "data-raw/nauplii-stage2-semilandmarks.tps")
devtools::use_data(stage2landmark, overwrite = TRUE)
