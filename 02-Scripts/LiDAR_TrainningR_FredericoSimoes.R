# TITLE:   SPATIAL ANALYSIS (MEDFOR+DATAFOREST) (META49405-2022)
#
# DATE:    2022-11-07
#
# AUTOR:  Frederico Tupinambá Simões
#         PhD Candidate at the University of Valladolid 
#         Early Stage Researcher at Skill-For-Action Project
#
# E-MAIL: frederico.tupinamba@uva.es
#
# Study Area: Jarandilla Marteloscopes
#
# Reference: https://r-lidar.github.io/lidRbook/index.html

#=============================================================================================================#
# Installing the PACMAN package if necessary:                                                              ----
#=============================================================================================================#
# install.packages("pacman")

#=============================================================================================================#
# Loaging the packages:                                                                                    ----
#=============================================================================================================#
library(pacman)
pacman::p_load(lidR, ggplot2, sf, sp, maptools, tmap, mapview, future, dplyr, rgl)

# Setting the number of cores
plan(multisession, workers = 4L)

set_lidr_threads(4L)

#=============================================================================================================#
# Functions:                                                                                               ----
#=============================================================================================================#
# Creating cros-section plot function
plot_crossection <- function(las,
                             p1 = c(min(las@data$X), mean(las@data$Y)),
                             p2 = c(max(las@data$X), mean(las@data$Y)),
                             width = 4, colour_by = NULL)
{
  colour_by <- enquo(colour_by)
  data_clip <- clip_transect(las, p1, p2, width)
  p <- ggplot(data_clip@data, aes(X,Z)) + geom_point(size = 0.5) + coord_equal() + theme_minimal()
  
  if (!is.null(colour_by))
    p <- p + aes(color = !!colour_by) + labs(color = "")
  
  return(p)
}

# Fill NAs and SMOOTH a CHM
fill.na <- function(x, i = 5) {if (is.na(x)[i]) {return(mean(x, na.rm = T))} else { return(x[i]) }}

# Local Maximum Filter with variable windows size
f <- function(x) {
  y <- 2.6 * (-(exp(-0.08*(x-2)) - 1)) + 3
  y[x < 2] <- 3
  y[x > 20] <- 5
  return(y)
}


#=============================================================================================================#
# Setting the LAS folder:                                                                                 -----
#=============================================================================================================#
# Here you must change only the main directory, where the downloaded file was saved.
setwd("E:/OneDrive - Universidad de Valladolid/Doutorado/GITHUB/LiDAR/LiDAR-Uva/LiDAR")
standDir <- file.path("01-DATA/01-Stand")
shpDir <- file.path("01-DATA/03-Boundary")


#####################################
###         STOP NOW              ###
# Thin the data at Cloud Compare    #
###                               ###
#####################################

#=============================================================================================================#
# Importing LAS file:  ++++                                                                                 
#=============================================================================================================#
lasLIST <- list.files(standDir, pattern = "*.las$", full.names = F)
las <- lidR::readLAS(file.path(standDir, lasLIST[1]))

# Setting the projection
crs(las) <- 25830

# Getting information about the point cloud
summary(las)
las_check(las)

#=============================================================================================================#
# Clipping the point cloud                                                                                 ----
#=============================================================================================================#
# Importing the Marteloscope (100 x 100 meters)
POI <- st_read(file.path(shpDir, "Quadrantes.shp"))

POI <- st_transform(POI, 25830) # Setting a UTM F30 projectionc
crs(POI)

# Printing the map
tmap_mode("view")                                                    # plot or view

tm_shape(POI, map = TRUE, map.type = "Esri.WorldImagery") + 
  tm_polygons(alpha = 0, lwd = 2, border.col = "black") +
  tm_legend(legend.outside = TRUE,
            legend.outside.position = "right",
            position = c("LEFT", "BOTTOM"),
            frame = F)

# Selecting only one small square (25 x 25 meters)
C3 <- subset(POI, Plot_ID == 8)                                  # Filtering a small area

# Clipping the Point Cloud with the small square
las_CLIPPED <- clip_roi(las, C3)

# Saving the LAS data
writeLAS(las_CLIPPED, file.path(standDir, paste0("01-Clipped_", las[1])))

# Saving memory
remove(las_01, las_CLIPPED, POI, C3)

#=============================================================================================================#
# Adjusting and checking some information in the clipped point cloud
#=============================================================================================================#
lasLIST <- list.files(standDir, pattern = "*.las$", full.names = F)
las_CLIPPED <- lidR::readLAS(file.path(standDir, lasLIST[1]))

# Checking for inconsistency in the data
summary(las_CLIPPED)
las_check(las_CLIPPED)

# Removing the duplicated points 
las_CLIPPED <- filter_duplicates(las_CLIPPED)

# Noise classification
las_02A <- classify_noise(las_CLIPPED, algorithm = sor(k = 10, m = 3, quantile = F))
las_02B <- classify_noise(las_CLIPPED, algorithm = ivf(res = 0.5, n = 50))
las_02C <- classify_noise(las_CLIPPED, algorithm = ivf(5, 2))

# Plotting the point cloud
plot(las_02A, bg = "white", axis = T, legend = T, color = "Classification")       # Commented because it takes some time do run!
DATA <- las_02B@data %>% as.data.frame() %>% select(X, Y, Z, Classification)
plot3d(DATA, zlim = min(DATA$Z), type = "p", col = DATA$Classification + 1, aspect = "iso")


# Saving the LAS Data
writeLAS(las_02A, file.path(noiseDir, paste0("01-NC_", sub("01-Clipped_", "", las[1]))))
writeLAS(las_02B, file.path(noiseDir, paste0("02-NC_", sub("01-Clipped_", "", las[1]))))
writeLAS(las_02C, file.path(noiseDir, paste0("03-NC_", sub("01-Clipped_", "", las[1]))))

# Saving memory
remove(las, las_02, las_02A, las_02B, las_02C)

#=============================================================================================================#
# Ground classification                                                                                    ----
#=============================================================================================================#

#=============================================================================================================#
# Working with LAS Catalog:                                                                                ----
#=============================================================================================================#
# Loading the CLIP file
ctgName <- list.files(file.path(noiseDir), pattern = "*.las")
ctg <- readLAScatalog(file.path(noiseDir, ctgName[2]))

# Optimization parameters
opt_chunk_size(ctg) <- 3
opt_chunk_buffer(ctg) <- 0
opt_laz_compression(ctg) <- T
opt_merge(ctg) <- T
opt_progress <- T
opt_stop_early(ctg) <- F
#opt_chunk_alignment(ctg) <- c(25, -4)
plot(ctg, chunk = T, map = T)

## Progressive Morphological Filter 
# PMF classification
myfunPMF = function(cluster, ws, th){
  las <- lidR::readLAS(cluster)
  if (is.empty(las)) return(NULL)
  las <- classify_ground(las, pmf(ws, th))
  return(las)
}

# Force the results to be written on disk
opt_output_files(ctg) <- file.path(groundDir, "01-Ground-PMF_{XLEFT}_{YBOTTOM}")

opt <- list(need_buffer = TRUE, automerge = TRUE)

out1 <- catalog_apply(ctg, myfunPMF, ws=1, th=1, .options = opt)

plot(out)

# Plotting the classification
plot(las_03A, color = "Classification", size = 3, bg = "white", axis = TRUE, legend = TRUE)


## Cloth Simulation Function
# CSF classification
myfunCSF = function(cluster){
  las <- lidR::readLAS(cluster)
  if (is.empty(las)) return(NULL)
  las <- classify_ground(las, csf(sloop_smooth = TRUE, class_threshold = 1, cloth_resolution = 1, time_step = 1), last_returns = TRUE)
  
  return(las)
}

# Force the results to be written on disk
opt_output_files(ctg) <- file.path(groundDir, "02-Ground-CSF_{XLEFT}_{YBOTTOM}")

opt <- list(need_buffer = TRUE, automerge = TRUE)

out2 <- catalog_apply(ctg, myfunCSF, .options = opt)

plot(out2)

# Saving memory
remove(las, las_02, las_02A, las_02B, las_02C, ctg, newctg, opt, out, out1, out2, newout)

#=============================================================================================================#
# Working with LAS file                                                                                    ----
#=============================================================================================================#
lasCG1 <- classify_ground(las_02B, algorithm = pmf(ws = 5, th = 3))
lasCG2 <- classify_ground(las_02B, myfunPMF, ws=1, th=1)

ws <- seq(3, 12, 3)
th <- seq(0.1, 1.5, length.out = length(ws))
lasCG3 <- classify_ground(las_02B, algorithm = pmf(ws = ws, th = th))

mycsf <- csf(sloop_smooth = TRUE, class_threshold = 1, cloth_resolution = 1, time_step = 1)
lasCG4 <- classify_ground(las_02B, mycsf)

plot(lasCG1, color = "Classification", size = 3, bg = "white") 

p1 <- c(276472, 4440347)
p2 <- c(276505, 4440380)
plot_crossection(lasCG4, p1 , p2, colour_by = factor(Classification))


gnd <- filter_ground(lasCG3)
plot(gnd, size = 3, bg = "white") 

#=============================================================================================================#
# Digital terrain model:                                                                                   ----
#=============================================================================================================#
las <- list.files(groundDir, pattern = " *.laz$", full.names = F)
las_04 <- lidR::readLAS(file.path(groundDir, las[1]))

## 1 - Generete a DTM model with the TIM algorithm
DTM_TIN <- rasterize_terrain(lasCG3, res = 1, algorithm = tin())
plot_dtm3d(DTM_TIN, bg = "white")

## 2 -  Invert Distance weighting (IDW)
DTM_IDW <- rasterize_terrain(lasCG3, algorithm = knnidw(k = 10L, p = 2))
plot_dtm3d(DTM_IDW, bg = "white")

## 3 - Kriging
DTM_KRI <- grid_terrain(lasCG3, algorithm = kriging(k = 40))
plot_dtm3d(DTM_KRI, bg = "white")

# Writing the raster
writeRaster(DTM_TIN, file.path(dtmDir, paste0(sub("*.laz", "", las[1]), "_DTM_TIN.tif")))
writeRaster(DTM_IDW, file.path(dtmDir, paste0(sub("*.laz", "", las[1]), "_DTM_IDW.tif")))
writeRaster(DTM_KRI, file.path(dtmDir, paste0(sub("*.laz", "", las[1]), "_DTM_KRI.tif")))

# Saving memory
remove(las, las_04, DTM_KRI, DTM_IDW)


#=============================================================================================================#
# Height normalization:                                                                                    ----
#=============================================================================================================#
las <- list.files(groundDir, pattern = "*.laz$", full.names = F)
las_05 <- lidR::readLAS(file.path(groundDir, las[1]))

# Point cloud normalization
NLAS <- normalize_height(lasCG3, knnidw())

NLAS2 <- normalize_height(lasCG3, tin(), dtm = DTM_TIN)

las_check(NLAS)
las_check(NLAS2)

hist(filter_ground(NLAS)$Z, breaks = seq(-1, 1, 0.01), main = "", xlab = "Elevation")
hist(filter_ground(NLAS2)$Z, breaks = seq(-2, 2, 0.01), main = "", xlab = "Elevation")


# Saving the LAS Data
writeLAS(NLAS, file.path(nlasDir, paste0(sub("*.laz", "", las[1]), "_NLAS.laz")))

# Saving memory
remove(las, las_05, NLAS, NLAS2)


#=============================================================================================================#
# Canopy Height model:                                                          -----
#=============================================================================================================#
las <- list.files(nlasDir, pattern = " *.laz$", full.names = F)
las_06 <- lidR::readLAS(file.path(nlasDir, las[1]))

## Point to raster
# A few different examples
chm <- grid_canopy(las_06, res = 1, algorithm = p2r()) # 1 meter resolution
col <- height.colors(50)
plot(chm, col = col)

chm <- grid_canopy(las_06, res = 0.5, algorithm = p2r()) # 0.5 meter resolution
plot(chm, col = col)

chm <- grid_canopy(las_06, res = 0.5, algorithm = p2r(subcircle = 0.15)) # 0.5 meter resolution, 0.15 cm the point size
plot(chm, col = col)

## Triangulation
chm <- grid_canopy(las_06, res = 0.5, algorithm = dsmtin(max_edge = 8))
plot(chm, col = col)

## Pit-free algorithm
chm <- grid_canopy(las_06, res = 0.5, pitfree(thresholds = c(0, 2, 5), max_edge = c(0, 1.5)))
plot(chm, col = col)

chm <- grid_canopy(las_06, res = 0.5, 
                   pitfree(thresholds = c(0, 2, 5), 
                           max_edge = c(0, 1.5),
                           subcircle = 0.15))
plot(chm, col = col)

## Post-processing a CHM
w <- matrix(1, 3, 3)

filled <- focal(chm, w, fun = fill.na)
smoothed <- focal(chm, w, fun = mean, na.rm = T)

chms <- stack(chm, filled, smoothed)
names(chms) <- c("Base", "Filled", "Smoothed")
plot(chms, col = col)


#=============================================================================================================#
# Individual tree detection                                                                                ----
#=============================================================================================================#
las <- list.files(nlasDir, pattern = " *.laz$", full.names = F)
las_06 <- lidR::readLAS(file.path(nlasDir, las[1]))

# Local Maximum Filter with fixed windows size
ttops_5 <- locate_trees(las_06, lmf(ws = 5))
ttops_3 <- locate_trees(las_06, lmf(ws = 3))

## Local Maximum Filter with variable windows size
ttops_F <- find_trees(las_06, lmf(f))

# Plotting the trees Top
x <- plot(las_06, bg = "white", size = 2)
add_treetops3d(x, ttops_F)

par(mfrow=c(1,3))
plot(chm, col = height.colors(50))
plot(ttops_5, add = T)
plot(chm, col = height.colors(50))
plot(ttops_3, add = T)
plot(chm, col = height.colors(50))
plot(ttops_F, add = T, size = 2, col = "BLACK", symbols = 5)


#=============================================================================================================#
# Local Maximum Filter on a CHM                                                                             ----
#=============================================================================================================#
# Point-to-raster 2 resolutions
chm_p2r_05 <- grid_canopy(las_06, 0.5, p2r(subcircle = 0.2))
chm_p2r_1 <- grid_canopy(las_06, 1, p2r(subcircle = 0.2))

# Pitfree with and without subcircle tweak
chm_pitfree_05_1 <- grid_canopy(las_06, 0.5, pitfree())
chm_pitfree_05_2 <- grid_canopy(las_06, 0.5, pitfree(subcircle = 0.2))

# Post-processing median filter
kernel <- matrix(1,3,3)
chm_p2r_05_smoothed <- raster::focal(chm_p2r_05, w = kernel, fun = median, na.rm = TRUE)
chm_p2r_1_smoothed <- raster::focal(chm_p2r_1, w = kernel, fun = median, na.rm = TRUE)

# Finding trees with chm
ttops_chm_p2r_05 <- find_trees(chm_p2r_05, lmf(5))
ttops_chm_p2r_1 <- find_trees(chm_p2r_1, lmf(5))
ttops_chm_pitfree_05_1 <- find_trees(chm_pitfree_05_1, lmf(5))
ttops_chm_pitfree_05_2 <- find_trees(chm_pitfree_05_2, lmf(5))
ttops_chm_p2r_05_smoothed <- find_trees(chm_p2r_05_smoothed, lmf(5))
ttops_chm_p2r_1_smoothed <- find_trees(chm_p2r_1_smoothed, lmf(5))

# Plotting the tree
par(mfrow=c(3,2))
col <- height.colors(50)
plot(chm_p2r_05, main = "CHM P2R 0.5", col = col); plot(ttops_chm_p2r_05, add = T)
plot(chm_p2r_1, main = "CHM P2R 1", col = col); plot(ttops_chm_p2r_1, add = T)
plot(chm_p2r_05_smoothed, main = "CHM P2R 0.5 smoothed", col = col); plot(ttops_chm_p2r_05_smoothed, add = T)
plot(chm_p2r_1_smoothed, main = "CHM P2R 1 smoothed", col = col); plot(ttops_chm_p2r_1_smoothed, add = T)
plot(chm_pitfree_05_1, main = "CHM PITFREE 1", col = col); plot(ttops_chm_pitfree_05_1, add = T)
plot(chm_pitfree_05_2, main = "CHM PITFREE 2", col = col); plot(ttops_chm_pitfree_05_2, add = T)


#=============================================================================================================#
# Individual tree detection                                                                                ----
#=============================================================================================================#
algo <- dalponte2016(chm_p2r_05_smoothed, ttops_chm_p2r_05_smoothed)
las <- segment_trees(las_06, algo)

# Plotting the segmeted trees point cloud
plot(las, bg = "white", size = 4, color = "treeID", pal = forest.colors(50))

# Plotting a single tree
tree <- filter_poi(las, treeID == 5)
plot(tree, size = 8, bg = "white")

#=============================================================================================================#
# Tree metrics                                                                                            ----
#=============================================================================================================#
metrics <- crown_metrics(las, ~list(z_max = max(Z), z_mean = mean(Z))) # calculate tree metrics
head(metrics)

# user-defined function
custom_crown_metrics <- function(z, i) { 
  metrics <- list(
    z_max = max(z),   # max height
    z_sd = sd(z),     # vertical variability of points
    i_mean = mean(i), # mean intensity
    i_max  = max(i)   # max intensity
  )
  return(metrics) # output
}

# Definning the function
ccm = ~custom_crown_metrics(z = Z, i = Intensity)

# Running the function
metrics <- crown_metrics(las, func = ccm, geom = "convex")

# Plotting the tree maximum height
plot(metrics["z_max"], pal = hcl.colors)

# Running the function
metrics <- crown_metrics(las, func = ccm, geom = "concave")
plot(metrics["z_max"], pal = hcl.colors)




tree.list.tls <- tree.detection.multiple(las.list = files,
                                         
                                         normalize.arguments = list(max.dist = 10,
                                                                    min.height = 0.25,
                                                                    max.height = 25,
                                                                    algorithm.dtm = "knnidw",
                                                                    res.dtm = 0.5),
                                         
                                         tree.detection.arguments = list(dbh.min = 7.5, dbh.max = 100,
                                                                         tls.resolution = list(point.dist = 7.67,
                                                                                               tls.dist = 10),
                                                                         breaks = 1.3),
                                         
                                         dir.data = dir.data, save.result = TRUE, dir.result = dir.result)