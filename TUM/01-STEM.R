#################################################################################
# Setting directory folder                                                   ----
#################################################################################
#setwd("Z:/T-LIDAR/TLS_Bonhorst_Anpassung_Waldbaum/Frederico")
setwd("E:/OneDrive - Universidad de Valladolid/Doutorado/02-TUM")
data_in <- "01-Data/01-LiDAR"
SINGLE_TREE_832 <- file.path(data_in, "Single Tree/832")
SINGLE_TREE_832_OUT <- file.path(SINGLE_TREE_832, "01-Stem")
dir.create(SINGLE_TREE_832_OUT)

SINGLE_TREE_833 <- file.path(data_in, "Single Tree/833")
SINGLE_TREE_833_OUT <- file.path(SINGLE_TREE_833, "01-Stem")
dir.create(SINGLE_TREE_833_OUT)

#################################################################################
# Install / Load packages                                                    ----
#################################################################################
library(pacman)

pacman::p_load("rlas", "geometry", "rgl", "alphashape3d", "lidR",
               "magrittr", "raster", "dbscan", "VoxR", "TreeLS", "conicfit",
               "dplyr", "readr", "qdapTools", "rlang", "forecast", "sp", "ggplot2",
               "alphahull", "tripack", "tidyr", "purrr", "deldir", "wesanderson", "conicfit")

#################################################################################
# Reading the data                                                           ----
#################################################################################
LASfile <- list.files(SINGLE_TREE_832, pattern = "las", full.names = T)

#################################################################################
# Data tree function                                                          ----
#################################################################################
j <- LASfile[1]

fts_stem <- function(data) {                                         # Enter the las file.list as data, and Z = the diameter height
  for (j in data) {
    # Reading and fixing da data
    .tree <- read.las(j)
    .tree <- .tree[,c("X", "Y", "Z")]
    .tree <- as.data.frame(.tree)
    .name <- sub(".*/", "", j)
    
    # Creating height class every 5 cm
    .tree$heightclass <- cut(.tree$Z,((max(.tree$Z)- min(.tree$Z))/0.05), labels=F) # Creating a 5 cm class interval
    
    # Clustering the data
    .stem <- matrix(ncol = 6)
    colnames(.stem) <- c("X", "Y", "Z", "heightclass", "dist", "cluster")
    
    # STEM Analysis
    .epsi <-  0.01
    .tree <- na.omit(.tree)
    for(zzz in unique(.tree$heightclass)) {
      .file_stem_section <- .tree[.tree$heightclass == zzz, ]
      #print(paste0(zzz, " - ", nrow(.file_stem_section)))
      if (nrow(.file_stem_section) <= 50) {next} else {
      .file_stem_section <- as.data.table(.file_stem_section)
      .file_stem_section <- filter_noise(.file_stem_section, k = 5, message = F)
      .center <- dbscan(.file_stem_section[ ,c('X','Y')],                 # DBSCAN to cluster the data
                        eps = .epsi,
                        minPts = nrow(.file_stem_section)*0.0025,
                        borderPoints = T) 
      
      .file_stem_section$cluster <- .center$cluster                                          # Adding the cluster data do the tree
      .prestem <- .file_stem_section[.file_stem_section$cluster != 0, ]                                   # Filtering unclassified Points
      .stem2 <- .prestem[.prestem$cluster == as.integer(which.max(table(.prestem$cluster))), ] # Filtering the class with maximum value
      .stem <- rbind(.stem, .stem2)
      }
    }
  }
  return(.stem)
}

stem <- fts_stem(LASfile[1])

plot3d(.tree, col = "black", add = F, size = 1); aspect3d("iso")
plot3d(stem, col = 'blue', add = T, size = 2)

write.csv(stem, file = file.path(SINGLE_TREE_832_OUT, paste0(sub(".las", "", .name), ".csv")))
