#################################################################################
# Setting directory folder                                                   ----
#################################################################################
setwd("E:/OneDrive - Universidad de Valladolid/Doutorado/02-TUM")
data_in <- "01-Data/01-LiDAR"
SINGLE_TREE_832 <- file.path(data_in, "Single Tree/832")
SINGLE_TREE_832_OUT <- file.path(SINGLE_TREE_832, "01-Stem")
#dir.create(SINGLE_TREE_832_OUT)         # To create the new directory

SINGLE_TREE_833 <- file.path(data_in, "Single Tree/833")
SINGLE_TREE_833_OUT <- file.path(SINGLE_TREE_833, "01-Stem")
# dir.create(SINGLE_TREE_833_OUT)        # To create the new directory

#################################################################################
# Install / Load packages                                                    ----
#################################################################################
library(pacman)

pacman::p_load("rlas", "geometry", "rgl", "alphashape3d", "lidR",
               "magrittr", "raster", "dbscan", "VoxR", "TreeLS", "conicfit",
               "dplyr", "readr", "qdapTools", "rlang", "forecast", "sp", "ggplot2",
               "alphahull", "tripack", "tidyr", "purrr", "deldir", "wesanderson", "conicfit", 
               "beepr")

#################################################################################
# Reading the data                                                           ----
#################################################################################
LASfile <- list.files(SINGLE_TREE_833, pattern = "las", full.names = T)

#################################################################################
# Data tree function                                                          ----
#################################################################################
j <- 9                    # Just for testing

fts_stem <- function(data) {       # Enter the las file.list as data, and Z = the diameter height
  # Initializes the progress bar
  pb1 <- txtProgressBar(title = "Total Progress", 
                        min = 0,      # Minimum value of the progress bar
                        max = length(data), # Maximum value of the progress bar
                        style = 3,    # Progress bar style (also available style = 1 and style = 2)
                        width = 50,   # Progress bar width. Defaults to getOption("width")
                        char = "=")   # Character used to create the bar
  
  
  
  for (j in 1:length(data)) {
    # Reading and fixing da data
    .tree <- read.las(data[j])
    .tree <- .tree[,c("X", "Y", "Z")]
    .tree <- as.data.frame(.tree)
    .name <- sub(".*/", "", data[j])
    
    # cat(paste0("Processing file: ", .name))
    
    # Compressing the Z
    .tree$Zw <- .tree$Z * 0.05
    
    # Creating height class every 5 cm
    .tree$heightclass <- cut(.tree$Z,((max(.tree$Z)- min(.tree$Z))/0.5), labels=F) # Creating a 5 cm class interval
    
    # Clustering the data
    .stem <- matrix(ncol = 7)
    colnames(.stem) <- c("X", "Y", "Z", "heightclass", "dist", "Zw", "cluster")
    
    # STEM Analysis
    .epsi <-  0.03
    .tree <- na.omit(.tree)
    .tree <- as.data.table(.tree)
    .clean <- filter_noise(.tree, k = 20, sigma = 1, store_noise = F, message = F)
    
    # Initializes the progress bar
    pb2 <- txtProgressBar(title = "Single Tree Process:", 
                          min = 0,
                          max = length(unique(.clean$heightclass)),
                          style = 3,
                          width = 50,
                          char = "=")
    
    for(zzz in unique(.clean$heightclass)) {
      .file_stem_section <- .clean[.clean$heightclass == zzz, ]
      # print(paste0(zzz, " - ", nrow(.file_stem_section)))
      if (nrow(.file_stem_section) <= 50) {next} else {
        .file_stem_section <- as.data.table(.file_stem_section)
        .center <- dbscan(.file_stem_section[ ,c('X','Y', 'Zw')],                 # DBSCAN to cluster the data
                          eps = .epsi,
                          minPts = nrow(.file_stem_section)*0.01,
                          borderPoints = T)
        
        .file_stem_section$cluster <- .center$cluster                                          # Adding the cluster data do the tree
        .prestem <- .file_stem_section[.file_stem_section$cluster != 0, ]                                   # Filtering unclassified Points
        .stem2 <- .prestem[.prestem$cluster == as.integer(which.max(table(.prestem$cluster))), ] # Filtering the class with maximum value
        .stem <- rbind(.stem, .stem2)
        
        #plot3d(.stem2, col = .stem$cluster + 1, add = F, size = 1)
      }
      setTxtProgressBar(pb2, pb2$getVal()+1)
    }
    setTxtProgressBar(pb1, pb1$getVal()+1)
    write.csv(.stem, file = file.path(SINGLE_TREE_833_OUT, paste0(sub(".las", "", .name), ".csv")))
  }
  #return(.stem)
  beep(3)
}

fts_stem(LASfile)

plot3d(.tree, col = "black", add = T, size = 1); aspect3d("iso")
plot3d(.clean, col = "red", add = T, size = 1); aspect3d("iso")
plot3d(.stem, col = "blue", add = T, size = 2)
plot3d(.prestem, col = .prestem$cluster, add = F, size = 1); aspect3d("iso")
plot3d(.prestem, col = .prestem$cluster+1,add=T,size = 0.1); aspect3d("iso")

write.csv(stem, file = file.path(SINGLE_TREE_832_OUT, paste0(sub(".las", "", .name), ".csv")))
