#################################################################################
# Setting directory folder                                                   ----
#################################################################################
setwd("E:/OneDrive - Universidad de Valladolid/Doutorado/02-TUM")
data_in <- "01-Data/01-LiDAR"

# Plot 832
SINGLE_TREE_832 <- file.path(data_in, "Single Tree/832")
STEM_TREE_832 <- file.path(data_in, "Single Tree/832/01-Stem")
D_TREE_832 <- file.path(data_in, "Single Tree/832/03-Diameter")
#dir.create(D_OUT)         # To create the new directory

#Plot 833
SINGLE_TREE_833 <- file.path(data_in, "Single Tree/833")
STEM_TREE_833 <- file.path(SINGLE_TREE_833, "01-Stem")
D_TREE_833 <- file.path(data_in, "Single Tree/833/03-Diameter")
#dir.create(D_TREE_833)        # To create the new directory

#################################################################################
# Install / Load packages                                                    ----
#################################################################################
library(pacman)

pacman::p_load("rlas", "geometry", "rgl", "alphashape3d", "lidR",
               "magrittr", "raster", "dbscan", "VoxR", "TreeLS", "conicfit",
               "dplyr", "readr", "qdapTools", "rlang", "forecast", "sp", "ggplot2",
               "alphahull", "tripack", "tidyr", "purrr", "deldir", "wesanderson", "conicfit",
               "stringr", "beepr")

#################################################################################
# Listing the data                                                           ----
#################################################################################
direcao <- SINGLE_TREE_832
i = 3
j = 33


fts_diameter <- function(direcao){
  .dados <- list.files(file.path(direcao, "01-Stem"), pattern = ".csv", full.names = T)
  
  # Creating the matrix for the results
  tapper <- matrix(ncol = 7)
  .tapper1 <- matrix(ncol = 7)
  .tapper2 <- matrix(ncol = 7)
  colnames(tapper) <- c("Tree", "Height", "X", "Y", "Z", "Cir_diam", "Ellip_diam")
  colnames(.tapper1) <- c("Tree", "Height", "X", "Y", "Z", "Cir_diam", "Ellip_diam")
  colnames(.tapper2) <- c("Tree", "Height", "X", "Y", "Z", "Cir_diam", "Ellip_diam")
  
  # Initializes the progress bar
  pb1 <- txtProgressBar(title = "Total Progress", 
                        min = 0,      # Minimum value of the progress bar
                        max = length(.dados), # Maximum value of the progress bar
                        style = 3,    # Progress bar style (also available style = 1 and style = 2)
                        width = 50,   # Progress bar width. Defaults to getOption("width")
                        char = "=")   # Character used to create the bar
  
  for (i in 1:length(.dados)) {
    .tree <- read.csv(.dados[i])
    .tree <- .tree[,c("X", "Y", "Z", "cluster")]
    .tree <- .tree %>% drop_na()
    .name <- sub(".*/", "", .dados[i]) %>% sub(".csv", "", .)
    
    # Tree height and tapper
    .Z <- max(.tree$Z) - min(.tree$Z)
    .tapper <- c(0.10, 0.5, 1, 1.3, 1.5)   # Initial diameter
    .incremento = 0.5
    while (.Z >= max(.tapper) + .incremento) {
      .tapper <- c(.tapper, max(.tapper) + .incremento)
    }
    
    for (j in 1:length(.tapper)) {
      # Stem classification
      .tree_diameter <- .tree[.tree$Z > (min(.tree$Z) + (.tapper[j] - 0.05)) & .tree$Z < (min(.tree$Z) + (.tapper[j] + 0.05)), ] # 10 cm from each tapper value 
      .tree_diameter <- .tree_diameter[!duplicated(.tree_diameter), ]
      .tree_diameter <- as.data.table(.tree_diameter)
      
      # Graphics
      png(file.path(direcao, "03-Diameter", paste0(.name, "_", .tapper[j], ".png")), width = 15, height = 15, units = 'in', res = 300)
      print(paste0(i, " - ", j," - N = ", nrow(.tree_diameter)))
      
      if (nrow(.tree_diameter) <= 300) {break} else { 
        .tree_diameter <- VoxR::distance_clustering(.tree_diameter, d_clust = 0.01, message = F)
        
        #plot(.tree_diameter[,1:2], type = "p", col = .tree_diameter$cluster)
        
        .tree_diameter2 <- .tree_diameter[.tree_diameter$cluster == names(sort(table(.tree_diameter$cluster), decreasing = TRUE)[1]), ]
        
        plot(.tree_diameter2[,1:2], type = "p", col = .tree_diameter2$cluster)
        
        ########################################################################
        # Getting the center point of the Stem
        .dd <- deldir(xyz.coords(.tree_diameter2))
        .cen.dd <- with(.dd$summary, sapply(list(x, y), weighted.mean, del.wts))
        .cen.dd[3] <- min(.tree$Z)+ .tapper[i]
        .center <- as.data.table(t(.cen.dd))
        colnames(.center) <- c("X", "Y", "Z")
        
        points(.center, col = 'red', pch = 16)
        
        
        ########################################################################
        # DBH measurements 
        .tree_diameter3 <- .tree_diameter2[,1:2] %>% distinct()
        
        # Calculating the Circle
        .Circle <- CircleFitBySpath(.tree_diameter3)
        .xyCircle <- calculateCircle(.Circle[1], .Circle[2], .Circle[3])
        .circle_pts <- data.frame(.xyCircle)
        .circle_pts$Z <- min(.tree$Z) + .tapper[i]
        .circle_diameter <- .Circle[3] * 200
        
        lines(.xyCircle[,1], .xyCircle[,2], col='blue')
        
        # Calculating the Ellipse
        .ellipDirect <- EllipseDirectFit(.tree_diameter3)
        .ellipDirectG <- AtoG(.ellipDirect)$ParG
        .xyDirect <- calculateEllipse(.ellipDirectG[1], .ellipDirectG[2], .ellipDirectG[3], 
                                     .ellipDirectG[4], 180/pi*.ellipDirectG[5])
        .xyDirect <- .xyDirect %>% as.data.frame() %>% round(digits = 3) %>% distinct()
        .ellip_pts <- ahull(.xyDirect, alpha = 0.5)
        .ellipDiameter <- (.ellip_pts$length / pi) * 100
        
        lines(.xyDirect[,1], .xyDirect[,2], col= 'red')
        dev.off()
      }
      .tapper1 <- c(Tree = .name, Height = .tapper[j], .center, Cir_diam = .circle_diameter, Ellip_diam = .ellipDiameter) %>% as.data.frame()
      .tapper2 <- rbind(.tapper2, .tapper1)
    }

    
    write.csv(.tapper2, file = file.path(direcao, "03-Diameter", paste0(.name, "_", .tapper[j], ".csv")))
    
    setTxtProgressBar(pb1, pb1$getVal()+1)
  }
  #return(tapper)
  beep(3)
}
  

fts_diameter(SINGLE_TREE_832)
