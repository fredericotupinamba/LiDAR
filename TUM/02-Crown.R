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
               "alphahull", "tripack", "tidyr", "purrr", "deldir", "wesanderson", "conicfit",
               "stringr")

#################################################################################
# Reading the data                                                           ----
#################################################################################
LASfile <- list.files(SINGLE_TREE_832, pattern = "las", full.names = T)

#################################################################################
# Reading the fiel data                                                      ----
#################################################################################
FIELDDATA_Increment <- readxl::read_excel(file.path("E:/OneDrive - Universidad de Valladolid/Doutorado/02-TUM/01-Data/03-Field Work", "FieldData_Final.xlsx" ),
                                         sheet = "Table1")

FIELDDATA_Diameter <- readxl::read_excel(file.path("E:/OneDrive - Universidad de Valladolid/Doutorado/02-TUM/01-Data/03-Field Work", "FieldData_Final.xlsx" ),
                                         sheet = "Table2")

FIELDdata <- FIELDDATA_Increment %>% select(Site, TreeID, Type, DBH, H = "H (tape)", CBH)

#################################################################################
# Data tree function                                                          ----
#################################################################################
j <- LASfile[1]

fts_crown <- function(data) {                                         # Enter the las file.list as data, and Z = the diameter height
  for (j in data) {
    # Reading and fixing da data
    .tree <- read.las(j)
    .tree <- .tree[,c("X", "Y", "Z")]
    .tree <- as.data.frame(.tree)
    .name <- sub(".*/", "", j)
    .Site <- substr(.name, 1, 3)
    .TreeID <- str_replace(.name, pattern = ".*_([^-]*)_.*", replacement = "\\1")
    CBH <- FIELDdata %>% filter(Site == .Site, TreeID == .TreeID) %>% select(CBH) %>% as.numeric()
    
    # Cutting the Crown
    .tree <- .tree[.tree$Z >= (min(.tree$Z) + CBH), ]
    
    # Creating height class every 5 cm
    .tree$heightclass <- cut(.tree$Z,((max(.tree$Z)- min(.tree$Z))/0.01), labels=F) # Creating a 5 cm class interval
    .tree$Zw <- .tree$Z * 0.05
    
    # STEM Analysis
    .epsi <-  0.05
    .tree <- na.omit(.tree)
    .center <- dbscan(.tree[ ,c('X', 'Y', "Zw")],  # DBSCAN to cluster the data
                      eps = .epsi,
                      minPts = 80,
                      borderPoints = T) 
    
    .tree$cluster <- .center$cluster                                          # Adding the cluster data do the tree
    .precrown <- .tree[.tree$cluster != 0, ]                                   # Filtering unclassified Points
    .crown <- .precrown[.precrown$cluster == as.integer(which.max(table(.precrown$cluster))), ] # Filtering the class with maximum value
    
    plot3d(.tree, col = "black", add = F, size = 1); aspect3d("iso")
    plot3d(.crown, col = .center$cluster + 1, add = T, size = 2)
    
  }
  return(.crown)
}

crown <- fts_crown(LASfile[1])

plot3d(.tree, col = "black", add = F, size = 1); aspect3d("iso")
plot3d(crown, col = .tree$cluster, add = T, size = 2)
