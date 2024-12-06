rm(list=ls())

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggplot2)
library(ggpubr)
library(viridisLite)
library(ggplot2)
library(tidyr)
library("tidyverse")
library("cowplot")
library("svglite")
library("pheatmap")



## This is the code I used to make it functional, above is trial and error

setwd("/Users/user/Documents/02_Andropogon/03_SPIPS_Metagenome/Hays_IL_mags")

data <- read.table("Mags_detection.csv", sep=",", header=TRUE, row.names=1)
data = data[,3:26]

#pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE)
metadata <- read.table("Mags_detection_metadata.csv", sep=",", header=TRUE, row.names=1)
pc = metadata[,2:3]
pc2 = metadata[,1:3]
pc3 = metadata[3]

pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=metadata, fontsize_col = 4, fontsize_row = 4)


svglite("HeatMap_Mags.svg")



my_colors <- c("#86c306", "#fc62e1","#d87fff", "#1bd3d8", "#d49900", "#1cd570")

# Create the heatmap
plot <- pheatmap(data,
                 cluster_rows = FALSE,
                 show_rownames = TRUE,
                 cluster_cols = TRUE,
                 annotation_col = pc2,
                 fontsize_col = 8,
                 fontsize_row = 10)

annotation_colors <- list(
  Ecotype = c("Dry" = "#d49900", "Wet" = "#1cd570"),  # Example for first column
  Population = c("CDB" = "#86c306", "WEB" = "#fc62e1", "DES" = "#d87fff","FUL" = "#1bd3d8")  # Example for second column
)


pheatmap(data,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         annotation_col = pc,
         fontsize_col = 8,
         fontsize_row = 10,
         annotation_colors = annotation_colors)



plot <-pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=pc2, fontsize_col = 8, fontsize_row = 10)

metadata2 <- read.table("Mags_detection_taxa.csv", sep=",", header=TRUE, row.names=1)
plot <- pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=metadata, annotation_row=metadata2, fontsize_col = 4, fontsize_row = 4)


svglite("HeatMap_Mags.svg")

