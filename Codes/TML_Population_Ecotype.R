#############################Clean up###########################################
rm(list=ls())

#############################Load libraries#####################################
library(lmerTest)
library(lme4)
library(nlme)
library(vegan)
library(labdsv)
library(MASS)
library(ggplot2)
library(car)
library(lme4)
library(lmerTest)
library(Rmisc)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(broom)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(viridisLite)
library(tidyr)
library("cowplot")
library("svglite")
library("pheatmap")


#######################Load Data - TML 5 MAGS ANOVA & Figures#################

setwd("/Users/user/Documents/02_Andropogon/03_SPIPS_Metagenome/00_DIABLO/5MAGS")
community<-read.csv(file="TML_Eco.csv", header = TRUE, sep = ",")

factors<-read.csv(file="5MAGs_Metadata_Eco.csv", header = TRUE, sep = ",")
as.factor (factors$Location)
as.factor (factors$Ecotype)
as.factor (factors$Population)
as.factor (factors$Row)

df = full_join(x = community, y = factors, by = "Group")
write.csv(df,"Merge.csv", row.names = TRUE)


df$Population <- ordered(df$Population,
                         levels = c("CDB","WEB", "FUL","DES"))



p <- ggplot(df, aes(x=Population, y=Trimethyllysine, color=Population), add = "jitter") + geom_boxplot(width=0.2) + theme_classic() +labs(x = "Population", colour = "Ecotype", y = "Trimethyllysine")

p + scale_x_discrete(limits=c("CDB","WEB", "DES", "FUL")) + scale_color_manual(values=c( "#E69F00","#E69F00", "#56B4E9", "#56B4E9"))



p <- ggplot(df, aes(x=Ecotype, y=Trimethyllysine, color=Ecotype), add = "jitter") + geom_boxplot(width=0.2) + theme_classic() +labs(x = "Ecotype", colour = "Ecotype", y = "Trimethyllysine")

p
two.way <- aov(Trimethyllysine ~ Ecotype+ Population +Ecotype*Population , data = df)
summary(two.way)
kruskal.test(Trimethyllysine ~ Ecotype, data = df)

kruskal.test(Trimethyllysine ~ Population, data = df)
pairwise.wilcox.test(df$Trimethyllysine, df$Population, p.adjust.method = "fdr")


#######################Load Data - 5 MAGS HEATMAP#################
setwd("/Users/user/Documents/02_Andropogon/03_SPIPS_Metagenome/00_DIABLO/5MAGS")

data <- read.table("5MAGs_Detection_HAYS.csv", sep=",", header=TRUE, row.names=1)

#pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE)
metadata <- read.table("5MAGs_Detection_HAYS_metadata.csv", sep=",", header=TRUE, row.names=1)
pc = metadata[,2:3]
pc2 = metadata[,1:3]
pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=pc, fontsize_col = 4, fontsize_row = 12)
pc3 = metadata[3]



svglite("HeatMap_Mags.svg")

plot <- pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=metadata, annotation_row=metadata, fontsize_col = 4, fontsize_row = 4)


my_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Create the heatmap
plot <- pheatmap(data,
                 cluster_rows = FALSE,
                 show_rownames = TRUE,
                 cluster_cols = TRUE,
                 annotation_col = pc2,
                 fontsize_col = 8,
                 fontsize_row = 10,
                 color = my_colors)



column_colors <- c("red", "blue", "green")  # Adjust based on your categories
names(column_colors) <- unique(pc2$YourCategory)


plot <-pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=pc2, fontsize_col = 8, fontsize_row = 10)