dev.off()
par(mfrow=c(1,1))
par(mar = rep(2, 4))
op <- par(oma=c(5,7,1,1))
par(op)
##Network analysis of the DIABLO results (same as circos plot but much more readable)
## To make this work you have to have the plot size as large as possible or you get a graphics or margins error
network(MyResult.diablo, color.node = c("grey70","tomato"), cutoff = 0.85)



#Metabolites_metaphaln 
#Metabolites_metaphaln 
#Metabolites_metaphaln 
#Metabolites_metaphaln 
######################################## START HERE ####################################################
#Metabolites_metaphaln 
#Anna 
#my working code
setwd("/Users/user/Documents/02_Andropogon/03_SPIPS_Metagenome/00_DIABLO")
rm(list = ls())
######################################## START HERE ####################################################
library(mixOmics)
Mix <- as.matrix(read.table("Mags_detections.csv", sep=",", header=TRUE, row.names=1))
Mix2 <- as.matrix(read.table("Metabolites_RNA_num.csv", sep=",", header=TRUE, row.names=1))
X <- list(MAGS = Mix,
          Metabolites = Mix2)
Y <- factor(c("Hyas","Hyas","Hyas","Hyas","Hyas","Hyas","Hyas","Hyas","Hyas","Hyas","Hyas","Hyas","IL","IL","IL","IL","IL","IL","IL","IL","IL","IL","IL","IL"))
MyResult.diablo <- block.splsda(X, Y)

plotIndiv(MyResult.diablo,          
          ind.names = FALSE, 
          legend=TRUE, cex=c(1,2),
          title = 'results')
chord <- circosPlot(MyResult.diablo, size.labels = 2 ,line = TRUE, size.legend = 0.9,size.variables = 0.7, cutoff=0.83, legend =TRUE)
plot(chord)

dev.off()
##Network analysis of the DIABLO results (same as circos plot but much more readable)
## To make this work you have to have the plot size as large as possible or you get a graphics or margins error
network.res <- network(MyResult.diablo, color.node = c("grey70","tomato"), save = 'jpeg', name.save = 'PLS_network_image', cutoff = 0.83) 
network.res <- network(MyResult.diablo, color.node = c("grey70","tomato"), name.save = 'PLS_network_image', cutoff = 0.83) 





