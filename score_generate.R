#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries and set up
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ClustGeo)
library(clusterCrit)
library(colorspace)
library(cowplot)
library(dbscan)
library(geoR)
library(MASS)
library(mvtnorm)
library(raster)
library(rgeoda)
library(sf)
library(spdep)
library(tidyverse)
library(tidyterra)
library(terra)

setwd("/Users/madelienelohman/Documents/GitHub/thesis_simulation")
source("order_clusts.R")

set.seed(2024)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define variables of interest
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

indices <- c("Calinski_Harabasz", "Davies_Bouldin", "Dunn", "SD_Scat", "Silhouette")
sig.names <- c("low", "medium", "high")

for(k in 1:length(indices)){
  index=indices[k]
  
  ### Load in algorithm results and simulated rasters
  load(paste0("cluster_results/cluster_res_",index,".RData"))
  low.hab <- readRDS("habitats/low_hab.rds")
  medium.hab <- readRDS("habitats/med_hab.rds")
  high.hab <- readRDS("habitats/high_hab.rds")
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Comparing internal scoring
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(i in 1:length(sig.names)){
    hab <- get(paste0(sig.names[i], ".hab"))
    int <- data.frame(Algorithm=c("Grid", "Voronoi", "AZP", "SKATER", "REDCAP"), 
                      Clusters = NA,
                      Score=NA, RS=NA, var.level=sig.names[i], Index=index)
    
    x <- matrix(values(hab), n*2, n*2, byrow=T)
    
    grid.clusts <- order.clusts(c(get(paste0(sig.names[i], ".grid"))$cluster))
    int$Score[1] <- unlist(intCriteria(x, grid.clusts$new,
                                       crit=index))
    v.clusts <- order.clusts(c(get(paste0(sig.names[i], ".v"))$cluster))
    int$Score[2] <- unlist(intCriteria(x, v.clusts$new,
                                       crit=index))
    azp.clusts <- order.clusts(get(paste0(sig.names[i], ".azp"))$clusters)
    int$Score[3] <- unlist(intCriteria(x, azp.clusts$new,
                                       crit=index))
    skater.clusts <- order.clusts(get(paste0(sig.names[i], ".skater"))$clusters)
    int$Score[4] <- unlist(intCriteria(x, skater.clusts$new, 
                                       crit=index))
    red.clusts <- order.clusts(get(paste0(sig.names[i], ".red"))$clusters)
    int$Score[5] <- unlist(intCriteria(x, red.clusts$new, 
                                       crit=index))
    
    int$RS[1] <- get(paste0(sig.names[i], ".grid"))$ssb.sst
    int$RS[2] <- get(paste0(sig.names[i], ".v"))$ssb.sst
    int$RS[3] <- get(paste0(sig.names[i], ".azp"))$ssb.sst
    int$RS[4] <- get(paste0(sig.names[i], ".skater"))$ssb.sst
    int$RS[5] <- get(paste0(sig.names[i], ".red"))$ssb.sst
    
    int$Clusters[1] <- length(unique(as.vector(as.integer(get(paste0(sig.names[i], ".grid"))$clusters))))
    int$Clusters[2] <- length(unique(get(paste0(sig.names[i], ".v"))$clusters))
    int$Clusters[5] <- length(unique(get(paste0(sig.names[i], ".azp"))$clusters))
    int$Clusters[3] <- length(unique(get(paste0(sig.names[i], ".skater"))$clusters))
    int$Clusters[4] <- length(unique(get(paste0(sig.names[i], ".red"))$clusters))
    
    
    #int[which(is.na(as.matrix(int)), arr.ind=T)] <- 1e+15
    
    int[,3:4] <- round(int[,3:4], 3)
    
    assign(paste0(sig.names[i], ".int.", index), int)
  }

}

rm(list=setdiff(ls(), ls()[sapply(ls(), str_detect, pattern=".int")]))

save.image("scoring_res.RData")


