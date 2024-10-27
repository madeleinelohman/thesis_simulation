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
set.seed(2024)


source("algs/grid_alg.R")
source("algs/voronoi_alg.R")
source("algs/azp_alg.R")
source("algs/SKATER_alg.R")
source("algs/REDCAP_alg.R")
source("order_clusts.R")

### Read in simulated rasters
low.hab <- readRDS("habitats/low_hab.rds")
medium.hab <- readRDS("habitats/med_hab.rds")
high.hab <- readRDS("habitats/high_hab.rds")


n = dim(low.hab)[1]/2
min.clust = 3
max.clust = 50

# "Calinski_Harabasz"
# "Davies_Bouldin"
# "Dunn"
# "SD_Scat"
# "Silhouette"

start.time <- Sys.time()
indices <- c("Calinski_Harabasz", "Davies_Bouldin", "Dunn", "SD_Scat", "Silhouette")

for(k in 1:length(indices)){
  
  index=indices[k]


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Grid
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  low.grid <- grid.alg(low.hab, "Low", n, index)
  low.grid$ssb.sst
  
  medium.grid <- grid.alg(medium.hab, "Medium", n, index)
  medium.grid$ssb.sst
  
  high.grid <- grid.alg(high.hab, "High", n, index)
  high.grid$ssb.sst
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Voronoi tesselation
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n.iter=25
  
  low.v <- v.alg(low.hab, "Low", min.clust, max.clust, n, index, n.iter)
  print(low.v$p)
  low.v$ssb.sst
  
  medium.v <- v.alg(medium.hab, "Medium", min.clust, max.clust, n, index, n.iter)
  print(medium.v$p)
  medium.v$ssb.sst
  
  high.v <- v.alg(high.hab, "High", min.clust, max.clust, n, index, n.iter)
  print(high.v$p)
  high.v$ssb.sst
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # AZP
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  low.azp <- azp.alg(low.hab, "Low", min.clust, max.clust, index)
  low.azp$ssb.sst
  
  medium.azp <- azp.alg(medium.hab, "Medium", min.clust, max.clust, index)
  medium.azp$ssb.sst
  
  high.azp <- azp.alg(high.hab, "High", min.clust, max.clust, index)
  high.azp$ssb.sst
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # SKATER
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  low.skater <- skater.alg(low.hab, "Low", min.clust, max.clust, index, n)
  low.skater$ssb.sst
  
  medium.skater <- skater.alg(medium.hab, "Medium", min.clust, max.clust, index, n)
  medium.skater$ssb.sst
  
  high.skater <- skater.alg(high.hab, "High", min.clust, max.clust, index, n)
  high.skater$ssb.sst
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # REDCAP
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  low.red <- redcap.alg(low.hab, "Low", min.clust, max.clust, index, n)
  low.red$ssb.sst
  
  medium.red <- redcap.alg(medium.hab, "Medium", min.clust, max.clust, index, n)
  medium.red$ssb.sst
  
  high.red <- redcap.alg(high.hab, "High", min.clust, max.clust, index, n)
  high.red$ssb.sst
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Save everything
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  save.image(paste0("cluster_results/cluster_res_",index,".RData"))
  
  
  end.time <- Sys.time()
  end.time - start.time

}