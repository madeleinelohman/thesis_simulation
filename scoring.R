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
set.seed(2025)

### Load in algorithm results
load("cluster_res.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Comparing internal scoring
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
crits <- c("Calinski_Harabasz", "Davies_Bouldin", "Silhouette")
sig.names <- c("low", "medium", "high")
for(i in 1:length(sig.names)){
  hab <- get(paste0(sig.names[i], ".hab"))
  int <- data.frame(Algorithm=c("Grid", "Voronoi", "SKATER", "REDCAP", "AZP"), 
                    Calinski_Harabasz=NA, Davies_Bouldin=NA,
                    Silhouette=NA)
  
  x <- matrix(values(hab), 50, 50, byrow=T)
  
  int[1,2:ncol(int)] <- intCriteria(x, as.vector(as.integer(get(paste0(sig.names[i], ".grid"))$clusters)), 
                                    crit=crits)
  int[2,2:ncol(int)] <- intCriteria(x, get(paste0(sig.names[i], ".v"))$cluster, 
                                    crit=crits)
  int[3,2:ncol(int)] <- intCriteria(x, get(paste0(sig.names[i], ".red"))$clusters, 
                                    crit=crits)
  int[4,2:ncol(int)] <- intCriteria(x, get(paste0(sig.names[i], ".azp"))$clusters, 
                                    crit=crits)
  
  int[which(is.na(as.matrix(int)), arr.ind=T)] <- 1e+15
  
  assign(paste0(sig.names[i], ".int"), int)
}

#~~~~~~~~~~~~~~~~~~~~
# Low heterogeneity habitat
#~~~~~~~~~~~~~~~~~~~~
low.int$Algorithm[bestCriterion(unlist(low.int$S_Dbw), "S_Dbw")]
low.int$Algorithm[bestCriterion(unlist(low.int$Calinski_Harabasz), "Calinski_Harabasz")]
low.int$Algorithm[bestCriterion(unlist(low.int$Davies_Bouldin), "Davies_Bouldin")]
low.int$Algorithm[bestCriterion(unlist(low.int$Davies_Bouldin), "Silhouette")]

#~~~~~~~~~~~~~~~~~~~~
# Medium heterogeneity habitat
#~~~~~~~~~~~~~~~~~~~~
medium.int$Algorithm[bestCriterion(unlist(medium.int$S_Dbw), "S_Dbw")]
medium.int$Algorithm[bestCriterion(unlist(medium.int$Calinski_Harabasz), "Calinski_Harabasz")]
medium.int$Algorithm[bestCriterion(unlist(medium.int$Davies_Bouldin), "Davies_Bouldin")]
medium.int$Algorithm[bestCriterion(unlist(medium.int$Davies_Bouldin), "Silhouette")]

#~~~~~~~~~~~~~~~~~~~~
# High heterogeneity habitat
#~~~~~~~~~~~~~~~~~~~~
high.int$Algorithm[bestCriterion(unlist(high.int$S_Dbw), "S_Dbw")]
high.int$Algorithm[bestCriterion(unlist(high.int$Calinski_Harabasz), "Calinski_Harabasz")]
high.int$Algorithm[bestCriterion(unlist(high.int$Davies_Bouldin), "Davies_Bouldin")]
high.int$Algorithm[bestCriterion(unlist(high.int$Davies_Bouldin), "Silhouette")]
