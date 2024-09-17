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


source("grid_alg.R")
source("voronoi_alg.R")
source("SKATER_alg.R")
source("REDCAP_alg.R")
source("azp_alg.R")
source("ward_alg.R")

### Read in simulated rasters
low.hab <- readRDS("low_hab.rds")
medium.hab <- readRDS("med_hab.rds")
high.hab <- readRDS("high_hab.rds")


n = dim(low.hab)[1]/2


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grid
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
low.grid <- grid.alg(low.hab, "Low", n)
low.grid$ssb.sst

medium.grid <- grid.alg(medium.hab, "Medium", n)
medium.grid$ssb.sst

high.grid <- grid.alg(high.hab, "High", n)
high.grid$ssb.sst


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Voronoi tesselation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
low.v <- v.alg(low.hab, "Low", n, n)
print(low.v$p)
low.v$ssb.sst

medium.v <- v.alg(medium.hab, "Medium", n, n)
print(med.v$p)
med.v$ssb.sst

high.v <- v.alg(high.hab, "High", n, n)
print(high.v$p)
high.v$ssb.sst


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AZP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
low.azp <- azp.alg(low.hab, "Low", 3, 30, "Calinski_Harabasz")
low.azp$ssb.sst

medium.azp <- azp.alg(medium.hab, "Medium", 3, 30, "Calinski_Harabasz")
medium.azp$ssb.sst

high.azp <- azp.alg(high.hab, "High", 3, 30, "Calinski_Harabasz")
high.azp$ssb.sst


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SKATER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
low.skater <- skater.alg(low.hab, "Low", 3, 30, "Calinski_Harabasz", n)
low.skater$ssb.sst

medium.skater <- skater.alg(medium.hab, "Medium", 3, 30, "Calinski_Harabasz", n)
medium.skater$ssb.sst

high.skater <- skater.alg(high.hab, "High", 3, 30, "Calinski_Harabasz", n)
high.skater$ssb.sst


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REDCAP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
low.red <- redcap.alg(low.hab, "Low", 3, 30, "Calinski_Harabasz", n)
low.red$ssb.sst
low.red$n.clust

medium.red <- redcap.alg(medium.hab, "Medium", 3, 30, "Calinski_Harabasz", n)
medium.red$ssb.sst

high.red <- redcap.alg(high.hab, "High", 3, 30, "Calinski_Harabasz", n)
high.red$ssb.sst


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save everything
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save.image("cluster_res.RData")





