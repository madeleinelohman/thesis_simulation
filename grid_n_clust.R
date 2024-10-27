

grid.n.clust <- function(x, index, n){
  ### Set up values for for loop
  rast.vals <- matrix(values(x), n*2, n*2) ### Make a matrix out of the raster values
  total.mean <- mean(values(x)) ### Find the mean of the raster
  
  ### What are potential window sizes?
  potential.windows <- seq(2, (n*2)-1)
  modulos <- (n*2) %% potential.windows
  grp.lvl <- potential.windows[which(modulos == 0)] # Window to group cells by (new grid dimension sizes)
  
  eval <- data.frame(clusters=grp.lvl, Index=NA)
  
  for(k in 1:length(grp.lvl)){
    ### Set up matrices for for loop
    grd <- matrix(NA, n*2/grp.lvl[k], n*2/grp.lvl[k]) # New grid
    clusters <- matrix(0, n*2, n*2) # Cluster assignment
    
    lats <- longs <- seq(grp.lvl[k], n*2, by=grp.lvl[k]) # Lat/long indices for each grid cell
    
    for(i in 1:length(lats)){
      for(j in 1:length(longs)){
        ### Identify the cell ranges to aggregate
        ind.lat <- lats[i]:(lats[i]-(grp.lvl[k]-1)) 
        ind.long <- longs[j]:(longs[j]-(grp.lvl[k]-1))
        
        ### Create grid
        mean.local.vals <- mean(rast.vals[ind.lat, ind.long]) # Find the mean from the original raster
        grd[i,j] <- mean.local.vals # Place the mean in the new grid
        clusters[ind.lat, ind.long] <- max(clusters) + 1 # Assign a cluster number
      }
    }
    new.clusts <- order.clusts(c(clusters))
    
    eval$Index[k] <- unlist(intCriteria(rast.vals, new.clusts$new, crit=index))
  }
  
  best.want <- which.max(eval$Index == eval$Index[bestCriterion(eval$Index[!is.nan(eval$Index)], index)])
  n.clust.want <- eval[best.want, "clusters"]
  
  return(n.clust.want)
}