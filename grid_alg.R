grid.alg <- function(x, level, n){
  #~~~~~~~~~~~~~~~~~~~~
  # Run method
  #~~~~~~~~~~~~~~~~~~~~
  rast.vals <- matrix(values(x), n*2, n*2) ### Make a matrix out of the raster values
  total.mean <- mean(values(x)) ### Find the mean of the raster
  
  ### Set up matrices for the for loop
  grp.lvl = 5 # Window to group cells by (new grid dimension sizes)
  grd <- matrix(NA, n*2/grp.lvl, n*2/grp.lvl) # New grid
  ssw <- matrix(NA, n*2/grp.lvl, n*2/grp.lvl) # Within-cluster sum of squares
  ssb <- matrix(NA, n*2/grp.lvl, n*2/grp.lvl) # Between-cluster sum of squares
  clusters <- matrix(0, n*2, n*2) # Cluster assignment
  
  lats <- longs <- seq(grp.lvl, n*2, by=grp.lvl) # Lat/long indices for each grid cell
  
  for(i in 1:length(lats)){
    for(j in 1:length(longs)){
      ### Identify the cell ranges to aggregate
      ind.lat <- lats[i]:(lats[i]-(grp.lvl-1)) 
      ind.long <- longs[j]:(longs[j]-(grp.lvl-1))
      
      ### Create grid
      mean.local.vals <- mean(rast.vals[ind.lat, ind.long]) # Find the mean from the original raster
      grd[i,j] <- mean.local.vals # Place the mean in the new grid
      clusters[ind.lat, ind.long] <- max(clusters) + 1 # Assign a cluster number
      
      ### Validate method
      ssw[i,j] <- sum((rast.vals[ind.lat, ind.long] - mean.local.vals)^2) # Calculate individual cell SSW
      ssb[i,j] <- (mean.local.vals - total.mean)^2 * grp.lvl^2 # Calculate individual cell SSB
    }
  }
  ss.total <- sum(ssb + ssw) # Total within-cluster sum of squares
  ssb.sst <- sum(ssb) / ss.total # The ratio of between to total sum of squares
  
  #~~~~~~~~~~~~~~~~~~~~
  # Rasterize results
  #~~~~~~~~~~~~~~~~~~~~
  new.r=rast(ncol=n*2/5,nrow=n*2/5,extent=c(-n,n,-n,n))
  new.r[]=0
  new.s.r=crds(new.r)
  new.hab=new.r
  new.hab[]=c(grd)
  
  #~~~~~~~~~~~~~~~~~~~~
  # Plot
  #~~~~~~~~~~~~~~~~~~~~
  # png("plots/grid_pred.png", 6, 6, "in", res=600)
  # plot(new.hab, main='Grided habitat')
  # dev.off()
  
  par(mfrow=c(1,2))
  plot(x, main=paste("Simulated habitat:", level))
  plot(new.hab, main=paste("Grided habitat:", level))
  par(mfrow=c(1,1))
  
  p <- ggplot() +
    geom_spatraster(data=new.hab) +
    scale_fill_binned_divergingx(palette="Geyser") +
    theme_classic() + 
    labs(title=paste(level, "variation habitat"), subtitle="Grid") +
    lims(x=c(-n, n), y=c(-n, n)) +
    guides(fill=guide_legend(title="Value"))
  
  png(paste0("plots/grid_plots/grid_",level,".png"), 6, 6, res=600, units="in")
  print(p)
  dev.off()
  
  ### Return objects
  return(list(new.hab=new.hab, ssw=ssw, ssb=ssb, sst=ss.total, ssb.sst=ssb.sst,
              clusters=clusters, p=p))
}
