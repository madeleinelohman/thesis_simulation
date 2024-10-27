skater.n.clust <- function(x, min.clust, max.clust, index){
  #~~~~~~~~~~~~~~~~~~~~
  # Run method
  #~~~~~~~~~~~~~~~~~~~~
  ### Transform raster to sf class
  # Need to go through sp class to get to sf class from SpatRaster 
  v <- as.polygons(x, aggregate=F) 
  s <- sf::st_as_sf(v)
  x1 <- matrix(values(x), nrow=50, ncol=50, byrow=T)
  
  ### Get input information
  queen_w <- queen_weights(s) # Queen weighted neighborhood
  data <- s[c('lyr.1')] # Data of interest
  
  eval <- data.frame(clusters=min.clust:max.clust, RBTSSE=NA, Index=NA)
  
  for(i in min.clust:max.clust){
    cr <- rgeoda::skater(i, queen_w, data, cpu_threads = 1) 
    eval$RBTSSE[i-(min.clust-1)] <- cr$`The ratio of between to total sum of squares`
    
    new.clusts <- order.clusts(cr$Clusters)
    eval$Index[i-(min.clust-1)] <- unlist(intCriteria(x1, new.clusts$new, crit=index))
  }
  
  
  eval <- eval[which(!is.nan(eval$Index) & eval$Index < 1*10^6),]
  # eval <- eval[which(!is.nan(eval$Index)),]
  best.want <- which.max(eval$Index == eval$Index[bestCriterion(eval$Index, index)])
  n.clust.want <- eval[best.want, "clusters"]
  
  return(n.clust.want)
}


