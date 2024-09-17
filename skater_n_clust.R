

x = low.hab
min.clust=3
max.clust=30
index="Calinski_Harabasz"

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
    eval$Index[i-(min.clust-1)] <- unlist(intCriteria(x1, as.vector(as.integer(cr$Clusters)), crit=index))
  }
  
  # rbsste.plot <- ggplot(eval, aes(x=clusters,y=RBTSSE)) +
  #   geom_line() +
  #   theme_bw() +
  #   labs(x='Clusters', y="Goodness of classification (RBTSSE)")
  # index.plot <- ggplot(eval, aes(x=clusters,y=Index)) +
  #   geom_line() +
  #   theme_bw() +
  #   labs(x='Clusters', y=index)
  
  eval <- eval[which(!is.nan(eval$Index)),]
  best.want <- which.max(eval$Index == eval$Index[bestCriterion(eval$Index[!is.nan(eval$Index)], index)])
  
  # n.clust.want <- eval[bestCriterion(eval$Index, index),]
  n.clust.want <- eval[best.want, "clusters"]
  
  # return(list(eval=eval, want=n.clust.want, rbsste.plot=rbsste.plot, 
  #             index.plot=index.plot))
  return(n.clust.want)
}


### Indices
# "S_Dbw"
# "Calinski_Harabasz"
# "Davies_Bouldin"
# low.red.n.clust <- redcap.n.clust(low.hab, 8, 30, "Calinski_Harabasz")
# print(low.red.n.clust$rbsste.plot)
# print(low.red.n.clust$index.plot)
# low.red.n.clust$eval
# low.red.n.clust$want
skater.n.clust(low.hab, 3, 30, "Calinski_Harabasz")
