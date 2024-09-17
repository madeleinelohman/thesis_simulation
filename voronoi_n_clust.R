

v.n.clust <- function(x, min.clust, max.clust, index){
  
  v <- as.polygons(x, aggregate=F) 
  s <- sf::st_as_sf(v)
  x1 <- matrix(values(x), nrow=50, ncol=50, byrow=T)
  
  ### Get input information
  queen_w <- queen_weights(s) # Queen weighted neighborhood
  data <- s[c('lyr.1')] # Data of interest
  
  eval <- data.frame(clusters=min.clust:max.clust, Index=NA)
  
  for(i in min.clust:max.clust){
    samps <- terra::spatSample(x, size=i, method="random", as.df=T, xy=T)
    samps <- st_as_sf(samps, coords=c('x', 'y'))
    
    
    g <- st_combine(st_geometry(samps)) 
    v <- st_voronoi(g)
    v <- st_collection_extract(v)
    v <- v[unlist(st_intersects(samps, v))]
    v <- st_as_sf(v)
    
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
