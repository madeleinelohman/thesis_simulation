

azp.n.clust <- function(x, min.clust, max.clust, index){
  
  v <- as.polygons(x, aggregate=F) 
  s <- sf::st_as_sf(v)
  x1 <- matrix(values(x), nrow=50, ncol=50, byrow=T)
  
  ### Get input information
  queen_w <- queen_weights(s) # Queen weighted neighborhood
  data <- s[c('lyr.1')] # Data of interest
  
  eval <- data.frame(clusters=min.clust:max.clust, RBTSSE=NA, Index=NA)
  
  for(i in min.clust:max.clust){
    cr <- azp_tabu(i, queen_w, data) 
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
  
  best.want <- which.max(eval$Index == eval$Index[bestCriterion(eval$Index[!is.nan(eval$Index)], index)])
  
  # n.clust.want <- eval[bestCriterion(eval$Index, index),]
  n.clust.want <- eval[best.want, "clusters"]
  
  # return(list(eval=eval, want=n.clust.want, rbsste.plot=rbsste.plot, 
  #             index.plot=index.plot))
  return(n.clust.want)
}