ward.alg <- function(x, level, n.clust, plotting=T){
  #~~~~~~~~~~~~~~~~~~~~
  # Run method
  #~~~~~~~~~~~~~~~~~~~~
  poly.x <- st_as_sf(as.polygons(x, round=F))
  list.nb <- poly2nb(poly.x) 
  A <- nb2mat(list.nb,style="B")
  diag(A) <- 1
  
  D1 <- as.dist(1-A)
  D0 <- dist(poly.x$lyr.1)
  
  
  range.alpha <- seq(0,1,0.05)
  cr <- choicealpha(D0, D1, range.alpha,
                    n.clust, graph=FALSE)
  # plot(cr)
  # plot(cr, norm = TRUE)
  
  best.alpha <- which.min(abs(cr$Qnorm[,1] - cr$Qnorm[,2]))
  
  tree <- hclustgeo(D0, D1, alpha=range.alpha[best.alpha])
  P5ter <- cutree(tree, n.clust)
  new <- x
  values(new) <- P5ter
  
  
  #~~~~~~~~~~~~~~~~~~~~
  # Plotting and plotting setup
  #~~~~~~~~~~~~~~~~~~~~
  ### Create a dataframe with old values of raster and new clusters
  pred.df <- data.frame(org=values(x), new=P5ter)
  names(pred.df)[1] <- "org"
  
  ### Summarize these clusters by their mean values
  summ.pred <- pred.df %>%
    group_by(new) %>%
    summarize(vals = mean(org))
  
  ### Associate this new mean with the clusters
  pred.df$new.mean <- NA
  clusts <- unique(summ.pred$new)
  for(i in 1:length(clusts)){
    pred.df[which(pred.df$new == clusts[i]), 3] <- summ.pred[which(summ.pred$new == clusts[i]), 2]
  }
  
  ### Create a new raster with these new means for each cluster
  new.hab <- x
  values(new.hab) <- pred.df$new.mean
  
  
  if(plotting == T){
    ### Plot!
    par(mfrow=c(1,2))
    plot(x, main=paste("Simulated habitat:", level))
    plot(new.hab, main=paste("Predicted clusters:", level))
    par(mfrow=c(1,1))
    
    p <- ggplot() +
      geom_spatraster(data=new.hab) +
      scale_fill_binned_divergingx(palette="Geyser") +
      theme_classic() + 
      labs(title=paste(level, "variation habitat"), subtitle="Ward's Clustering") +
      lims(x=c(-n, n), y=c(-n, n)) +
      guides(fill=guide_legend(title="Value"))
    
    png(paste0("plots/ward_plots/ward_",level,".png"), 6, 6, res=600, units="in")
    print(p)
    dev.off()
  }
  
  #~~~~~~~~~~~~~~~~~~~~
  # Validate method
  #~~~~~~~~~~~~~~~~~~~~
  ### The following code is similar to the grid algorithm validation
  rast.vals <- matrix(values(x), n*2, n*2)
  total.mean <- mean(rast.vals)
  pred <- matrix(values(new), nrow=n*2)
  
  ssw <- NA
  ssb <- NA
  for(i in 1:n.clust){
    ind.want <- which(pred == i, arr.ind=T)
    m.want <- mean(rast.vals[ind.want])
    
    ssw[i] <- sum((rast.vals[ind.want] - m.want)^2)
    ssb[i] <- (m.want - total.mean)^2 * nrow(ind.want)
  }
  ss.total <- sum(ssb + ssw)
  ssb.sst <- sum(ssb) / ss.total
  
  ### Return objects
  return(list(new.hab=new.hab, ssw=ssw, ssb=ssb, sst=ss.total, ssb.sst=ssb.sst,
              clusters=P5ter))
}





