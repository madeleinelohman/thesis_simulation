
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DBSCAN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(colorspace)
library(dbscan)



xy <- st_coordinates(rec.bands)
mp <- 4

kNNdistplot(xy, minPts=mp) # or 3 for minPts

e <- 35000

abline(h=e,col="red",lty=2)



res <- dbscan(xy, eps=e, minPts=mp) # eps potentially as 2 or 3
res



rec.bands$cluster <- as.factor(res$cluster)
ggplot(rec.bands, aes(color=cluster)) +
  scale_color_discrete_divergingx("Geyser") +
  geom_sf()

