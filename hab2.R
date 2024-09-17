

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate habitat
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n = 25

sig1 = 3
phi1 = 6

sig2 = 1
phi2 = 3


r=rast(ncol=n*2,nrow=n*2,extent=c(-n,n,-n,n))
r[]=0
#s.r=crds(r)
rf1=grf(n, xyFromCell(r,1:ncell(r)), cov.pars=c(sig1, phi1))
rf2=grf(n, xyFromCell(r,1:ncell(r)), cov.pars=c(sig2, phi2))


hab1=r
hab1[]=rf1$data
hab2=r
hab2[]=rf2$data
hab <- c(hab1, hab2)

#plot(hab, main='Simulated habitat')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sample habitat
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ### Extract center points from raster
xy <- xyFromCell(hab, 1:ncell(hab))
xy <- as.data.frame(xy)

### Assign values from the raster to those associated center points
xy$values <- values(hab)
xy <- as.matrix(xy)

want <- matrix(values(hab), n*2, n*2)
total.mean <- mean(c(want))
  



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Redcap
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~
# Run method
#~~~~~~~~~~~~~~~~~~~~
v <- as.polygons(hab, aggregate=F)
s <- sf::st_as_sf(v)

n.clust = 10

queen_w <- queen_weights(s) 
rook_w <- rook_weights(s) 
data <- s[c('lyr.1')] 
cr <- redcap(n.clust, rook_w, data, "fullorder-completelinkage") 
cr


pred <- matrix(cr$Clusters, nrow=n*2)
s$pred <- cr$Clusters

plot(s["pred"])



pred.clust <- hab
values(pred.clust) <- cr$Clusters

gross <- cbind(values(hab), cr$Clusters)
colnames(gross)[1:2] <- c("org1", "org2")

ack <- gross %>%
  group_by(new) %>%
  summarize(vals = mean(org))

gross$m <- NA
g <- unique(ack$new)

for(i in 1:length(g)){
  gross[which(gross$new == g[i]), 3] <- ack[which(ack$new == g[i]), 2]
}

hab.new <- hab
values(hab.new) <- gross$m

par(mfrow=c(1,2))
plot(hab, main='Simulated habitat')
plot(hab.new, main="Predicted clusters")
par(mfrow=c(1,2))

png("plots/redcap_pred.png", 6, 6, "in", res=600)
plot(hab.cat4, main="Redcap predicted habitat")
dev.off()


#~~~~~~~~~~~~~~~~~~~~
# Validate method
#~~~~~~~~~~~~~~~~~~~~
ssw <- NA
for(i in 1:n.clust){
  ind.want <- which(pred == i, arr.ind=T)
  ssw[i] <- sum((want[ind.want] - mean(want[ind.want]))^2) 
}
ssw
sum(ssw)

ssb <- NA
for(i in 1:n.clust){
  ind.want <- which(pred == i, arr.ind=T)
  m.want <- mean(want[ind.want])
  ssb[i] <- (m.want - total.mean)^2 * nrow(ind.want)
}
ssb
sum(ssb)

ss.total <- sum(ssb + ssw)


redcap.prop <- sum(ssb) / ss.total
redcap.prop 



