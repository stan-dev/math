# Checking the neighbors...
rm(list=ls())
setwd("/Users/luzhang/Documents/Biostats/research/bitbucket/stan_math/nngp_dev/data")
n <- 200
m <- 8
coords <- matrix(read.table("coords.txt", quote="\"", comment.char="")$V1, 
                 nrow = n, ncol = 2)
nearind <- matrix(read.table("nearind.txt", quote="\"", comment.char="")$V1,
                  nrow = (n-1), ncol = m, byrow = T) 

par(mfrow = c(1, 1))
for(i in (n-10):n){
  plot(coords) # all points
  points(coords[1:i,,drop=FALSE], col="grey", pch=19) # points 1:i
  points(coords[i,,drop=FALSE], col="blue", pch=19) # the target point
  
  # neighborhood
  if (i < m) {dim = i} else {dim = m}
  for (j in 1:dim){
    points(coords[nearind[i-1,j],,drop=FALSE], col="orange", 
           pch=19)
    points(coords[nearind[i-1,j],,drop=FALSE], col="orange", 
           pch=as.character(j), cex = 2)
  }
  readline(prompt = "Pause. Press <Enter> to continue...")
}









