# set working directory into the data folder #
setwd("./data")
rm(list=ls())

#----------------------------    Get Data    ----------------------------------#
rmvn <- function(n, mu = 0, V = matrix(1)){
  
  # function for generating simulation of multivariate normal distribution 
  #  mu: mean vector;  V: Covariance matrix
  
  p <- length(mu)
  if(any(is.na(match(dim(V), p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

set.seed(123)
n <- 200                                # number of observations 
coords <- cbind(runif(n), runif(n))     # get random locations in a unit square
coords <- coords[order(coords[, 1]), ]  # sort by the first coordinate
X <- as.matrix(cbind(1, rnorm(n)))      # design matrix, explanatory variable 
                                        # are from N(0,1)

B <- as.matrix(c(1, 5))                 # intercept and slope
sigma.sq <- 1                           # partial sill
tau.sq <- 0.1                           # nuggets
phi <- 12                               # 1/range

D <- as.matrix(dist(coords))
R <- exp(- phi * D)
w <- rmvn(1, rep(0, n), sigma.sq * R)   # spatial residual
y <- rnorm(n, X%*%B + w, sqrt(tau.sq))  # response

#------------------------    Get Extra Input    -------------------------------#
# get the distance matrix and index matriox of the nearest neighbor 
# for each observed location
# Only need sorted coordinate (coord) and number of neighbors (m)

library(fields)
m = 8                                   # number of neighbors

i_index <- function(i, s, m) {
  
  # get the neighbor index for each location
  #    s: location matrix, each row records coordinates of one point
  #    m: number of neighbor
  
  if(m >= (i - 1)) {im <- 1:(i - 1)}
  else 	{
    dist <- rdist(s[c(1, i), ], s[c(1:(i - 1)), ])[-1, ]
    im <- sort(order(dist)[1:m])
  }
  return(im)
}

#### distance matrix for location i and its neighbors ####
i_dist <- function(i, neighbor_index, s)
  dist(s[c(i, neighbor_index[[i - 1]]), ])

get_index_dist <- function(s, m) {
  
  n = nrow(s)
  m = min(m, n - 1)
  
  # get index of neighborhood
  neighbor_index <- sapply(2:n, i_index, s, m)
  
  # get distance vector of each i and its neighbors
  neighbor_dist <- sapply(2:n, i_dist, neighbor_index, s)
  
  return(list(i = neighbor_index, d = neighbor_dist))
}

ind_distM <- get_index_dist(coords, m)

#-------- get the input for NNGP ----------#
get_neardistM <- function (ind, ind_distM_d) {
  if (ind < m ){l = ind } else {l = m} 
  M_i <- rep(0, m * m);
  M_i[1: (l * l)]<-  c(as.matrix(ind_distM_d[[ind]])[-1, -1])
  return(M_i)
}

get_neardist <- function (ind, ind_distM_d) {
  if (ind < m ){l = ind } else {l = m} 
  D_i <- rep(0, m)
  D_i[1:l] <- c(ind_distM_d[[ind]])[1:l]
  return(D_i)
}

get_nearind <- function (ind, ind_distM_i) {
  if (ind < m ){l = ind } else {l = m} 
  D_i <- rep(0, m)
  D_i[1:l] <- c(ind_distM_i[[ind]])[1:l]
  return(D_i)
}

# distance matrix among neighbors 
neardistM <- sapply(1:(n-1), get_neardistM, ind_distM$d) 
# distance of i th location and its' neighbor
neardist <- sapply(1:(n-1), get_neardist, ind_distM$d)  
# neighbor index
nearind <- sapply(1: (n-1), get_nearind, ind_distM$i)   


#------------------------     Save the data     -------------------------------#

# save in TXT
write.table(c(coords), "coords.txt", row.names=F, col.names=F, sep="\t")
write.table(c(X), "X.txt", row.names=F, col.names=F, sep="\t")
write.table(c(y), "y.txt", row.names=F, col.names=F, sep="\t")
write.table(c(w), "w.txt", row.names=F, col.names=F, sep="\t")
write.table(c(D), "dist.txt", row.names=F, col.names=F, sep="\t")

write.table(c(neardistM), "neardistM.txt", row.names=F, col.names=F, sep="\t")
write.table(c(neardist), "neardist.txt", row.names=F, col.names=F, sep="\t")
write.table(c(nearind), "nearind.txt", row.names=F, col.names=F, sep="\t")

# save R data
neardistM <- t(neardistM)
neardist <- t(neardist)
nearind <- t(nearind)

save(list = ls(all.names = TRUE), file = "nngp.RData", envir = .GlobalEnv)




