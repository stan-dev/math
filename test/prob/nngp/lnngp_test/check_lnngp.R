# check the nngp loglikehood #
# set the working directory to test folder "nngp" #
setwd("./nngp")
rm(list=ls())
load("./data/nngp.RData")

nngp_log <- function(y, X, beta, sigmasq, tausq, phi, neardist, 
                     neardistM, nearind, M){
  kappa_p_1 = tausq/sigmasq +1;
  
  V <- rep(kappa_p_1, n)
  Uw <- y - X%*%beta
  temp_w <- Uw
  
  for (i in 1:(n-1)){
    d <- ifelse(i < M, i, M)
    temp_distM <- exp(-phi * matrix(
      c(neardistM[i, 1:(d * d)]), nrow = d, ncol = d)) + 
      (tausq/sigmasq)*diag(d)
    
    U <- chol(temp_distM)
    u <- exp(-phi * neardist[i, 1: d ]) 
    
    v <- forwardsolve(t(U), u)
    
    V[i+1] = kappa_p_1 - sum(v^2)
    
    v2 <- backsolve(U, v)
    
    for (j in 1: d){
      Uw[i+1] = Uw[i+1] - v2[j] * temp_w[nearind[i, j]]
    }
  }
  out = -0.5 * (sum(log(V)) + 1 / sigmasq * sum(Uw^2 / V) + n * log(sigmasq))
  
  return(out)
}

beta = c(1.0, 5.0)
sigma = 1
tau = sqrt(0.1)
phi = 12
nngp_log(y, X, beta, sigmasq, tausq, phi, neardist, neardistM, nearind, m)

#-64.62981












