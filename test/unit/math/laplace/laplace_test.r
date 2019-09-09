# Unit test for Laplace approximation

theta <- c(-0.191826, 0.191826)

W_root <- sqrt(diag(exp(theta) / (1 + exp(theta))^2))
K <- matrix(c(1, 0.575785, 0.575785, 1), nrow = 2)
B <- diag(c(1, 1)) + W_root %*% K %*% W_root
L <- t(chol(B))
Z <- W_root %*% solve(t(L), solve(L, W_root))

C <- solve(L, W_root %*% K)
third_tensor <- exp(theta) * (exp(theta) - 1) / (exp(theta) + 1)^3
s2 <- - 0.5 * (diag(K) - diag(t(C) %*% C)) * third_tensor

# derivatives of covariance matrix
C_0 <- matrix(c(2, 1.15157, 1.15157, 2), ncol = 2)
a <- c(-0.45219, 0.45219)
s1 <- 0.5 * t(a) %*% C_0 %*% a - 0.5 * sum(diag(Z %*% C_0))

gradient <- c(NA, NA)
gradient[1] <- 0 - 1 / (1 + exp(-theta[1]))
gradient[2] <- 1 - 1 / (1 + exp(-theta[2]))

b <- C_0 %*% gradient
s3 <- b - K %*% Z %*% b

adj <- s1 + t(s2) %*% s3

## Evaluate a
W = W_root^2
b = W %*% theta + gradient
a = b - W_root %*% solve(t(L), solve(L, W_root %*% K %*% b))

x <- c(1, 1)

#####################################################################
# log density and derivatives for poisson log model with exposure 
# term.

theta <- c(1, 1)
n_samples <- c(1, 1)
sums <- c(1, 0)
exposure <- c(0.5, 2)
eta <- log(exposure)

# log density
sum(-log(factorial(sums)) + (theta + eta) * sums 
      - n_samples * exp(theta + eta))
# -7.488852

# gradient
sums - n_samples * exp(theta + eta)
# -0.3591409 -5.4365637

# diagonal hessian
- n_samples * exp(theta + eta)
# -1.359141 -5.436564

# diagonal third-order derivative tenssor
- n_samples * exp(theta + eta)
# -1.359141 -5.436564



