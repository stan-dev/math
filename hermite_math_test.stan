functions {
  tuple(matrix, vector) interp1_cubic_hermite_setup(vector xk, vector yk, vector dydxk);
  real interp1_cubic_hermite_eval(real x, tuple(matrix, vector) W);
}
data {
  int<lower=2> K;
  vector[K] xk; // Lookup-table x values; must be strictly increasing.
  vector[K] yk; // Lookup-table y values.
  vector[K] dydxk; // Lookup-table y values.
  int<lower=1> N_eval;
  vector[N_eval] x_eval;
}
transformed data {
  tuple(matrix[K - 1, 4], vector[K]) W = interp1_cubic_hermite_setup(xk, yk,dydxk);
}
generated quantities {
  // vector[K] pchip_slope = pchip_slopes(xk, yk);
  matrix[K - 1, 4] hermite_coef = W.1;
  vector[N_eval] y_eval;
  
  for (n in 1 : N_eval) {
    y_eval[n] = interp1_cubic_hermite_eval(x_eval[n], W);
  }
}
