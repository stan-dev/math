functions {
  vector dz_dt(real t, vector z, real alpha, real beta, real gamma, real delta) {
    real u = z[1];
    real v = z[2];

    real du_dt = (alpha - beta * v) * u;
    real dv_dt = (-gamma + delta * u) * v;

    vector[2] dydt;
    dydt[1] = du_dt;
    dydt[2] = dv_dt;
    return dydt;
  }
}
data {
  real rtol;
  real atol;
  int<lower = 0> N;          // number of measurement times
  real ts[N];                // measurement times > 0
  real y_init[2];            // initial measured populations
  real<lower = 0> y[N, 2];   // measured populations
}

parameters {
  real<lower = 0> theta[4];   // { alpha, beta, gamma, delta }
  vector<lower = 0>[2] z_init;  // initial population
  real<lower = 0> sigma[2];   // measurement errors
}
transformed parameters {
  vector<lower = 0>[2] z[N] = ode_bdf_tol(dz_dt, z_init, 0.0, ts, rtol, atol, 500, theta[1], theta[2], theta[3], theta[4]);
}
model {
  theta[{1, 3}] ~ normal(1, 0.5);
  theta[{2, 4}] ~ normal(0.05, 0.05);
  sigma ~ lognormal(-1, 1);
  z_init ~ lognormal(log(10), 1);
  // for (k in 1:2) {
  //   y_init[k] ~ lognormal(log(z_init[k]), sigma[k]);
  //   y[ , k] ~ lognormal(log(z[, k]), sigma[k]);
  // }
}
/* generated quantities { */
/*   real y_init_rep[2]; */
/*   real y_rep[N, 2]; */
/*   for (k in 1:2) { */
/*     y_init_rep[k] = lognormal_rng(log(z_init[k]), sigma[k]); */
/*     for (n in 1:N) */
/*       y_rep[n, k] = lognormal_rng(log(z[n, k]), sigma[k]); */
/*   } */
/* } */
