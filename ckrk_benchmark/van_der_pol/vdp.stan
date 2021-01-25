// Van Der Pol oscillator
// theta[0] is 1/mu^2 where mu is from standard definition of the system.
functions {
  real[] vdp(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    real dydt[2];
    dydt[1] = -y[2];
    dydt[2] = (-y[1] + y[2] * (1 - y[1] * y[1])) / theta[1];
    return dydt;
  }
}

data {
  int use_ckrk;
  int<lower=0> N_t;
  real t[N_t];
  real y1[N_t];
}

transformed data {
  real t0 = 0;
  real y0[2] = {2.0, 2.0/3.0};
  real x_r[0];
  int x_i[0];
}

parameters {
  real<lower=0> eps;
  real<lower=0> sigma;
}

transformed parameters {
  real y[N_t, 2];
  real y1hat[N_t];
  real theta[1] = {eps};
  if(use_ckrk == 1) {
    y = integrate_ode_adams(vdp, y0, t0, t, theta, x_r, x_i, 1e-6, 1e-8, 100000);
  } else {
    y = integrate_ode_rk45(vdp, y0, t0, t, theta, x_r, x_i, 1e-6, 1e-8, 100000);
  }
  y1hat = y[, 1];
}

model {
  eps ~ lognormal(log(0.1), 0.2);
  sigma ~ cauchy(0, 0.5);
  y1 ~ normal(y1hat, sigma);
}

// generated quantities {
//   real y_pred[N_t, 2];  
//   real y1_pred[N_t];  
//   real theta_gen[1] = {eps};
//   y_pred = integrate_ode_bdf(vdp, y0, t0, t, theta_gen, x_r, x_i, 1e-5, 1e-8, 5000);
//   y1_pred = normal_rng(y_pred[, 1], sigma);
// }
