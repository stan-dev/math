functions {

  vector linked_mass_flow(real t, vector y,
                          real kp, real ks, real e50, real k12, real k21) {

    int num_main_states = num_elements(y)/2;
    int num_states = 2 * num_main_states;

    vector[num_states] dydt;
    vector[num_main_states] mass_flow;

    for (i in 1:num_main_states) {
      int m = 2 * (i-1) + 1;  // main state
      int a = m + 1;  // auxilary state
      mass_flow[i] = (kp + ks  / (y[m] + e50)) * y[m];

      dydt[m] = -mass_flow[i] - k12 * y[m] + k21 * y[a];
      dydt[a] = +k12 * y[m] - k21 * y[a];

      if (i != 1) {
        dydt[m] += mass_flow[i-1];
      }
    }
    return dydt;

  }

  vector sample_vector_rng(real m, real s, int N) {
    vector[N] sample;
    for(i in 1:N)
      sample[i] = normal_rng(m, s);
    return sample;
  }

  vector[] simulate_mean(real t0, vector log_a0, real[] ts, int adjoint_integrator,
                         real rel_tol, real abs_tol, int max_num_steps,
                         int num_checkpoints,
                         int interpolation_polynomial,
                         int solver_f, int solver_b,
                         real log_kp, real log_ks, real log_e50, real log_k12, real log_k21) {
    int num_sim = num_elements(ts);
    int num_states = num_elements(log_a0);
    int system_size = num_states/2;
    vector[num_states] y[num_sim];
    real kp = exp(log_kp);
    real ks = exp(log_ks);
    real e50 = exp(log_e50);
    real k12 = exp(log_k12);
    real k21 = exp(log_k21);
    vector[num_states] a0 = exp(log_a0);

    if(adjoint_integrator) {
      y = ode_adjoint_tol_ctl(linked_mass_flow, a0, t0, ts, 
                              rel_tol, rep_vector(abs_tol/100.0, num_states), // forward
                              rel_tol, rep_vector(abs_tol/10.0, num_states), // backward
                              rel_tol, abs_tol, // quadrature
                              max_num_steps,
                              num_checkpoints, // number of steps between checkpoints
                              interpolation_polynomial,  // polynomials
                              solver_f,  // bdf forward
                              solver_b,  // bdf backward
                              kp, ks, e50, k12, k21);
      /* simplified interface
      y = ode_adjoint_tol(linked_mass_flow, a0, t0, ts, 
                              rel_tol, abs_tol,
                              max_num_steps,
                              kp, ks, e50, k12, k21);
      */
    } else {
      y = ode_bdf_tol(linked_mass_flow, a0, t0, ts, rel_tol, abs_tol, max_num_steps,
                      kp, ks, e50, k12, k21);
    }

    return y;
  }
  
}
data {
  real<lower=0> rel_tol;
  real<lower=0> abs_tol;
  int<lower=0,upper=1> adjoint_integrator;
  int<lower=1> max_num_steps;
  int<lower=1> num_checkpoints;
  int<lower=1,upper=2> interpolation_polynomial;
  int<lower=1,upper=2> solver_f;
  int<lower=1,upper=2> solver_b;
  int<lower=1> system_size;
  int<lower=1> num_obs;
  real<lower=0> sigma_sim;
  real<lower=0> sigma_y;
}
transformed data {
  int num_states = 2 * system_size;
  real param_scale = sigma_sim / sqrt(num_obs/4);

  real log_kp_ = normal_rng(0.0, sigma_sim);
  real log_ks_ = normal_rng(0.0, sigma_sim);
  real log_e50_ = normal_rng(0.0, sigma_sim);
  real log_k12_ = normal_rng(0.0, sigma_sim);
  real log_k21_ = normal_rng(0.0, sigma_sim);
  vector[system_size] log_sigma_y_ = sample_vector_rng(0.0, sigma_y, system_size);
  vector[num_states] log_a0_ = sample_vector_rng(0.0, sigma_sim, num_states);
  vector[num_states] y_[num_obs];
  real ts[num_obs];
  real t0 = 0.0;

  ts[1] = 1.0;
  for(i in 2:num_obs) {
    ts[i] = 1.5 * ts[i-1];
  }

  y_ = simulate_mean(t0, log_a0_, ts, 0, rel_tol, abs_tol, max_num_steps,
                     num_checkpoints,
                     interpolation_polynomial,
                     solver_f, solver_b,
                     log_kp_, log_ks_, log_e50_, log_k12_, log_k21_);
  
  for(i in 1:num_obs) {
    for(j in 1:system_size) {
      int m = 2*(j-1) + 1;
      y_[i,m] = lognormal_rng(log(y_[i,m]+1E-3), exp(log_sigma_y_[j]));
    }
  }

  if(adjoint_integrator) {
    print("Using adjoint integrator.");
  } else {
    print("Using bdf integrator.");
  }
  print("relative tolerance: ", rel_tol);
  print("absolute tolerance: ", abs_tol);
  print("maximum number of steps: ", max_num_steps);
  print("number of checkpoints: ", num_checkpoints);
  print("interpolation polynomial: ", interpolation_polynomial);
  print("solver forward: ", solver_f);
  print("solver backward: ", solver_b);
  print("number of time points: ", num_obs);
  print("system size: ", system_size);
  print("time points: ", ts);
  print("y_: ", y_);
  print("ODE states N: ", num_states);
  print("ODE parameters varying M: ", 5 + num_states);
}
parameters {
  real<multiplier=param_scale,offset=log_kp_> log_kp;
  real<multiplier=param_scale,offset=log_ks_> log_ks;
  real<multiplier=param_scale,offset=log_e50_> log_e50;
  real<multiplier=param_scale,offset=log_k12_> log_k12;
  real<multiplier=param_scale,offset=log_k21_> log_k21;
  vector<offset=log_sigma_y_>[system_size] log_sigma_y;
  vector<multiplier=param_scale,offset=log_a0_>[num_states] log_a0;
}
transformed parameters {
}
model {
  vector[num_states] mu[num_obs];

  profile("ode") {
    mu = simulate_mean(t0, log_a0, ts, adjoint_integrator, rel_tol, abs_tol, max_num_steps,
                       num_checkpoints,
                       interpolation_polynomial,
                       solver_f, solver_b,
                       log_kp, log_ks, log_e50, log_k12, log_k21);
  }
  
  target += normal_lpdf(log_kp| 0.0, sigma_sim);
  target += normal_lpdf(log_ks| 0.0, sigma_sim);
  target += normal_lpdf(log_e50| 0.0, sigma_sim);
  target += normal_lpdf(log_k12| 0.0, sigma_sim);
  target += normal_lpdf(log_k21| 0.0, sigma_sim);
  target += normal_lpdf(log_a0| 0.0, sigma_sim);
  target += normal_lpdf(log_sigma_y| 0.0, sigma_y);

  for(j in 1:system_size) {
    int m = 2*(j-1) + 1;
    target += lognormal_lpdf(to_vector(y_[:,m])| log(to_vector(mu[:,m])+1E-3), exp(log_sigma_y[j]));
  }
}
generated quantities {
  int rank_log_kp;
  int rank_log_ks;
  int rank_log_e50;
  int rank_log_k12;
  int rank_log_k21;
  int rank_log_sigma_y[system_size];
  int rank_log_a0[system_size];
  real bias_log_kp = log_kp - log_kp_;
  real bias_log_ks = log_ks - log_ks_;
  real bias_log_e50 = log_e50 - log_e50_;
  real bias_log_k12 = log_k12 - log_k12_;
  real bias_log_k21 = log_k21 - log_k21_;
  vector[system_size] bias_log_sigma_y = log_sigma_y - log_sigma_y_;
  vector[system_size] bias_log_a0 = log_a0 - log_a0_;

  rank_log_kp = (log_kp > log_kp_ ? 1 : 0);
  rank_log_ks = (log_ks > log_ks_ ? 1 : 0);
  rank_log_e50 = (log_e50 > log_e50_ ? 1 : 0);
  rank_log_k12 = (log_k12 > log_k12_ ? 1 : 0);
  rank_log_k21 = (log_k21 > log_k21_ ? 1 : 0);
  
  for(i in 1:system_size) {
    rank_log_sigma_y[i] = (log_sigma_y[i] > log_sigma_y_[i] ? 1 : 0);
    rank_log_a0[i] = (log_a0[i] > log_a0_[i] ? 1 : 0);
  }
  
}
