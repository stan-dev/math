functions {

  vector linked_mass_flow(real t, vector y,
                          real[] kp, real[] ks, real[] e50, real[] k12, real[] k21) {

    int num_main_states = num_elements(kp);
    int num_states = 2 * num_main_states;

    vector[num_states] dydt;
    vector[num_main_states] mass_flow;

    for (i in 1:num_main_states) {
      int m = 2 * (i-1) + 1;  // main state
      int a = m + 1;  // auxilary state
      mass_flow[i] = (kp[i] + ks[i]  / (y[m] + e50[i])) * y[m];

      dydt[m] = -mass_flow[i] - k12[i] * y[m] + k21[i] * y[a];
      dydt[a] = +k12[i] * y[m] - k21[i] * y[a];

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
                         real rel_tol_f, real abs_tol_f,
                         real rel_tol_b, real abs_tol_b,
                         real rel_tol_q, real abs_tol_q,
                         int max_num_steps,
                         int num_checkpoints,
                         int interpolation_polynomial,
                         int solver_f, int solver_b,
                         vector log_kp, vector log_ks, vector log_e50, vector log_k12, vector log_k21) {
    int num_sim = num_elements(ts);
    int num_states = num_elements(log_a0);
    int system_size = num_elements(log_kp);
    vector[num_states] y[num_sim];
    real kp[system_size];
    real ks[system_size];
    real e50[system_size];
    real k12[system_size];
    real k21[system_size];
    vector[num_states] a0;

    for(i in 1:system_size) {
      kp[i] = exp(log_kp[i]);
      ks[i] = exp(log_ks[i]);
      e50[i] = exp(log_e50[i]);
      k12[i] = exp(log_k12[i]);
      k21[i] = exp(log_k21[i]);
    }

    for(i in 1:num_states) {
        a0[i] = exp(log_a0[i]);
    }

    if(adjoint_integrator) {
      y = ode_adjoint_tol_ctl(linked_mass_flow, a0, t0, ts, 
                              rel_tol_f, rep_vector(abs_tol_f, num_states), // forward
                              rel_tol_b, rep_vector(abs_tol_b, num_states), // backward
                              rel_tol_q, abs_tol_q, // quadrature
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
      y = ode_bdf_tol(linked_mass_flow, a0, t0, ts,
                      rel_tol_f, abs_tol_f,
                      max_num_steps,
                      kp, ks, e50, k12, k21);
    }

    return y;
  }
  
}
data {
  real<lower=0> rel_tol_f;
  real<lower=0> abs_tol_f;
  real<lower=0> rel_tol_b;
  real<lower=0> abs_tol_b;
  real<lower=0> rel_tol_q;
  real<lower=0> abs_tol_q;
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

  vector[system_size] log_kp_ = sample_vector_rng(0.0, sigma_sim, system_size);
  vector[system_size] log_ks_ = sample_vector_rng(0.0, sigma_sim, system_size);
  vector[system_size] log_e50_ = sample_vector_rng(0.0, sigma_sim, system_size);
  vector[system_size] log_k12_ = sample_vector_rng(0.0, sigma_sim, system_size);
  vector[system_size] log_k21_ = sample_vector_rng(0.0, sigma_sim, system_size);
  vector[system_size] log_sigma_y_ = sample_vector_rng(0.0, sigma_y, system_size);
  vector[num_states] log_a0_ = sample_vector_rng(0.0, sigma_sim, num_states);
  vector[num_states] y_[num_obs];
  real ts[num_obs];
  real t0 = 0.0;

  ts[1] = 1.0;
  for(i in 2:num_obs) {
    ts[i] = 1.5 * ts[i-1];
  }

  y_ = simulate_mean(t0, log_a0_, ts, 0,
                     rel_tol_f, abs_tol_f,
                     rel_tol_b, abs_tol_b,
                     rel_tol_q, abs_tol_q,
                     max_num_steps,
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
  print("relative tolerance forward: ", rel_tol_f);
  print("absolute tolerance forward: ", abs_tol_f);
  print("relative tolerance backward: ", rel_tol_b);
  print("absolute tolerance backward: ", abs_tol_b);
  print("relative tolerance quadrature: ", rel_tol_q);
  print("absolute tolerance quadrature: ", abs_tol_q);
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
  print("ODE parameters varying M: ", 5*system_size + num_states);
}
parameters {
  vector<multiplier=param_scale,offset=log_kp_>[system_size] log_kp;
  vector<multiplier=param_scale,offset=log_ks_>[system_size] log_ks;
  vector<multiplier=param_scale,offset=log_e50_>[system_size] log_e50;
  vector<multiplier=param_scale,offset=log_k12_>[system_size] log_k12;
  vector<multiplier=param_scale,offset=log_k21_>[system_size] log_k21;
  vector<offset=log_sigma_y_>[system_size] log_sigma_y;
  vector<multiplier=param_scale,offset=log_a0_>[num_states] log_a0;
}
transformed parameters {
}
model {
  vector[num_states] mu[num_obs];

  profile("ode") {
    mu = simulate_mean(t0, log_a0, ts, adjoint_integrator,
                       rel_tol_f, abs_tol_f,
                       rel_tol_b, abs_tol_b,
                       rel_tol_q, abs_tol_q,
                       max_num_steps,
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
  int rank_log_kp[system_size];
  int rank_log_ks[system_size];
  int rank_log_e50[system_size];
  int rank_log_k12[system_size];
  int rank_log_k21[system_size];
  int rank_log_sigma_y[system_size];
  int rank_log_a0[system_size];
  vector[system_size] bias_log_kp = log_kp - log_kp_;
  vector[system_size] bias_log_ks = log_ks - log_ks_;
  vector[system_size] bias_log_e50 = log_e50 - log_e50_;
  vector[system_size] bias_log_k12 = log_k12 - log_k12_;
  vector[system_size] bias_log_k21 = log_k21 - log_k21_;
  vector[system_size] bias_log_sigma_y = log_sigma_y - log_sigma_y_;
  vector[system_size] bias_log_a0 = log_a0 - log_a0_;

  for(i in 1:system_size) {
    rank_log_kp[i] = (log_kp[i] > log_kp_[i] ? 1 : 0);
    rank_log_ks[i] = (log_ks[i] > log_ks_[i] ? 1 : 0);
    rank_log_e50[i] = (log_e50[i] > log_e50_[i] ? 1 : 0);
    rank_log_k12[i] = (log_k12[i] > log_k12_[i] ? 1 : 0);
    rank_log_k21[i] = (log_k21[i] > log_k21_[i] ? 1 : 0);
    rank_log_sigma_y[i] = (log_sigma_y[i] > log_sigma_y_[i] ? 1 : 0);
    rank_log_a0[i] = (log_a0[i] > log_a0_[i] ? 1 : 0);
  }
  
}
