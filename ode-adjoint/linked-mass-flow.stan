functions {

  vector linked_mass_flow(real t, vector y,
                          real[] kt, real[] e50, real[] k12, real[] k21) {

    int num_main_states = num_elements(kt);
    int num_states = 2 * num_main_states;

    vector[num_states] dydt;
    vector[num_main_states] ksat;

    for (i in 1:num_main_states) {
      int m = 2 * (i-1) + 1;  // main state
      int a = m + 1;  // auxilary state
      ksat[i] = kt[i] * y[m] / (y[m] + e50[i]);

      dydt[m] = -ksat[i] * y[m] - k12[i] * y[m] + k21[i] * y[a];
      dydt[a] = +k12[i] * y[m] - k21[i] * y[a];

      if (i != 1) {
        dydt[m] += ksat[i-1] * y[2 * (i - 1) + 1];
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
                         vector log_kt, vector log_e50, vector log_k12, vector log_k21) {
    int num_sim = num_elements(ts);
    int num_states = num_elements(log_a0);
    int system_size = num_elements(log_kt);
    vector[num_states] y[num_sim];
    real kt[system_size];
    real e50[system_size];
    real k12[system_size];
    real k21[system_size];
    vector[num_states] a0;

    for(i in 1:system_size) {
      kt[i] = exp(log_kt[i]);
      e50[i] = exp(log_e50[i]);
      k12[i] = exp(log_k12[i]);
      k21[i] = exp(log_k21[i]);
    }

    for(i in 1:num_states) {
        a0[i] = exp(log_a0[i]);
    }

    if(adjoint_integrator) {
      y = ode_adams_tol(linked_mass_flow, a0, t0, ts, rel_tol, abs_tol, max_num_steps,
                        kt, e50, k12, k21);
    } else {
      y = ode_bdf_tol(linked_mass_flow, a0, t0, ts, rel_tol, abs_tol, max_num_steps,
                      kt, e50, k12, k21);
    }

    return y;
  }
  
}
data {
  real<lower=0> rel_tol;
  real<lower=0> abs_tol;
  int<lower=0,upper=1> adjoint_integrator;
  int<lower=1> max_num_steps;
  int<lower=1> system_size;
  int<lower=1> num_obs;
  real<lower=0> sigma_sim;
  real<lower=0> sigma_y;
}
transformed data {
  int num_states = 2 * system_size;
  real param_scale = sigma_sim / sqrt(num_obs/4);

  vector[system_size] log_kt_ = sample_vector_rng(0.0, sigma_sim, system_size);
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

  y_ = simulate_mean(t0, log_a0_, ts, 0, rel_tol, abs_tol, max_num_steps,
                     log_kt_, log_e50_, log_k12_, log_k21_);
  
  for(i in 1:num_obs) {
    for(j in 1:system_size) {
      int m = 2*(j-1) + 1;
      y_[i,m] = lognormal_rng(log(y_[i,m]+1E-3), exp(log_sigma_y_[j]));
    }
  }

  if(adjoint_integrator) {
    print("Using bdf_adjoint integrator.");
  } else {
    print("Using bdf integrator.");
  }
  print("relative tolerance: ", rel_tol);
  print("absolute tolerance: ", abs_tol);
  print("maximum number of steps: ", max_num_steps);
  print("number of time points: ", num_obs);
  print("system size: ", system_size);
  print("time points: ", ts);
  print("y_: ", y_);
}
parameters {
  vector<multiplier=param_scale,offset=log_kt_>[system_size] log_kt;
  vector<multiplier=param_scale,offset=log_e50_>[system_size] log_e50;
  vector<multiplier=param_scale,offset=log_k12_>[system_size] log_k12;
  vector<multiplier=param_scale,offset=log_k21_>[system_size] log_k21;
  vector<offset=log_sigma_y_>[system_size] log_sigma_y;
  vector<multiplier=param_scale,offset=log_a0_>[num_states] log_a0;
}
transformed parameters {
}
model {
  vector[num_states] mu[num_obs] = simulate_mean(t0, log_a0, ts, adjoint_integrator, rel_tol, abs_tol, max_num_steps,
                                                 log_kt, log_e50, log_k12, log_k21);
  
  target += normal_lpdf(log_kt| 0.0, sigma_sim);
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
  int rank_log_kt[system_size];
  int rank_log_e50[system_size];
  int rank_log_k12[system_size];
  int rank_log_k21[system_size];
  int rank_log_sigma_y[system_size];
  int rank_log_a0[system_size];
  vector[system_size] bias_log_kt = log_kt - log_kt_;
  vector[system_size] bias_log_e50 = log_e50 - log_e50_;
  vector[system_size] bias_log_k12 = log_k12 - log_k12_;
  vector[system_size] bias_log_k21 = log_k21 - log_k21_;
  vector[system_size] bias_log_sigma_y = log_sigma_y - log_sigma_y_;
  vector[system_size] bias_log_a0 = log_a0 - log_a0_;

  for(i in 1:system_size) {
    rank_log_kt[i] = (log_kt[i] > log_kt_[i] ? 1 : 0);
    rank_log_e50[i] = (log_e50[i] > log_e50_[i] ? 1 : 0);
    rank_log_k12[i] = (log_k12[i] > log_k12_[i] ? 1 : 0);
    rank_log_k21[i] = (log_k21[i] > log_k21_[i] ? 1 : 0);
    rank_log_sigma_y[i] = (log_sigma_y[i] > log_sigma_y_[i] ? 1 : 0);
    rank_log_a0[i] = (log_a0[i] > log_a0_[i] ? 1 : 0);
  }
  
}
