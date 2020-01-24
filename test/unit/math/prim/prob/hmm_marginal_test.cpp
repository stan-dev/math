// #include <stan/math/prim.hpp>
// #include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
// #include <test/unit/math/prim/prob/util.hpp>
// #include <boost/math/distributions.hpp>
// #include <boost/random/mersenne_twister.hpp>
// #include <limits>
// #include <vector>

#include <stan/math/prim/prob/hmm_marginal_lpdf.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <gtest/gtest.h>

// In the proposed example, the latent sate x determines
// the observational distribution:
//  (i) normal(mu, sigma)
//  (ii) normal(-mu, sigma)

double state_lpdf(double y, double abs_mu, double sigma, int state) {
  int x  = (-2 * state + 1);
  double chi =  (y + x * abs_mu) / sigma;
  return - 0.5 * chi * chi
         - 0.5 * std::log(283185307179586)
         - std::log(sigma);
}

// simulate an observation (emission)
double state_simu(double normal_variate, double abs_mu, double sigma,
                  int state) {
  int x = (-2 * state + 1);
  return x * abs_mu + sigma * normal_variate;
}

TEST(hmm_marginal_lpdf, two_state) {
  using stan::math::hmm_marginal_lpdf;
  using stan::math::var;

  int n_states = 2,
      n_transitions = 10;
  double abs_mu = 1,
         sigma = 1,
         p1_init = 0.65,
         gamma1 = 0.7,
         gamma2 = 0.45;

  // Simulate data
  // CHECK -- make sure we recover the same results.
  // If the rng is machine dependent, sub in with fixed data.
  Eigen::VectorXd rho(n_states);
  rho << p1_init, 1 - p1_init;

  Eigen::MatrixXd Gamma(n_states, n_states);
  Gamma << gamma1, gamma2, 1 - gamma1, 1 - gamma2;

  Eigen::VectorXd obs_data(n_transitions + 1);

  boost::random::mt19937 prng(1954);
  boost::random::discrete_distribution<> cat_init{rho[0], rho[1]};
  boost::random::discrete_distribution<>
    cat_zero{Gamma.col(0)[0], Gamma.col(0)[1]};
  boost::random::discrete_distribution<>
    cat_one{Gamma.col(1)[0], Gamma.col(1)[1]};

  boost::random::normal_distribution<> unit_normal(0, 1);

  int state = cat_init(prng);
  obs_data[0] = state_simu(unit_normal(prng), abs_mu, sigma, state);

  for (int n = 0; n < n_transitions; n++) {
    if (state == 0) state = cat_zero(prng);
    else state = cat_one(prng);

    obs_data[n + 1] = state_simu(unit_normal(prng), abs_mu, sigma, state);
  }

  // Compute observational densities
  Eigen::MatrixXd log_omegas(n_states, n_transitions + 1);
  for (int n = 0; n < n_transitions + 1; n++) {
    log_omegas.col(n)[0] = state_lpdf(obs_data[n], abs_mu, sigma, 0);
    log_omegas.col(n)[1] = state_lpdf(obs_data[n], abs_mu, sigma, 1);
  }

  // Compute marginal density
  Eigen::MatrixXd alpha_save(n_states, n_transitions + 1);
  Eigen::VectorXd alpha_log_norms_save(n_transitions + 1);
  Eigen::MatrixXd omegas_save;

  double density = hmm_marginal_lpdf(log_omegas, Gamma, rho,
                                     alpha_save, alpha_log_norms_save,
                                     omegas_save);

  // CHECK this density -- e.g. compare to Michael's code.
  std::cout << density << std::endl;  // 1.42004

  // Run the code with partials and operands.
  density = hmm_marginal_lpdf(log_omegas, Gamma, rho);

  std::cout << density << std::endl;

  // TEST GRADIENT EVALUATIONS
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
    log_omegas_v = log_omegas;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> Gamma_v = Gamma;
  Eigen::Matrix<var, Eigen::Dynamic, 1> rho_v = rho;

  var density_v = hmm_marginal_lpdf(log_omegas_v, Gamma_v, rho_v);
}
