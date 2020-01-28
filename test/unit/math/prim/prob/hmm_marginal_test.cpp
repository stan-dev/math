#include <stan/math/prim/prob/hmm_marginal_lpdf.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <test/unit/math/test_ad.hpp>
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
  // NOTE: This first test is exploratory, and will eventually be deleted.
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
  // std::cout << density << std::endl;  // 1.42004

  // Run the code with partials and operands.
  density = hmm_marginal_lpdf(log_omegas, Gamma, rho);

  // std::cout << density << std::endl;

  // TEST GRADIENT EVALUATIONS
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
    log_omegas_v = log_omegas;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> Gamma_v = Gamma;
  Eigen::Matrix<var, Eigen::Dynamic, 1> rho_v = rho;

  var density_v = hmm_marginal_lpdf(log_omegas_v, Gamma_v, rho_v);

  // Preliminary finite diff computation
  // REMARK: not quite correct, since we constrain rho and so forth
  // to be simplices.
  // double diff = 1e-6;
  // Eigen::MatrixXd Gamma_l = Gamma, Gamma_u = Gamma;
  // Gamma_l(0, 0) = Gamma(0, 0) - diff;
  // Gamma_u(0, 0) = Gamma(0, 0) + diff;
  // double gamma_diff
  //   = (hmm_marginal_lpdf(log_omegas, Gamma_u, rho)
  //     - hmm_marginal_lpdf(log_omegas, Gamma_l, rho)) / (2 * diff);

  // std::cout << "gamma_diff: " << gamma_diff << std::endl;

  // Eigen::VectorXd rho_l = rho, rho_u = rho;
  // rho_l(0) = rho(0) - diff;
  // rho_u(0) = rho(0) + diff;
  // double rho_diff
  // = (hmm_marginal_lpdf(log_omegas, Gamma, rho_u)
  //   - hmm_marginal_lpdf(log_omegas, Gamma, rho_l)) / (2 * diff);

  // std::cout << "rho diff: " << rho_diff << std::endl;

  // Eigen::MatrixXd log_omegas_u = log_omegas, log_omegas_l = log_omegas;
  // log_omegas_l(0, 0) = log_omegas(0, 0) - diff;
  // log_omegas_u(0, 0) = log_omegas(0, 0) + diff;
  // double omegas_diff
  //   = (hmm_marginal_lpdf(log_omegas_u, Gamma, rho)
  //     - hmm_marginal_lpdf(log_omegas_l, Gamma, rho)) / (2 * diff);

  // std::cout << "omegas_diff: " << omegas_diff << std::endl;
}

/**
 * Wrapper around hmm_marginal_density which passes rho and
 * Gamma without the last element of each column. We recover
 * the last element using the fact each column sums to 1.
 * The purpose of this function is to do finite diff benchmarking,
 * without breaking the simplex constraint.
 */
template <typename T_omega, typename T_Gamma, typename T_rho>
inline stan::return_type_t<T_omega, T_Gamma, T_rho>
hmm_marginal_test_wrapper (
  const Eigen::Matrix<T_omega, Eigen::Dynamic, Eigen::Dynamic>& log_omegas,
  const Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic>&
    Gamma_unconstrained,
  const std::vector<T_rho>& rho_unconstrained) {
  using stan::math::sum;
  using stan::math::col;
  int n_states = log_omegas.rows();

  Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic>
    Gamma(n_states, n_states);
  for (int j = 0; j < n_states; j++) {
    Gamma(n_states - 1, j) = 1 - sum(col(Gamma_unconstrained, j + 1));
    for (int i = 0; i < n_states - 1; i++) {
      Gamma(i, j) = Gamma_unconstrained(i, j);
    }
  }
  Eigen::Matrix<T_rho, Eigen::Dynamic, 1> rho(n_states);
  rho(1) = 1 - sum(rho_unconstrained);
  for (int i = 0; i < n_states - 1; i++) rho(i) = rho_unconstrained[i];

  return stan::math::hmm_marginal_lpdf(log_omegas, Gamma, rho);
 }


TEST(hmm_marginal_lpdf, autodiff) {
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

  // Construct "uncontrained" versions of rho and Gamma, without
  // the final element which can be determnied using the fact
  // the columns sum to 1. This allows us to do finite diff tests,
  // without violating the simplex constraint of rho and Gamma.
  std::vector<double> rho_unconstrained(n_states - 1);
  for (int i = 0; i < rho.size() - 1; i++)
    rho_unconstrained[i] = rho(i);

  Eigen::MatrixXd
    Gamma_unconstrained = Gamma.block(0, 0, n_states - 1, n_states);

  auto hmm_functor = [](const auto& log_omegas,
                        const auto& Gamma_unconstrained,
                        const auto& rho_unconstrained) {
    return hmm_marginal_test_wrapper(log_omegas, Gamma_unconstrained,
                                     rho_unconstrained);
  };

  stan::test::ad_tolerances tols;
  double infinity = std::numeric_limits<double>::infinity();
  tols.hessian_val_ = infinity;
  tols.hessian_grad_ = infinity;
  tols.hessian_hessian_ = infinity;
  tols.hessian_fvar_val_ = infinity;
  tols.hessian_fvar_grad_ = infinity;
  tols.hessian_fvar_hessian_ = infinity;
  tols.grad_hessian_val_ = infinity;
  tols.grad_hessian_hessian_ = infinity;
  tols.grad_hessian_grad_hessian_ = infinity;

  stan::test::expect_ad(tols, hmm_functor, log_omegas,
                        Gamma_unconstrained, rho_unconstrained);
}
