#include <stan/math/prim/prob/hmm_marginal_lpdf.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

/**
 * In the proposed example, the latent sate x determines
 * the observational distribution:
 *  (i) normal(mu, sigma)
 *  (ii) normal(-mu, sigma)
 */
double state_lpdf(double y, double abs_mu, double sigma, int state) {
  int x  = (-2 * state + 1);
  double chi =  (y - x * abs_mu) / sigma;
  return - 0.5 * chi * chi
         - 0.5 * std::log(6.283185307179586)
         - std::log(sigma);
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
  using stan::math::row;
  int n_states = log_omegas.rows();

  Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic>
    Gamma(n_states, n_states);
  for (int i = 0; i < n_states; i++) {
    Gamma(i, n_states - 1) = 1 - sum(row(Gamma_unconstrained, i + 1));
    for (int j = 0; j < n_states - 1; j++) {
      Gamma(i, j) = Gamma_unconstrained(i, j);
    }
  }

  Eigen::Matrix<T_rho, Eigen::Dynamic, 1> rho(n_states);
  rho(1) = 1 - sum(rho_unconstrained);
  for (int i = 0; i < n_states - 1; i++) rho(i) = rho_unconstrained[i];

  return stan::math::hmm_marginal_lpdf(log_omegas, Gamma, rho);
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
  Eigen::VectorXd rho(n_states);
  rho << p1_init, 1 - p1_init;

  Eigen::MatrixXd Gamma(n_states, n_states);
  Gamma << gamma1, 1 - gamma1, gamma2, 1 - gamma2;

  Eigen::VectorXd obs_data(n_transitions + 1);
  obs_data << -0.3315914, -0.1655340, -0.7984021, 0.2364608, -0.4489722,
    2.1831438, -1.4778675, 0.8717423, -1.0370874, 0.1370296, 1.9786208;

  // Compute observational densities
  Eigen::MatrixXd log_omegas(n_states, n_transitions + 1);
  for (int n = 0; n < n_transitions + 1; n++) {
    log_omegas.col(n)[0] = state_lpdf(obs_data[n], abs_mu, sigma, 0);
    log_omegas.col(n)[1] = state_lpdf(obs_data[n], abs_mu, sigma, 1);
  }

  // CHECK -- EXPECT_EQ returns an error.
  EXPECT_FLOAT_EQ(-18.37417, hmm_marginal_lpdf(log_omegas, Gamma, rho));

  // Construct "uncontrained" versions of rho and Gamma, without
  // the final element which can be determnied using the fact
  // the columns sum to 1. This allows us to do finite diff tests,
  // without violating the simplex constraint of rho and Gamma.
  std::vector<double> rho_unconstrained(n_states - 1);
  for (int i = 0; i < rho.size() - 1; i++)
    rho_unconstrained[i] = rho(i);

  Eigen::MatrixXd
    Gamma_unconstrained = Gamma.block(0, 0, n_states, n_states - 1);

  // Differentiation tests
  auto hmm_functor = [](const auto& log_omegas,
                        const auto& Gamma_unconstrained,
                        const auto& rho_unconstrained) {
    return hmm_marginal_test_wrapper(log_omegas, Gamma_unconstrained,
                                     rho_unconstrained);
  };

  double density = hmm_marginal_test_wrapper(log_omegas,
    Gamma_unconstrained, rho_unconstrained);

  stan::test::ad_tolerances tols;

  // Adjust tolerances to effectively not test hessians.
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

TEST(hmm_marginal_lpdf, exceptions) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::hmm_marginal_lpdf;

  int n_states = 2;
  int n_transitions = 2;
  MatrixXd log_omegas(n_states, n_transitions + 1);
  MatrixXd Gamma(n_states, n_states);
  VectorXd rho(n_states);

  for (int i = 0; i < n_states; i++)
    for (int j = 0; j < n_transitions + 1; j++)
      log_omegas(i, j) = 1;

  rho(0) = 0.65;
  rho(1) = 0.35;
  Gamma << 0.8, 0.2, 0.6, 0.4;

  // Gamma is not square.
  MatrixXd Gamma_rec(n_states, n_states + 1);
  EXPECT_THROW_MSG(
    hmm_marginal_lpdf(log_omegas, Gamma_rec, rho),
    std::invalid_argument,
    "hmm_marginal_lpdf: Expecting a square matrix; rows of Gamma (2) "
    "and columns of Gamma (3) must match in size");

  // Gamma has a column that is not a simplex.
  MatrixXd Gamma_bad = Gamma;
  Gamma_bad(0, 0) = Gamma(0, 0) + 1;
  EXPECT_THROW_MSG(
    hmm_marginal_lpdf(log_omegas, Gamma_bad, rho),
    std::domain_error,
    "hmm_marginal_lpdf: Gamma[i, ] is not a valid simplex. "
    "sum(Gamma[i, ]) = 2, but should be 1"
  )

  // The size of Gamma is inconsistent with that of log_omega
  MatrixXd Gamma_wrong_size(n_states + 1, n_states + 1);

  EXPECT_THROW_MSG(
    hmm_marginal_lpdf(log_omegas, Gamma_wrong_size, rho),
    std::invalid_argument,
    "hmm_marginal_lpdf: Gamma has dimension = 3, expecting dimension = 2;"
    " a function was called with arguments of different scalar,"
    " array, vector, or matrix types, and they were not consistently sized;"
    "  all arguments must be scalars or multidimensional values of"
    " the same shape."
  )

  // rho is not a simplex.
  VectorXd rho_bad = rho;
  rho_bad(0) = rho(0) + 1;
  EXPECT_THROW_MSG(
    hmm_marginal_lpdf(log_omegas, Gamma, rho_bad),
    std::domain_error,
    "hmm_marginal_lpdf: rho is not a valid simplex. "
    "sum(rho) = 2, but should be 1"
  )

  // The size of rho is inconsistent with that of log_omega
  VectorXd rho_wrong_size(n_states + 1);
  EXPECT_THROW_MSG(
    hmm_marginal_lpdf(log_omegas, Gamma, rho_wrong_size),
    std::invalid_argument,
    "hmm_marginal_lpdf: rho has dimension = 3, expecting dimension = 2;"
    " a function was called with arguments of different scalar,"
    " array, vector, or matrix types, and they were not consistently sized;"
    "  all arguments must be scalars or multidimensional values of"
    " the same shape."
  )
}
