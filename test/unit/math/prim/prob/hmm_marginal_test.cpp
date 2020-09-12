#include <test/unit/math/prim/prob/hmm_util.hpp>
#include <stan/math/prim/prob/hmm_marginal.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

// For evaluation of the density, the C++ code is benchmarked against
// a forward algorithm written in R.
// TODO(charlesm93): Add public repo link with R script.
TEST_F(hmm_test, ten_transitions) {
  using stan::math::hmm_marginal;

  EXPECT_FLOAT_EQ(-18.37417, hmm_marginal(log_omegas_, Gamma_, rho_));

  // Differentiation tests
  auto hmm_functor = [](const auto& log_omegas, const auto& Gamma_unconstrained,
                        const auto& rho_unconstrained) {
    return hmm_marginal_test_wrapper(log_omegas, Gamma_unconstrained,
                                     rho_unconstrained);
  };

  stan::test::expect_ad(tols_, hmm_functor, log_omegas_, Gamma_unconstrained_,
                        rho_unconstrained_);
}

TEST_F(hmm_test, zero_transitions) {
  using stan::math::hmm_marginal;

  EXPECT_FLOAT_EQ(-1.520827, hmm_marginal(log_omegas_zero_, Gamma_, rho_));

  // Differentiation tests
  auto hmm_functor = [](const auto& log_omegas, const auto& Gamma_unconstrained,
                        const auto& rho_unconstrained) {
    return hmm_marginal_test_wrapper(log_omegas, Gamma_unconstrained,
                                     rho_unconstrained);
  };

  stan::test::expect_ad(tols_, hmm_functor, log_omegas_zero_,
                        Gamma_unconstrained_, rho_unconstrained_);
}

TEST(hmm_marginal, one_state) {
  using stan::math::hmm_marginal;
  int n_states = 1, p1_init = 1, gamma1 = 1, n_transitions = 10, abs_mu = 1,
      sigma = 1;
  Eigen::VectorXd rho(n_states);
  rho << p1_init;
  Eigen::MatrixXd Gamma(n_states, n_states);
  Gamma << gamma1;
  Eigen::VectorXd obs_data(n_transitions + 1);
  obs_data << -0.9692032, 1.6367754, 1.0339449, 0.9798393, 0.4829358, 2.7508704,
      0.3122448, 1.8316583, 1.6327319, 1.2097332, 0.4087620;
  Eigen::MatrixXd log_omegas(n_states, n_transitions + 1);
  for (int n = 0; n < n_transitions + 1; n++)
    log_omegas.col(n)[0] = state_lpdf(obs_data[n], abs_mu, sigma, 0);

  EXPECT_FLOAT_EQ(-14.89646, hmm_marginal(log_omegas, Gamma, rho));

  // Differentiation tests
  // In the case where we have one state, Gamma and rho
  // are fixed (i.e = 1)
  auto hmm_functor = [](const auto& log_omegas) {
    Eigen::MatrixXd Gamma(1, 1);
    Gamma << 1;
    Eigen::VectorXd rho(1);
    rho << 1;

    return hmm_marginal(log_omegas, Gamma, rho);
  };

  stan::test::ad_tolerances tols;
  stan::test::expect_ad(tols, hmm_functor, log_omegas);
}

TEST(hmm_marginal, exceptions) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::hmm_marginal;

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
  EXPECT_THROW_MSG(hmm_marginal(log_omegas, Gamma_rec, rho),
                   std::invalid_argument,
                   "hmm_marginal: Expecting a square matrix; rows of Gamma (2) "
                   "and columns of Gamma (3) must match in size");

  // Gamma has a column that is not a simplex.
  MatrixXd Gamma_bad = Gamma;
  Gamma_bad(0, 0) = Gamma(0, 0) + 1;
  EXPECT_THROW_MSG(hmm_marginal(log_omegas, Gamma_bad, rho), std::domain_error,
                   "hmm_marginal: Gamma[i, ] is not a valid simplex. "
                   "sum(Gamma[i, ]) = 2, but should be 1")

  // The size of Gamma is 0, even though there is at least one transition
  MatrixXd Gamma_empty(0, 0);
  EXPECT_THROW_MSG(
      hmm_marginal(log_omegas, Gamma_empty, rho), std::invalid_argument,
      "hmm_marginal: Gamma has size 0, but must have a non-zero size")

  // The size of Gamma is inconsistent with that of log_omega
  MatrixXd Gamma_wrong_size(n_states + 1, n_states + 1);

  EXPECT_THROW_MSG(hmm_marginal(log_omegas, Gamma_wrong_size, rho),
                   std::invalid_argument,
                   "hmm_marginal: Columns of Gamma (3)"
                   " and Rows of log_omegas (2) must match in size")

  // rho is not a simplex.
  VectorXd rho_bad = rho;
  rho_bad(0) = rho(0) + 1;
  EXPECT_THROW_MSG(hmm_marginal(log_omegas, Gamma, rho_bad), std::domain_error,
                   "hmm_marginal: rho is not a valid simplex. "
                   "sum(rho) = 2, but should be 1")

  // The size of rho is inconsistent with that of log_omega
  VectorXd rho_wrong_size(n_states + 1);
  EXPECT_THROW_MSG(
      hmm_marginal(log_omegas, Gamma, rho_wrong_size), std::invalid_argument,
      "hmm_marginal: rho has dimension = 3, expecting dimension = 2;"
      " a function was called with arguments of different scalar,"
      " array, vector, or matrix types, and they were not consistently sized;"
      "  all arguments must be scalars or multidimensional values of"
      " the same shape.")
}
