#include <stan/math/prim/err/hmm_check.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(err, hmm_check) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::hmm_check;

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
      hmm_check(log_omegas, Gamma_rec, rho, "hmm_marginal_lpdf"),
      std::invalid_argument,
      "hmm_marginal_lpdf: Expecting a square matrix; rows of Gamma (2) "
      "and columns of Gamma (3) must match in size")

  // Gamma has a column that is not a simplex.
  MatrixXd Gamma_bad = Gamma;
  Gamma_bad(0, 0) = Gamma(0, 0) + 1;
  EXPECT_THROW_MSG(hmm_check(log_omegas, Gamma_bad, rho, "hmm_marginal_lpdf"),
                   std::domain_error,
                   "hmm_marginal_lpdf: Gamma[i, ] is not a valid simplex. "
                   "sum(Gamma[i, ]) = 2, but should be 1")

  // The size of Gamma is 0, even though there is at least one transition
  MatrixXd Gamma_empty(0, 0);
  EXPECT_THROW_MSG(
      hmm_check(log_omegas, Gamma_empty, rho, "hmm_marginal_lpdf"),
      std::invalid_argument,
      "hmm_marginal_lpdf: Gamma has size 0, but must have a non-zero size")

  // The size of Gamma is inconsistent with that of log_omega
  MatrixXd Gamma_wrong_size(n_states + 1, n_states + 1);

  EXPECT_THROW_MSG(
      hmm_check(log_omegas, Gamma_wrong_size, rho, "hmm_marginal_lpdf"),
      std::invalid_argument,
      "hmm_marginal_lpdf: Columns of Gamma (3)"
      " and Rows of log_omegas (2) must match in size")

  // rho is not a simplex.
  VectorXd rho_bad = rho;
  rho_bad(0) = rho(0) + 1;
  EXPECT_THROW_MSG(hmm_check(log_omegas, Gamma, rho_bad, "hmm_marginal_lpdf"),
                   std::domain_error,
                   "hmm_marginal_lpdf: rho is not a valid simplex. "
                   "sum(rho) = 2, but should be 1")

  // The size of rho is inconsistent with that of log_omega
  VectorXd rho_wrong_size(n_states + 1);
  EXPECT_THROW_MSG(
      hmm_check(log_omegas, Gamma, rho_wrong_size, "hmm_marginal_lpdf"),
      std::invalid_argument,
      "hmm_marginal_lpdf: rho has dimension = 3, expecting dimension = 2;"
      " a function was called with arguments of different scalar,"
      " array, vector, or matrix types, and they were not consistently sized;"
      "  all arguments must be scalars or multidimensional values of"
      " the same shape.")
}
