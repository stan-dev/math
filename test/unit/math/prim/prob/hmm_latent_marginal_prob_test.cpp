#include <test/unit/math/prim/prob/hmm_marginal_test.cpp>
#include <stan/math/prim/prob/hmm_latent_marginal_prob.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>


TEST_F(hmm_test, latent_prob) {
  using stan::math::hmm_latent_marginal_prob;

  Eigen::MatrixXd prob = hmm_latent_marginal_prob(log_omegas_, Gamma_, rho_);
  std::cout << prob << std::endl;
}
/*
TEST_F(hmm_check_test, latent_prob_exceptions) {
  using stan::math::hmm_latent_marginal_prob;

  // Gamma is not square.
  EXPECT_THROW_MSG(
      hmm_latent_marginal_prob(log_omegas_, Gamma_rec_, rho_),
      std::invalid_argument,
      "hmm_latent_marginal_prob: Expecting a square matrix; rows of Gamma (2) "
      "and columns of Gamma (3) must match in size");

  // Gamma has a column that is not a simplex.
  EXPECT_THROW_MSG(hmm_latent_marginal_prob(log_omegas_, Gamma_bad_, rho_),
                   std::domain_error,
               "hmm_latent_marginal_prob: Gamma[i, ] is not a valid simplex. "
                   "sum(Gamma[i, ]) = 2, but should be 1")

  // The size of Gamma is inconsistent with that of log_omega
  EXPECT_THROW_MSG(
      hmm_latent_marginal_prob(log_omegas_, Gamma_wrong_size_, rho_),
      std::invalid_argument,
  "hmm_latent_marginal_prob: Gamma has dimension = 3, expecting dimension = 2;"
      " a function was called with arguments of different scalar,"
      " array, vector, or matrix types, and they were not consistently sized;"
      "  all arguments must be scalars or multidimensional values of"
      " the same shape.")

  // rho is not a simplex.
  EXPECT_THROW_MSG(hmm_latent_marginal_prob(log_omegas_, Gamma_, rho_bad_),
                   std::domain_error,
                   "hmm_latent_marginal_prob: rho is not a valid simplex. "
                   "sum(rho) = 2, but should be 1")

  // The size of rho is inconsistent with that of log_omega
  EXPECT_THROW_MSG(
      hmm_latent_marginal_prob(log_omegas_, Gamma_, rho_wrong_size_),
      std::invalid_argument,
    "hmm_latent_marginal_prob: rho has dimension = 3, expecting dimension = 2;"
      " a function was called with arguments of different scalar,"
      " array, vector, or matrix types, and they were not consistently sized;"
      "  all arguments must be scalars or multidimensional values of"
      " the same shape.")
}  */
