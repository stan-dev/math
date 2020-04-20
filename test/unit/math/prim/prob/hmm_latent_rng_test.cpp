#include <test/unit/math/prim/prob/hmm_marginal_test.cpp>
#include <stan/math/prim/prob/hmm_latent_rng.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>


TEST_F(hmm_marginal_lpdf_test, rng) {
  using stan::math::hmm_latent_rng;
  boost::random::mt19937 rng;

  std::vector<int> x = hmm_latent_rng(log_omegas, Gamma, rho, rng);
  for (size_t i = 0; i < x.size(); ++i)
    std::cout << x[i] << " ";
  std::cout << std::endl;
}

TEST_F(hmm_check_test, rng_exceptions) {
  using stan::math::hmm_latent_rng;

  boost::random::mt19937 rng;

  // Gamma is not square.
  EXPECT_THROW_MSG(
      hmm_latent_rng(log_omegas, Gamma_rec, rho, rng), std::invalid_argument,
      "hmm_latent_rng: Expecting a square matrix; rows of Gamma (2) "
      "and columns of Gamma (3) must match in size");

  // Gamma has a column that is not a simplex.
  EXPECT_THROW_MSG(hmm_latent_rng(log_omegas, Gamma_bad, rho, rng),
                   std::domain_error,
                   "hmm_latent_rng: Gamma[i, ] is not a valid simplex. "
                   "sum(Gamma[i, ]) = 2, but should be 1")

  // The size of Gamma is inconsistent with that of log_omega
  EXPECT_THROW_MSG(
      hmm_latent_rng(log_omegas, Gamma_wrong_size, rho, rng),
      std::invalid_argument,
      "hmm_latent_rng: Gamma has dimension = 3, expecting dimension = 2;"
      " a function was called with arguments of different scalar,"
      " array, vector, or matrix types, and they were not consistently sized;"
      "  all arguments must be scalars or multidimensional values of"
      " the same shape.")

  // rho is not a simplex.
  EXPECT_THROW_MSG(hmm_latent_rng(log_omegas, Gamma, rho_bad, rng),
                   std::domain_error,
                   "hmm_latent_rng: rho is not a valid simplex. "
                   "sum(rho) = 2, but should be 1")

  // The size of rho is inconsistent with that of log_omega
  EXPECT_THROW_MSG(
      hmm_latent_rng(log_omegas, Gamma, rho_wrong_size, rng),
      std::invalid_argument,
      "hmm_latent_rng: rho has dimension = 3, expecting dimension = 2;"
      " a function was called with arguments of different scalar,"
      " array, vector, or matrix types, and they were not consistently sized;"
      "  all arguments must be scalars or multidimensional values of"
      " the same shape.")
}
