#include <test/unit/math/prim/prob/hmm_marginal_test.cpp>
#include <stan/math/prim/prob/hmm_latent_marginal_prob.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>


TEST_F(hmm_marginal_lpdf_test, rng) {
  using stan::math::hmm_latent_marginal_prob;

  Eigen::MatrixXd prob = hmm_latent_marginal_prob(log_omegas, Gamma, rho);
  std::cout << prob << std::endl;
}
