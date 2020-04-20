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
    std::cout << x[i] << std::endl;

}
