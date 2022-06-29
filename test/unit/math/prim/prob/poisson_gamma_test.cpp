#include <stan/math/prim.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ProbDistributionsPoissonGamma, values) {
  using stan::math::poisson_gamma_lpmf;
  int y = 10;
  double alpha = 4;
  double beta = 2;

  double lpmf = poisson_gamma_lpmf(y, alpha, beta);

  EXPECT_FLOAT_EQ(lpmf, -6.95199150829391);
}
