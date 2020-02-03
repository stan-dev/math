#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbDoubleExponential, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::double_exponential_lccdf(y, mu, sigma)),
                  (stan::math::double_exponential_ccdf_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::double_exponential_lccdf<double, double, double>(y, mu,
                                                                    sigma)),
      (stan::math::double_exponential_ccdf_log<double, double, double>(y, mu,
                                                                       sigma)));
}
