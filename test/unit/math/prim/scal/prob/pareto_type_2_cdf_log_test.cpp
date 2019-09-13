#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbParetoType2, cdf_log_matches_lcdf) {
  double y = 0.8;
  double mu = -2;
  double lambda = 2;
  double alpha = 2.3;

  EXPECT_FLOAT_EQ((stan::math::pareto_type_2_lcdf(y, mu, lambda, alpha)),
                  (stan::math::pareto_type_2_cdf_log(y, mu, lambda, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::pareto_type_2_lcdf<double, double, double, double>(
          y, mu, lambda, alpha)),
      (stan::math::pareto_type_2_cdf_log<double, double, double, double>(
          y, mu, lambda, alpha)));
}
