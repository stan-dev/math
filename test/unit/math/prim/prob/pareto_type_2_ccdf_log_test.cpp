#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbParetoType2, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double mu = -2;
  double lambda = 2;
  double alpha = 2.3;

  EXPECT_FLOAT_EQ((stan::math::pareto_type_2_lccdf(y, mu, lambda, alpha)),
                  (stan::math::pareto_type_2_ccdf_log(y, mu, lambda, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::pareto_type_2_lccdf<double, double, double, double>(
          y, mu, lambda, alpha)),
      (stan::math::pareto_type_2_ccdf_log<double, double, double, double>(
          y, mu, lambda, alpha)));
}
