#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbParetoType2, log_matches_lpdf) {
  double y = 0.8;
  double mu = -2;
  double lambda = 2;
  double alpha = 2.3;

  EXPECT_FLOAT_EQ((stan::math::pareto_type_2_lpdf(y, mu, lambda, alpha)),
                  (stan::math::pareto_type_2_log(y, mu, lambda, alpha)));
  EXPECT_FLOAT_EQ((stan::math::pareto_type_2_lpdf<true>(y, mu, lambda, alpha)),
                  (stan::math::pareto_type_2_log<true>(y, mu, lambda, alpha)));
  EXPECT_FLOAT_EQ((stan::math::pareto_type_2_lpdf<false>(y, mu, lambda, alpha)),
                  (stan::math::pareto_type_2_log<false>(y, mu, lambda, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::pareto_type_2_lpdf<true, double, double, double, double>(
          y, mu, lambda, alpha)),
      (stan::math::pareto_type_2_log<true, double, double, double, double>(
          y, mu, lambda, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::pareto_type_2_lpdf<false, double, double, double, double>(
          y, mu, lambda, alpha)),
      (stan::math::pareto_type_2_log<false, double, double, double, double>(
          y, mu, lambda, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::pareto_type_2_lpdf<double, double, double, double>(
          y, mu, lambda, alpha)),
      (stan::math::pareto_type_2_log<double, double, double, double>(
          y, mu, lambda, alpha)));
}
