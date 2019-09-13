#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbPareto, log_matches_lpdf) {
  double y = 2.0;
  double y_min = 1.1;
  double alpha = 2.3;

  EXPECT_FLOAT_EQ((stan::math::pareto_lpdf(y, y_min, alpha)),
                  (stan::math::pareto_log(y, y_min, alpha)));
  EXPECT_FLOAT_EQ((stan::math::pareto_lpdf<true>(y, y_min, alpha)),
                  (stan::math::pareto_log<true>(y, y_min, alpha)));
  EXPECT_FLOAT_EQ((stan::math::pareto_lpdf<false>(y, y_min, alpha)),
                  (stan::math::pareto_log<false>(y, y_min, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::pareto_lpdf<true, double, double, double>(y, y_min, alpha)),
      (stan::math::pareto_log<true, double, double, double>(y, y_min, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::pareto_lpdf<false, double, double, double>(y, y_min, alpha)),
      (stan::math::pareto_log<false, double, double, double>(y, y_min, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::pareto_lpdf<double, double, double>(y, y_min, alpha)),
      (stan::math::pareto_log<double, double, double>(y, y_min, alpha)));
}
