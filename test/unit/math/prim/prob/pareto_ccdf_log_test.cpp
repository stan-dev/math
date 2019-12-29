#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbPareto, ccdf_log_matches_lccdf) {
  double y = 2.0;
  double y_min = 1.1;
  double alpha = 2.3;

  EXPECT_FLOAT_EQ((stan::math::pareto_lccdf(y, y_min, alpha)),
                  (stan::math::pareto_ccdf_log(y, y_min, alpha)));
  EXPECT_FLOAT_EQ(
      (stan::math::pareto_lccdf<double, double, double>(y, y_min, alpha)),
      (stan::math::pareto_ccdf_log<double, double, double>(y, y_min, alpha)));
}
