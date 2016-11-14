#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbPareto, cdf_log_matches_lcdf) {
  double y = 2.0;
  double y_min = 1.1;
  double alpha = 2.3;

  EXPECT_FLOAT_EQ((stan::math::pareto_lcdf(y, y_min, alpha)),
                  (stan::math::pareto_cdf_log(y, y_min, alpha)));
  EXPECT_FLOAT_EQ((stan::math::pareto_lcdf<double, double, double>(y, y_min, alpha)),
                  (stan::math::pareto_cdf_log<double, double, double>(y, y_min, alpha)));
}
