#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbStudentT, cdf_log_matches_lcdf) {
  double y = 0.8;
  double nu = 4.8;
  double mu = 2;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::student_t_lcdf(y, nu, mu, sigma)),
                  (stan::math::student_t_cdf_log(y, nu, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::student_t_lcdf<double, double, double, double>(y, nu, mu,
                                                                  sigma)),
      (stan::math::student_t_cdf_log<double, double, double, double>(y, nu, mu,
                                                                     sigma)));
}
