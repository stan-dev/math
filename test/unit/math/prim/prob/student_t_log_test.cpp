#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbStudentT, log_matches_lpdf) {
  double y = 0.8;
  double nu = 4.8;
  double mu = 2;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::student_t_lpdf(y, nu, mu, sigma)),
                  (stan::math::student_t_log(y, nu, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::student_t_lpdf<true>(y, nu, mu, sigma)),
                  (stan::math::student_t_log<true>(y, nu, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::student_t_lpdf<false>(y, nu, mu, sigma)),
                  (stan::math::student_t_log<false>(y, nu, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::student_t_lpdf<true, double, double, double, double>(
          y, nu, mu, sigma)),
      (stan::math::student_t_log<true, double, double, double, double>(
          y, nu, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::student_t_lpdf<false, double, double, double, double>(
          y, nu, mu, sigma)),
      (stan::math::student_t_log<false, double, double, double, double>(
          y, nu, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::student_t_lpdf<double, double, double, double>(
                      y, nu, mu, sigma)),
                  (stan::math::student_t_log<double, double, double, double>(
                      y, nu, mu, sigma)));
}
