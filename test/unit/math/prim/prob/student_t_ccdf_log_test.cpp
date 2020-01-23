#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbStudentT, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double nu = 4.8;
  double mu = 2;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::student_t_lccdf(y, nu, mu, sigma)),
                  (stan::math::student_t_ccdf_log(y, nu, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::student_t_lccdf<double, double, double, double>(y, nu, mu,
                                                                   sigma)),
      (stan::math::student_t_ccdf_log<double, double, double, double>(y, nu, mu,
                                                                      sigma)));
}
